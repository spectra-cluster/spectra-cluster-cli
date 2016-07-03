package uk.ac.ebi.pride.spectracluster.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.FileFilterUtils;
import uk.ac.ebi.pride.spectracluster.binning.BinningSpectrumConverter;
import uk.ac.ebi.pride.spectracluster.cdf.CdfLearner;
import uk.ac.ebi.pride.spectracluster.cdf.CdfResult;
import uk.ac.ebi.pride.spectracluster.cdf.CumulativeDistributionFunction;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusterer;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusteringCallable;
import uk.ac.ebi.pride.spectracluster.conversion.MergingCGFConverter;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.implementation.SpectraClusterStandalone;
import uk.ac.ebi.pride.spectracluster.io.CGFSpectrumIterable;
import uk.ac.ebi.pride.spectracluster.io.DotClusterClusterAppender;
import uk.ac.ebi.pride.spectracluster.merging.BinaryFileMergingClusterer;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.util.*;
import uk.ac.ebi.pride.spectracluster.util.function.Functions;
import uk.ac.ebi.pride.spectracluster.util.function.peak.HighestNPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveReporterIonPeaksFunction;

import java.io.*;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: jg
 * Date: 9/15/13
 * Time: 11:19 AM
 * To change this template use File | Settings | File Templates.
 */
public class SpectraClusterCliMain implements IProgressListener {
    public static final boolean DELETE_TEMPORARY_CLUSTERING_RESULTS = true;

    private boolean verbose;

    public static void main(String[] args) {
        SpectraClusterCliMain instance = new SpectraClusterCliMain();
        instance.run(args);
    }

    private void run(String[] args) {
        CommandLineParser parser = new GnuParser();

        try {
            CommandLine commandLine = parser.parse(CliOptions.getOptions(), args);

            // HELP
            if (commandLine.hasOption(CliOptions.OPTIONS.HELP.getValue())) {
                printUsage();
                return;
            }

            // create the clustering standalone object
            SpectraClusterStandalone spectraClusterStandalone = new SpectraClusterStandalone();

            // RESULT FILE PATH
            if (!commandLine.hasOption(CliOptions.OPTIONS.OUTPUT_PATH.getValue()))
                throw new MissingParameterException("Missing required option " + CliOptions.OPTIONS.OUTPUT_PATH.getValue());
            File finalResultFile = new File(commandLine.getOptionValue(CliOptions.OPTIONS.OUTPUT_PATH.getValue()));

            if (finalResultFile.exists())
                throw new Exception("Result file " + finalResultFile + " already exists");

            // NUMBER OF JOBS
            int paralellJobs = 4;
            if (commandLine.hasOption(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue())) {
                paralellJobs = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue()));
                spectraClusterStandalone.setParallelJobs(paralellJobs);
            }

            // NUMBER OF ROUNDS
            int rounds = 4;
            if (commandLine.hasOption(CliOptions.OPTIONS.ROUNDS.getValue()))
                rounds = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ROUNDS.getValue()));

            // START THRESHOLD
            float startThreshold = 0.999F;
            if (commandLine.hasOption(CliOptions.OPTIONS.START_THRESHOLD.getValue()))
                startThreshold = Float.parseFloat(commandLine.getOptionValue(CliOptions.OPTIONS.START_THRESHOLD.getValue()));

            // END THRESHOLD
            float endThreshold = 0.99F;
            if (commandLine.hasOption(CliOptions.OPTIONS.END_THRESHOLD.getValue()))
                endThreshold = Float.parseFloat(commandLine.getOptionValue(CliOptions.OPTIONS.END_THRESHOLD.getValue()));

            List<Float> thresholds = spectraClusterStandalone.generateClusteringThresholds(startThreshold, endThreshold, rounds);

            // PRECURSOR TOLERANCE
            if (commandLine.hasOption(CliOptions.OPTIONS.PRECURSOR_TOLERANCE.getValue())) {
                float precursorTolerance = Float.parseFloat(commandLine.getOptionValue(CliOptions.OPTIONS.PRECURSOR_TOLERANCE.getValue()));
                Defaults.setDefaultPrecursorIonTolerance(precursorTolerance);
            }

            // FRAGMENT ION TOLERANCE
            if (commandLine.hasOption(CliOptions.OPTIONS.FRAGMENT_TOLERANCE.getValue())) {
                float fragmentTolerance = Float.parseFloat(commandLine.getOptionValue(CliOptions.OPTIONS.FRAGMENT_TOLERANCE.getValue()));
                Defaults.setFragmentIonTolerance(fragmentTolerance);
            }

            // BINARY TMP DIR
            if (commandLine.hasOption(CliOptions.OPTIONS.BINARY_TMP_DIR.getValue())) {
                File binaryTmpDirectory = new File(commandLine.getOptionValue(CliOptions.OPTIONS.BINARY_TMP_DIR.getValue()));
                spectraClusterStandalone.setTemporaryDirectory(binaryTmpDirectory);
            }

            // KEEP BINARY FILES
            if (commandLine.hasOption(CliOptions.OPTIONS.KEEP_BINARY_FILE.getValue())) {
                spectraClusterStandalone.setKeepBinaryFiles(true);
            }

            // FAST MODE
            if (commandLine.hasOption(CliOptions.OPTIONS.FAST_MODE.getValue())) {
                spectraClusterStandalone.setUseFastMode(true);
            }

            // VERBOSE
            if (commandLine.hasOption(CliOptions.OPTIONS.VERBOSE.getValue())) {
                spectraClusterStandalone.setVerbose(true);
            }

            // REMOVE QUANT PEAKS
            if (commandLine.hasOption(CliOptions.OPTIONS.REMOVE_REPORTER_PEAKS.getValue())) {
                String quantTypeString = commandLine.getOptionValue(CliOptions.OPTIONS.REMOVE_REPORTER_PEAKS.getValue());
                RemoveReporterIonPeaksFunction.REPORTER_TYPE reporterType;

                if (quantTypeString.toLowerCase().equals("itraq")) {
                    reporterType = RemoveReporterIonPeaksFunction.REPORTER_TYPE.ITRAQ;
                }
                else if (quantTypeString.toLowerCase().equals("tmt")) {
                    reporterType = RemoveReporterIonPeaksFunction.REPORTER_TYPE.TMT;
                }
                else if (quantTypeString.toLowerCase().equals("all")) {
                    reporterType = RemoveReporterIonPeaksFunction.REPORTER_TYPE.ALL;
                }
                else {
                    throw new MissingParameterException("Invalid reporter type defined. Valid values are " +
                            "'ITRAQ', 'TMT', and 'ALL'.");
                }

                spectraClusterStandalone.setReporterType(reporterType);
            }

            // FILES TO PROCESS
            String[] peaklistFilenames = commandLine.getArgs();

            // RE-USE BINARY FILES
            boolean reUseBinaryFiles = commandLine.hasOption(CliOptions.OPTIONS.REUSE_BINARY_FILES.getValue());

            // if re-use is set, binaryTmpDirectory is required and merging is impossible
            if (reUseBinaryFiles && !commandLine.hasOption(CliOptions.OPTIONS.BINARY_TMP_DIR.getValue()))
                throw new MissingParameterException("Missing required option '" + CliOptions.OPTIONS.BINARY_TMP_DIR.getValue() + "' with " + CliOptions.OPTIONS.REUSE_BINARY_FILES.getValue());

            if (reUseBinaryFiles && peaklistFilenames.length > 0)
                System.out.println("WARNING: " + CliOptions.OPTIONS.REUSE_BINARY_FILES.getValue() + " set, input files will be ignored");

            // make sure input files were set
            if (!reUseBinaryFiles && peaklistFilenames.length < 1)
                throw new MissingParameterException("No spectrum files passed. Please list the peak list files to process after the command.");

            /**
             * Advanced options
             */
            // MIN NUMBER COMPARISONS
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_MIN_NUMBER_COMPARISONS.getValue())) {
                int minComparisons = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_MIN_NUMBER_COMPARISONS.getValue()));
                Defaults.setMinNumberComparisons(minComparisons);
            }

            // N HIGHEST PEAKS
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_NUMBER_PREFILTERED_PEAKS.getValue())) {
                int nHighestPeaks = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_NUMBER_PREFILTERED_PEAKS.getValue()));
                ClusteringSettings.setLoadingSpectrumFilter(new HighestNPeakFunction(nHighestPeaks));
            }

            /**
             * ------ Learn the CDF if set --------
             */
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_LEARN_CDF.getValue())) {
                String cdfOuputFilename = commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_LEARN_CDF.getValue());
                File cdfOutputFile = new File(cdfOuputFilename);

                if (cdfOutputFile.exists()) {
                    throw new Exception("CDF output file " + cdfOuputFilename + " already exists.");
                }

                CdfLearner cdfLearner = new CdfLearner();
                System.out.println("Learning CDF...");
                CdfResult cdfResult = cdfLearner.learnCumulativeDistribution(peaklistFilenames, paralellJobs);

                // write it to the file
                FileWriter writer = new FileWriter(cdfOutputFile);
                writer.write(cdfResult.toString());
                writer.close();

                System.out.println("CDF successfully written to " + cdfOuputFilename);
                return;
            }

            /**
             * ------ Load the CDF from file -------
             */
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_LOAD_CDF_FILE.getValue())) {
                BufferedReader reader = new BufferedReader(
                        new FileReader(
                                commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_LOAD_CDF_FILE.getValue())));

                StringBuilder cdfString = new StringBuilder();
                String line;
                while ((line = reader.readLine()) != null) {
                    cdfString.append(line);
                }
                reader.close();

                Defaults.setCumulativeDistributionFunction(CumulativeDistributionFunction.fromString(cdfString.toString()));
            }

            /**
             * ------- THE ACTUAL LOGIC STARTS HERE -----------
             */
            printSettings(finalResultFile, paralellJobs, startThreshold, endThreshold, rounds, spectraClusterStandalone.isKeepBinaryFiles(),
                    spectraClusterStandalone.getTemporaryDirectory(), peaklistFilenames, reUseBinaryFiles, spectraClusterStandalone.isUseFastMode());

            spectraClusterStandalone.addProgressListener(this);

            if (reUseBinaryFiles) {
                spectraClusterStandalone.clusterExistingBinaryFiles(
                        spectraClusterStandalone.getTemporaryDirectory(), thresholds, finalResultFile);
            }
            else {
                List<File> peaklistFiles = new ArrayList<File>(peaklistFilenames.length);
                for (String filename : peaklistFilenames) {
                    peaklistFiles.add(new File(filename));
                }

                spectraClusterStandalone.clusterPeaklistFiles(peaklistFiles, thresholds, finalResultFile);
            }
        } catch (MissingParameterException e) {
            System.out.println("Error: " + e.getMessage() + "\n\n");
            printUsage();

            System.exit(1);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Error: " + e.getMessage());

            System.exit(1);
        }
    }

    private void printSettings(File finalResultFile, int nMajorPeakJobs, float startThreshold,
                                      float endThreshold, int rounds, boolean keepBinaryFiles, File binaryTmpDirectory,
                                      String[] peaklistFilenames, boolean reUseBinaryFiles, boolean fastMode) {
        System.out.println("spectra-cluster API Version 1.0");
        System.out.println("Created by Rui Wang & Johannes Griss\n");

        System.out.println("-- Settings --");
        System.out.println("Number of threads: " + String.valueOf(nMajorPeakJobs));
        System.out.println("Thresholds: " + String.valueOf(startThreshold) + " - " + String.valueOf(endThreshold) + " in " + rounds + " rounds");
        System.out.println("Keeping binary files: " + (keepBinaryFiles ? "true" : "false"));
        System.out.println("Binary file directory: " + binaryTmpDirectory);
        System.out.println("Result file: " + finalResultFile);
        System.out.println("Reuse binary files: " + (reUseBinaryFiles ? "true" : "false"));
        System.out.println("Input files: " + peaklistFilenames.length);
        System.out.println("Using fast mode: " + (fastMode ? "yes" : "no"));

        System.out.println("\nOther settings:");
        System.out.println("Precursor tolerance: " + Defaults.getDefaultPrecursorIonTolerance());
        System.out.println("Fragment ion tolerance: " + Defaults.getFragmentIonTolerance());

        // only show certain settings if they were changed
        if (Defaults.getMinNumberComparisons() != Defaults.DEFAULT_MIN_NUMBER_COMPARISONS)
            System.out.println("Minimum number of comparisons: " + Defaults.getMinNumberComparisons());

        System.out.println();
    }

    private void printUsage() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Spectra Cluster - Clusterer",
                "Clusters the spectra found in an MGF file and writes the results in a text-based file.\n",
                CliOptions.getOptions(), "\n\n", true);
    }

    @Override
    public void onProgressUpdate(ProgressUpdate progressUpdate) {
        System.out.println(progressUpdate.getMessage());
    }
}
