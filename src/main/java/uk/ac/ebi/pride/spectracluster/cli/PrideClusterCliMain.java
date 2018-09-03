package uk.ac.ebi.pride.spectracluster.cli;

import org.apache.commons.cli.*;
import uk.ac.ebi.pride.spectracluster.binning.BinningSpectrumConverter;
import uk.ac.ebi.pride.spectracluster.binning.FixedReferenceMzBinner;
import uk.ac.ebi.pride.spectracluster.cdf.*;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusterer;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.implementation.ScoreCalculator;
import uk.ac.ebi.pride.spectracluster.implementation.SpectraClusterStandalone;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterIterable;
import uk.ac.ebi.pride.spectracluster.io.DotClusterClusterAppender;
import uk.ac.ebi.pride.spectracluster.merging.BinaryFileMergingClusterer;
import uk.ac.ebi.pride.spectracluster.util.*;
import uk.ac.ebi.pride.spectracluster.util.function.peak.BinnedHighestNPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.HighestNPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveReporterIonPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.cluster.ClusterOnlyIdentifiedPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.cluster.ClusterOnlyUnidentifiedPredicate;

import java.io.*;
import java.security.InvalidParameterException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created with IntelliJ IDEA.
 * User: jg
 * Date: 9/15/13
 * Time: 11:19 AM
 * To change this template use File | Settings | File Templates.
 */
public class PrideClusterCliMain implements IProgressListener {
    public enum OPTIONS {
        WRITE_BINARY_FILE("write_binary"),
        MERGE_BINARY_FILE("merge_binary"),
        CLUSTER_BINARY_FILE("cluster_binary"),
        CONVERT_BINARY_FILE("convert_binary"),
        CLUSTER_MERGE_BINARY_FILES("cluster_merge_binary");

        private String value;

        OPTIONS(String value) {
            this.value = value;
        }

        public String getValue() {
            return value;
        }

        @Override
        public String toString() {
            return value;
        }
    }

    public static final boolean DELETE_TEMPORARY_CLUSTERING_RESULTS = true;

    private boolean verbose;

    public static void main(String[] args) {
        PrideClusterCliMain instance = new PrideClusterCliMain();
        instance.run(args);
    }

    private void run(String[] args) {
        CommandLineParser parser = new GnuParser();

        try {
            // Add custom PRIDE Cluster command line options
            Options orgOptions = CliOptions.getOptions();

            orgOptions.addOption(OptionBuilder.withDescription("Create binary files").create(OPTIONS.WRITE_BINARY_FILE.getValue()));
            orgOptions.addOption(OptionBuilder
                    .withDescription("Merge binary files. Input files are treated as binary result files. Output path " +
                            "is used as the name of the merged file")
                    .create(OPTIONS.MERGE_BINARY_FILE.getValue()));
            orgOptions.addOption(OptionBuilder
                    .withDescription("Cluster binary files")
                    .create(OPTIONS.CLUSTER_BINARY_FILE.getValue()));
            orgOptions.addOption(OptionBuilder
                    .withDescription("Converts the passed binary files to the specified result file in the .clustering format")
                    .create(OPTIONS.CONVERT_BINARY_FILE.getValue()));
            orgOptions.addOption(OptionBuilder
                    .withDescription("Performs the clustering merging step on the passed binary files.")
                    .create(OPTIONS.CLUSTER_MERGE_BINARY_FILES.getValue()));

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

            // NUMBER OF JOBS
            int paralellJobs = 4;
            if (commandLine.hasOption(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue())) {
                paralellJobs = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue()));
                spectraClusterStandalone.setParallelJobs(paralellJobs);
            }

            // NUMBER OF ROUNDS
            int rounds = 5;
            if (commandLine.hasOption(CliOptions.OPTIONS.ROUNDS.getValue()))
                rounds = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ROUNDS.getValue()));

            // START THRESHOLD
            float startThreshold = 1F;
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

            // PRECURSOR TOLERANCE UNIT
            if (commandLine.hasOption(CliOptions.OPTIONS.PRECURSOR_TOELRANCE_UNIT.getValue())) {
                String unit = commandLine.getOptionValue(CliOptions.OPTIONS.PRECURSOR_TOELRANCE_UNIT.getValue()).toLowerCase();

                if ("ppm".equals(unit)) {
                    // adapt the precursor tolerance
                    float userTolerance = Defaults.getDefaultPrecursorIonTolerance();
                    float ppmFraction = userTolerance / 1000000;

                    // set the "precursor" tolerance based on a high m/z
                    Defaults.setDefaultPrecursorIonTolerance(ppmFraction * 3500);

                    // set the actual ppm tolerance
                    ClusteringSettings.ppmThreshold = userTolerance;
                }
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

            // LOADING MODE
            if (commandLine.hasOption(CliOptions.OPTIONS.ONLY_IDENTIFIED.getValue()) && commandLine.hasOption(CliOptions.OPTIONS.ONLY_UNIDENTIFIED.getValue())) {
                throw new Exception(CliOptions.OPTIONS.ONLY_IDENTIFIED.getValue() + " and " + CliOptions.OPTIONS.ONLY_UNIDENTIFIED.getValue() +
                        " cannot be used together");
            }
            if (commandLine.hasOption(CliOptions.OPTIONS.ONLY_UNIDENTIFIED.getValue())) {
                ClusteringSettings.setLoadingMode(ClusteringSettings.LOADING_MODE.ONLY_IDENTIFIED);
            }
            if (commandLine.hasOption(CliOptions.OPTIONS.ONLY_UNIDENTIFIED.getValue())) {
                ClusteringSettings.setLoadingMode(ClusteringSettings.LOADING_MODE.ONLY_UNIDENTIFIED);
            }

            // VERBOSE
            if (commandLine.hasOption(CliOptions.OPTIONS.VERBOSE.getValue())) {
                spectraClusterStandalone.setVerbose(true);
                Defaults.setSaveDebugInformation(true);
                Defaults.setSaveAddingScore(true);
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

            // add the filters
            List<String> addedFilters = new ArrayList<String>();
            if (commandLine.hasOption(CliOptions.OPTIONS.FILTER.getValue())) {
                for (String filterName : commandLine.getOptionValues(CliOptions.OPTIONS.FILTER.getValue())) {
                    ClusteringSettings.SPECTRUM_FILTER filter = ClusteringSettings.SPECTRUM_FILTER.getFilterForName(filterName);

                    if (filter == null) {
                        throw new InvalidParameterException("Error: Unknown filter name passed: '" + filterName + "'");
                    }

                    ClusteringSettings.addIntitalSpectrumFilter(filter.filter);
                    addedFilters.add(filterName);
                }
            }

            /**
             * Advanced options
             */
            // MIN NUMBER COMPARISONS
            // By default, use the SpectraPerBinNumberComparisonAssessor
            // if the command line option is set, use the MinNumberComparisonAssessor instead.
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_MIN_NUMBER_COMPARISONS.getValue())) {
                int minComparisons = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_MIN_NUMBER_COMPARISONS.getValue()));
                Defaults.setNumberOfComparisonAssessor(new MinNumberComparisonsAssessor(minComparisons));
            }
            // adaptive version with a min number of spectra set
            else if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_MIN_ADAPTIVE_COMPARISONS.getValue())) {
                int minComparisons = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_MIN_ADAPTIVE_COMPARISONS.getValue()));
                Defaults.setNumberOfComparisonAssessor(new SpectraPerBinNumberComparisonAssessor(Defaults.getDefaultPrecursorIonTolerance(),
                        minComparisons));
            } else {
                Defaults.setNumberOfComparisonAssessor(new SpectraPerBinNumberComparisonAssessor(Defaults.getDefaultPrecursorIonTolerance()));
            }

            // N HIGHEST PEAKS
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_NUMBER_PREFILTERED_PEAKS.getValue())) {
                int nHighestPeaks = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_NUMBER_PREFILTERED_PEAKS.getValue()));
                ClusteringSettings.setLoadingSpectrumFilter(new HighestNPeakFunction(nHighestPeaks));
            }

            // FILTER PEAKS PER MZ
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_FILTER_PEAKS_PER_MZ.getValue())) {
                ClusteringSettings.setLoadingSpectrumFilter(new BinnedHighestNPeakFunction(10, 100, 30));
            }

            // MIN CONSENSUS PEAKS TO KEEP
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_MIN_CONSENSUS_PEAKS_TO_KEEP.getValue())) {
                Defaults.setDefaultConsensusMinPeaks(Integer.parseInt(
                        commandLine.getOptionValue(CliOptions.OPTIONS.ADVANCED_MIN_CONSENSUS_PEAKS_TO_KEEP.getValue())));
            }

            // MGF COMMENT SUPPORT
            ClusteringSettings.disableMGFCommentSupport = commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_DISABLE_MGF_COMMENTS.getValue());

            // MERGE BINARY FILES
            boolean mergeBinaryFilesMode = commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_MERGE_BINARY_FILES.getValue());

            /**
             * ------ Convert CLS mode ----
             */
            if (commandLine.hasOption(CliOptions.OPTIONS.ADVANCED_CONVERT_CGF.getValue())) {
                convertToCgf(peaklistFilenames, finalResultFile);
                return;
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
             * ------- PRIDE Cluster specific methods ---------
             */
            if (commandLine.hasOption(OPTIONS.WRITE_BINARY_FILE.getValue())) {
                System.out.println("Creating binary files...");
                createBinaryFiles(peaklistFilenames, finalResultFile, paralellJobs);
                System.exit(0);
            }

            if (commandLine.hasOption(OPTIONS.MERGE_BINARY_FILE.getValue())) {
                System.out.println("Merging binary files...");
                mergeBinaryFiles(peaklistFilenames, finalResultFile);
                System.exit(0);
            }

            if (commandLine.hasOption(OPTIONS.CLUSTER_BINARY_FILE.getValue())) {
                System.out.println("Clustering files...");
                clusterBinaryFile(peaklistFilenames, finalResultFile, paralellJobs, thresholds, spectraClusterStandalone.isVerbose());
                System.exit(0);
            }

            if (commandLine.hasOption(OPTIONS.CONVERT_BINARY_FILE.getValue())) {
                System.out.println("Converting binary files...");
                convertBinaryFiles(peaklistFilenames, finalResultFile);
                System.exit(0);
            }

            if (commandLine.hasOption(OPTIONS.CLUSTER_MERGE_BINARY_FILES.getValue())) {
                System.out.println("Clustering and merging binary files...");
                mergeClusterBinaryFiles(peaklistFilenames, finalResultFile, paralellJobs, thresholds, spectraClusterStandalone.isVerbose());
                System.exit(0);
            }

            /**
             * ------- THE ACTUAL LOGIC STARTS HERE -----------
             */
            printSettings(finalResultFile, paralellJobs, startThreshold, endThreshold, rounds, spectraClusterStandalone.isKeepBinaryFiles(),
                    spectraClusterStandalone.getTemporaryDirectory(), peaklistFilenames, reUseBinaryFiles, spectraClusterStandalone.isUseFastMode(),
                    addedFilters);

            spectraClusterStandalone.addProgressListener(this);

            // merge mode
            if (mergeBinaryFilesMode) {
                System.out.println("Binary file merging mode set.");

                if (commandLine.hasOption(CliOptions.OPTIONS.ADD_SCORES.getValue())) {
                    System.out.println("Error: Scores cannot be added in binary file merging mode");
                    System.exit(1);
                }

                spectraClusterStandalone.mergeBinaryFiles(peaklistFilenames, thresholds, finalResultFile);
                System.exit(0);
            }

            // make sure binary files exist
            if (reUseBinaryFiles) {
                File binaryFileDirectory = new File(spectraClusterStandalone.getTemporaryDirectory(), "spectra");
                if (!binaryFileDirectory.isDirectory()) {
                    reUseBinaryFiles = false;
                    System.out.println("No binary files found. Re-creating them...");
                }
            }

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

            // add the scores if set
            if (commandLine.hasOption(CliOptions.OPTIONS.ADD_SCORES.getValue())) {
                // get the directories of the MGF files
                Set<File> directories = new HashSet<>();
                for (String peakListFile : peaklistFilenames) {
                    directories.add(new File(peakListFile).getParentFile());
                }
                File[] mgfDirs = new File[directories.size()];
                directories.toArray(mgfDirs);

                System.out.println("Adding scores to output file...");
                ScoreCalculator scoreCalculator = new ScoreCalculator();
                scoreCalculator.processClusteringResult(finalResultFile, mgfDirs);
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

    /**
     * Merges and cluster the passed binary files.
     * @param peaklistFilenames Names of the binary files to process
     * @param finalResultFile Path to the final result file
     * @param paralellJobs Number of parallel jobs
     * @param thresholds List of thresholds to use
     * @param verbose Indicates whether verbose output should be used.
     */
    private void mergeClusterBinaryFiles(String[] peaklistFilenames, File finalResultFile, int paralellJobs, List<Float> thresholds, boolean verbose) throws Exception {
        if (!finalResultFile.isDirectory()) {
            throw new Exception("Error: Result file must be a directory.");
        }

        // add max verbosity to the result files
        Defaults.setSaveAddingScore(true);
        Defaults.setSaveDebugInformation(true);

        // create the temporary directory
        File tmpDir = SpectraClusterStandalone.createTemporaryDirectory("spectra_cluster_cli");

        // count the spectra per bin while scanning the files
        SpectraPerBinNumberComparisonAssessor spectraPerBinNumberComparisonAssessor = null;
        if (Defaults.getNumberOfComparisonAssessor().getClass() == SpectraPerBinNumberComparisonAssessor.class) {
            spectraPerBinNumberComparisonAssessor = (SpectraPerBinNumberComparisonAssessor) Defaults.getNumberOfComparisonAssessor();
        }

        // scan the binary files
        List<File> filenames = Arrays.stream(peaklistFilenames).map(File::new).collect(Collectors.toList());
        List<BinaryClusterFileReference> clusterReferences = BinaryFileScanner.scanBinaryFiles(
                spectraPerBinNumberComparisonAssessor,
                null,
                filenames.toArray(new File[filenames.size()]));

        // Create the merger
        BinaryFileMergingClusterer merger = new BinaryFileMergingClusterer(paralellJobs, finalResultFile, thresholds,
                false, Defaults.getDefaultPrecursorIonTolerance(), false, tmpDir);

        if (verbose)
            merger.addProgressListener(this);

        // launch the merging
        merger.clusterFiles(clusterReferences);

        System.out.println("Result files written to " + finalResultFile.toString());
    }

    /**
     * Converts the passed binary files to the cgf format
     * @param peaklistFilenames
     * @param finalResultFile
     */
    private void convertBinaryFiles(String[] peaklistFilenames, File finalResultFile) throws Exception {
        FileWriter writer = new FileWriter(finalResultFile);

        String[] options = {"PrideCluster 2.0"};

        DotClusterClusterAppender.INSTANCE.appendStart(writer, options);

        for (String binaryFile : peaklistFilenames) {
            System.out.println("Converting " + binaryFile + "...");

            BinaryClusterIterable binaryClusterIterable = new BinaryClusterIterable(new ObjectInputStream(new FileInputStream(binaryFile)));

            for (ICluster cluster : binaryClusterIterable) {
                DotClusterClusterAppender.INSTANCE.appendCluster(writer, cluster);
            }
        }

        DotClusterClusterAppender.INSTANCE.appendEnd(writer);

        System.out.println("Result written to " + finalResultFile.toString());
    }

    /**
     * Only cluster the passed binary files without merging them.
     * @param peaklistFilenames
     * @param finalResultFile
     */
    private void clusterBinaryFile(String[] peaklistFilenames, File finalResultFile, int nJobs, List<Float> clusteringThresholds, boolean verbose) throws Exception {
        if (!finalResultFile.isDirectory()) {
            throw new Exception("Error: Result file must be a directory.");
        }

        // add max verbosity to the result files
        Defaults.setSaveAddingScore(true);
        Defaults.setSaveDebugInformation(true);

        // create the temporary directory
        File tmpDir = SpectraClusterStandalone.createTemporaryDirectory("spectra_cluster_cli");

        // count the spectra per bin while scanning the files
        SpectraPerBinNumberComparisonAssessor spectraPerBinNumberComparisonAssessor = null;
        if (Defaults.getNumberOfComparisonAssessor().getClass() == SpectraPerBinNumberComparisonAssessor.class) {
            spectraPerBinNumberComparisonAssessor = (SpectraPerBinNumberComparisonAssessor) Defaults.getNumberOfComparisonAssessor();
        }

        // scan the binary files
        List<File> filenames = Arrays.stream(peaklistFilenames).map(File::new).collect(Collectors.toList());
        List<BinaryClusterFileReference> clusterReferences = BinaryFileScanner.scanBinaryFiles(
                spectraPerBinNumberComparisonAssessor,
                null,
                filenames.toArray(new File[filenames.size()]));

        // set the cluster predicate
        IPredicate<ICluster> clusterPredicate;

        switch (ClusteringSettings.getLoadingMode()) {
            case ONLY_IDENTIFIED:
                clusterPredicate = new ClusterOnlyIdentifiedPredicate();
                break;
            case ONLY_UNIDENTIFIED:
                clusterPredicate = new ClusterOnlyUnidentifiedPredicate();
                break;
            default:
                clusterPredicate = null;
                break;
        }

        BinaryFileClusterer binaryFileClusterer = new BinaryFileClusterer(nJobs, finalResultFile,
                clusteringThresholds, false, tmpDir, clusterPredicate);

        System.out.println(String.format("Clustering %d binary files...", peaklistFilenames.length));

        if (verbose) {
            System.out.println("Verbose mode active");
            binaryFileClusterer.addProgressListener(this);
        }

        // start the clustering
        binaryFileClusterer.clusterFiles(clusterReferences);

        // TODO: delete temporary directory
    }

    /**
     * Merge the passed binary result files into a single binary file.
     * @param peaklistFilenames
     * @param finalResultFile
     */
    private void mergeBinaryFiles(String[] peaklistFilenames, File finalResultFile) throws Exception {
        // Create the output file
        ObjectOutputStream outputStream = new ObjectOutputStream(new FileOutputStream(finalResultFile));

        // open the files
        List<Iterator<ICluster>> clusterIterables = new ArrayList<>(peaklistFilenames.length);

        for (String filename : peaklistFilenames) {
            clusterIterables.add(new BinaryClusterIterable(new ObjectInputStream(new FileInputStream(filename))).iterator());
        }

        // read the first clusters from the files
        Map<Integer, ICluster> currentClusters = new HashMap<>(peaklistFilenames.length);

        for (int i = 0; i < clusterIterables.size(); i++) {
            Iterator<ICluster> iterator = clusterIterables.get(i);
            currentClusters.put(i, iterator.next());
        }

        // always write the lowest cluster to the output file
        while (true) {
            int currentLowestIndex = -1;
            float currentLowestMz = Float.MAX_VALUE;

            for (int i = 0; i < currentClusters.size(); i++) {
                ICluster cluster = currentClusters.get(i);

                if (cluster == null) {
                    continue;
                }

                if (cluster.getPrecursorMz() < currentLowestMz) {
                    currentLowestIndex = i;
                    currentLowestMz = cluster.getPrecursorMz();
                }
            }

            // if the currentLowestIndex wasn't updated, all iterators are done
            if (currentLowestIndex == -1) {
                break;
            }

            // write the lowest cluster and read the next cluster from that file
            BinaryClusterAppender.INSTANCE.appendCluster(outputStream, currentClusters.get(currentLowestIndex));

            Iterator<ICluster> iterator = clusterIterables.get(currentLowestIndex);
            currentClusters.put(currentLowestIndex, iterator.hasNext() ? iterator.next() : null);
        }

        BinaryClusterAppender.INSTANCE.appendEnd(outputStream);
        outputStream.close();
    }

    /**
     * Convert MGF files into the binary format (first stage of clustering)
     * @param peaklistFilenames
     * @param temporaryDirectory
     */
    private void createBinaryFiles(String[] peaklistFilenames, File temporaryDirectory, int nJobs) throws Exception {
        if (!temporaryDirectory.isDirectory()) {
            throw new Exception("Error: Output path must be set to a directory");
        }

        // if the binary spectra directory doesn't exist, create it
        if (!temporaryDirectory.exists()) {
            if (!temporaryDirectory.mkdir()) {
                throw new Exception("Failed to create temporary binary directory: " + temporaryDirectory);
            }
        }

        // always bin to 2 m/z wide bins
        int windowSize = 2;

        BinningSpectrumConverter binningSpectrumConverter = new BinningSpectrumConverter(temporaryDirectory,
                nJobs, false, new FixedReferenceMzBinner(windowSize));

        binningSpectrumConverter.processPeaklistFiles(peaklistFilenames);

        System.out.println("Binary files written to " + temporaryDirectory.toString());
    }

    /**
     * Convert the passed files from the .cgf format to the .clustering format.
     * @param inputFiles
     * @param outputFile
     */
    private void convertToCgf(String[] inputFiles, File outputFile) throws Exception {
        SpectraClusterStandalone spectraClusterStandalone = new SpectraClusterStandalone();
        spectraClusterStandalone.addProgressListener(this);

        if (inputFiles.length != 1) {
            System.out.println("Error: Only one cgf file can be converted to one .clustering file");
            return;
        }

        spectraClusterStandalone.convertCgfToClustering(new File(inputFiles[0]), outputFile, 0.99F);
    }

    private void printSettings(File finalResultFile, int nMajorPeakJobs, float startThreshold,
                                      float endThreshold, int rounds, boolean keepBinaryFiles, File binaryTmpDirectory,
                                      String[] peaklistFilenames, boolean reUseBinaryFiles, boolean fastMode,
                                      List<String> addedFilters) {
        System.out.println("spectra-cluster API Version 1.0.11");
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
        if (ClusteringSettings.ppmThreshold != null)
            System.out.println("Precursor tolerance: " + ClusteringSettings.ppmThreshold + " ppm");
        else
            System.out.println("Precursor tolerance: " + Defaults.getDefaultPrecursorIonTolerance() + " m/z");

        System.out.println("Fragment ion tolerance: " + Defaults.getFragmentIonTolerance());

        // loading peak filter
        if (ClusteringSettings.getLoadingSpectrumFilter().getClass() == HighestNPeakFunction.class) {
            System.out.println("Loading filter: Filtering top N peaks per spectrum");
        }
        if (ClusteringSettings.getLoadingSpectrumFilter().getClass() == BinnedHighestNPeakFunction.class) {
            System.out.println("Loading filter: Filtering top 10 peaks / 100 m/z");
        }

        // used filters
        System.out.print("Added filters: ");
        for (int i = 0; i < addedFilters.size(); i++) {
            if (i > 0) {
                System.out.print(", ");
            }
            System.out.print(addedFilters.get(i));
        }
        System.out.println("");

        // number of comparisons
        if (Defaults.getNumberOfComparisonAssessor().getClass() == MinNumberComparisonsAssessor.class) {
            MinNumberComparisonsAssessor assessor = (MinNumberComparisonsAssessor) Defaults.getNumberOfComparisonAssessor();
            System.out.println("Minimum number of comparisons: " + assessor.getMinNumberComparisons());
        } else {
            System.out.println("Minimum number of comparisons: adaptive");
        }

        System.out.println();
    }

    private void printUsage() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Spectra Cluster - Clusterer",
                "Clusters the spectra found in a MGF files or .clustering files and writes the results in a text-based file.\n",
                CliOptions.getOptions(), "\n\n", true);
    }

    @Override
    public void onProgressUpdate(ProgressUpdate progressUpdate) {
        System.out.println(progressUpdate.getMessage());
    }
}
