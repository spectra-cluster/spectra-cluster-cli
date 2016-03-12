package uk.ac.ebi.pride.spectracluster.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.FileFilterUtils;
import uk.ac.ebi.pride.spectracluster.binning.BinarySpectrumReferenceWriterCallable;
import uk.ac.ebi.pride.spectracluster.binning.ISpectrumReferenceBinner;
import uk.ac.ebi.pride.spectracluster.binning.ReferenceMzBinner;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusteringCallable;
import uk.ac.ebi.pride.spectracluster.clustering.ClusteringProcessLauncher;
import uk.ac.ebi.pride.spectracluster.conversion.MergingCGFConverter;
import uk.ac.ebi.pride.spectracluster.engine.IIncrementalClusteringEngine;
import uk.ac.ebi.pride.spectracluster.io.*;
import uk.ac.ebi.pride.spectracluster.merging.LoadingSimilarClusteringEngine;
import uk.ac.ebi.pride.spectracluster.spectra_list.*;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Created with IntelliJ IDEA.
 * User: jg
 * Date: 9/15/13
 * Time: 11:19 AM
 * To change this template use File | Settings | File Templates.
 */
public class SpectraClusterCliMain {
    public final static int MAJOR_PEAK_CLUSTERING_JOBS = 4;
    public static final boolean DELETE_TEMPORARY_CLUSTERING_RESULTS = true;

    public static void main(String[] args) {
        CommandLineParser parser = new GnuParser();

        try {
            CommandLine commandLine = parser.parse(CliOptions.getOptions(), args);

            // HELP
            if (commandLine.hasOption(CliOptions.OPTIONS.HELP.getValue())) {
                printUsage();
                return;
            }

            // RESULT FILE PATH
            if (!commandLine.hasOption(CliOptions.OPTIONS.OUTPUT_PATH.getValue()))
                throw new Exception("Missing required option " + CliOptions.OPTIONS.OUTPUT_PATH.getValue());
            File finalResultFile = new File(commandLine.getOptionValue(CliOptions.OPTIONS.OUTPUT_PATH.getValue()));

            if (finalResultFile.exists())
                throw new Exception("Result file " + finalResultFile + " already exists");

            // NUMBER OF JOBS
            int nMajorPeakJobs = MAJOR_PEAK_CLUSTERING_JOBS;
            if (commandLine.hasOption(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue())) {
                nMajorPeakJobs = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue()));
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

            List<Float> thresholds = generateThreshold(startThreshold, endThreshold, rounds);

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

            // MERGE DUPLICATE
            boolean mergeDuplicate = commandLine.hasOption(CliOptions.OPTIONS.MERGE_DUPLICATE.getValue());

            // BINARY TMP DIR
            File binaryTmpDirectory;

            if (commandLine.hasOption(CliOptions.OPTIONS.BINARY_TMP_DIR.getValue()))
                binaryTmpDirectory = new File(commandLine.getOptionValue(CliOptions.OPTIONS.BINARY_TMP_DIR.getValue()));
            else
                binaryTmpDirectory = createTemporaryDirectory("binary_converted_spectra");

            // KEEP BINARY FILES
            boolean keepBinaryFiles = commandLine.hasOption(CliOptions.OPTIONS.KEEP_BINARY_FILE.getValue());

            // RE-USE BINARY FILES
            boolean reUseBinaryFiles = commandLine.hasOption(CliOptions.OPTIONS.REUSE_BINARY_FILES.getValue());

            // FAST MODE
            boolean fastMode = commandLine.hasOption(CliOptions.OPTIONS.FAST_MODE.getValue());

            // FILES TO PROCESS
            String[] peaklistFilenames = commandLine.getArgs();

            // if re-use is set, binaryTmpDirectory is required and merging is impossible
            if (reUseBinaryFiles && !commandLine.hasOption(CliOptions.OPTIONS.BINARY_TMP_DIR.getValue()))
                throw new Exception("Missing required option '" + CliOptions.OPTIONS.BINARY_TMP_DIR.getValue() + "' with " + CliOptions.OPTIONS.REUSE_BINARY_FILES.getValue());

            if (reUseBinaryFiles && mergeDuplicate)
                throw new Exception("Merging of results is not possible if binary files are being re-used");

            if (reUseBinaryFiles && peaklistFilenames.length > 0)
                System.out.println("WARNING: " + CliOptions.OPTIONS.REUSE_BINARY_FILES.getValue() + " set, input files will be ignored");

            /**
             * SPECIAL MODES
             */
            // cluster binary file
            if (commandLine.hasOption(CliOptions.OPTIONS.CLUSTER_BINARY_FILE.getValue())) {
                clusterBinaryFile(commandLine.getOptionValue(CliOptions.OPTIONS.CLUSTER_BINARY_FILE.getValue()), finalResultFile, thresholds, fastMode);
                return;
            }

            // merge binary files
            if (commandLine.hasOption(CliOptions.OPTIONS.MERGE_BINARY_RESULTS.getValue())) {
                mergeBinaryFiles(commandLine.getArgs(), finalResultFile, mergeDuplicate);
                return;
            }

            // convert cgf
            if (commandLine.hasOption(CliOptions.OPTIONS.CONVERT_CGF.getValue())) {
                if (commandLine.getArgs().length > 1)
                    throw new Exception("Can only convert a single file at the time.");

                convertClusters(new File(commandLine.getArgs()[0]), finalResultFile, endThreshold);
                return;
            }

            /**
             * ------- THE ACTUAL LOGIC STARTS HERE -----------
             */
            printSettings(finalResultFile, nMajorPeakJobs, startThreshold, endThreshold, rounds, mergeDuplicate, keepBinaryFiles, binaryTmpDirectory, peaklistFilenames, reUseBinaryFiles, fastMode);

            List<File> binaryFiles;
            BinningSpectrumConverter binningSpectrumConverter = null;

            if (!reUseBinaryFiles) {
                System.out.print("Writing binary files...");
                long start = System.currentTimeMillis();

                binningSpectrumConverter = new BinningSpectrumConverter(binaryTmpDirectory, nMajorPeakJobs, fastMode);
                binningSpectrumConverter.processPeaklistFiles(peaklistFilenames);
                binaryFiles = binningSpectrumConverter.getWrittenFiles();

                String message = String.format("Done. Found %d spectra", binningSpectrumConverter.getSpectrumReferences().size());
                printDone(start, message);
            }
            else {
                // get the list of files
                File[] existingBinaryFiles = binaryTmpDirectory.listFiles((FilenameFilter) FileFilterUtils.suffixFileFilter(".cls"));
                binaryFiles = new ArrayList<File>(existingBinaryFiles.length);
                for (File file : existingBinaryFiles)
                    binaryFiles.add(file);

                System.out.println("Found " + binaryFiles.size() + " existing binary files.");
            }

            // create a temporary directory for the clustering results
            File tmpClusteringResults = createTemporaryDirectory("clustering_results");

            // cluster the binary files and immediately convert the results
            BinaryFileClusterer binaryFileClusterer = new BinaryFileClusterer(nMajorPeakJobs, tmpClusteringResults, thresholds, fastMode);

            File combinedResultFile = File.createTempFile("combined_clustering_results", ".cgf");
            MergingCGFConverter mergingCGFConverter = new MergingCGFConverter(combinedResultFile, DELETE_TEMPORARY_CLUSTERING_RESULTS, !keepBinaryFiles, binaryTmpDirectory);
            binaryFileClusterer.addListener(mergingCGFConverter);

            System.out.println("Clustering " + binaryFiles.size() + " binary files...");
            long start = System.currentTimeMillis();

            binaryFileClusterer.clusterFiles(binaryFiles);

            printDone(start, "Completed clustering.");

            // delete the temporary directories if set
            if (!keepBinaryFiles) {
                if (!binaryTmpDirectory.delete())
                    System.out.println("Warning: Failed to delete " + binaryTmpDirectory);
            }

            if (DELETE_TEMPORARY_CLUSTERING_RESULTS) {
                if (!tmpClusteringResults.delete())
                    System.out.println("Warning: Failed to delete " + tmpClusteringResults);
            }

            // TODO: this implementation is still missing the merging job

            // create the output file
            if (mergeDuplicate)
                mergeDuplicateClusters(combinedResultFile, mergingCGFConverter.getClusterReferences(), finalResultFile, getSpectrumReferencesPerId(binningSpectrumConverter.getSpectrumReferences()), peaklistFilenames, endThreshold, binningSpectrumConverter.getFileIndices());
            else
                convertClusters(combinedResultFile, finalResultFile, endThreshold);

            if (!combinedResultFile.delete())
                System.out.println("Warning: Failed to delete " + combinedResultFile);

            System.out.println("Clustering completed. Results written to " + finalResultFile);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Error: " + e.getMessage());

            System.exit(1);
        }
    }

    private static void mergeBinaryFiles(String[] binaryFilenames, File finalResultFile, boolean mergeDuplicate) throws Exception {
        File tmpResultFile = File.createTempFile("combined_results", ".cgf");

        MergingCGFConverter mergingCGFConverter = new MergingCGFConverter(tmpResultFile, false, false, null); // do not delete any files

        System.out.print("Merging " + binaryFilenames.length + " binary files...");
        long start = System.currentTimeMillis();

        for (String binaryFilename : binaryFilenames) {
            mergingCGFConverter.onNewResultFile(new File(binaryFilename));
        }

        printDone(start);

        // copy the temporary file to the final destination
        start = System.currentTimeMillis();
        System.out.print("Copying result to " + finalResultFile + "...");

        FileUtils.copyFile(tmpResultFile, finalResultFile);

        printDone(start);
        tmpResultFile.delete();
    }

    /**
     * Clusters the passed binary file in a single thread and writes the result to "binaryFilename"
     * @param binaryFilename
     * @param finalResultFile
     * @param thresholds
     */
    private static void clusterBinaryFile(String binaryFilename, File finalResultFile, List<Float> thresholds, boolean fastMode) throws Exception {
        System.out.println("spectra-cluster API Version 1.0");
        System.out.println("Created by Rui Wang & Johannes Griss\n");

        System.out.println("Clustering single binary file: " + binaryFilename);
        System.out.println("Result file: " + finalResultFile);

        // write to a (local) temporary file
        File tmpResultFile = File.createTempFile("clustering_result", ".cls");

        long start = System.currentTimeMillis();
        System.out.print("Clustering file...");

        BinaryFileClusteringCallable binaryFileClusteringCallable = new BinaryFileClusteringCallable(tmpResultFile, new File(binaryFilename), thresholds, fastMode);
        binaryFileClusteringCallable.call();

        printDone(start);

        // copy the file
        System.out.println("Copying result file to " + finalResultFile);
        FileUtils.copyFile(tmpResultFile, finalResultFile);

        tmpResultFile.delete();
    }

    private static void printSettings(File finalResultFile, int nMajorPeakJobs, float startThreshold, float endThreshold, int rounds, boolean mergeDuplicate, boolean keepBinaryFiles, File binaryTmpDirectory, String[] peaklistFilenames, boolean reUseBinaryFiles, boolean fastMode) {
        System.out.println("spectra-cluster API Version 1.0");
        System.out.println("Created by Rui Wang & Johannes Griss\n");

        System.out.println("-- Settings --");
        System.out.println("Number of threads: " + String.valueOf(nMajorPeakJobs));
        System.out.println("Thresholds: " + String.valueOf(startThreshold) + " - " + String.valueOf(endThreshold) + " in " + rounds + " rounds");
        System.out.println("Merging duplicate: " + (mergeDuplicate ? "true" : "false"));
        System.out.println("Keeping binary files: " + (keepBinaryFiles ? "true" : "false"));
        System.out.println("Binary file directory: " + binaryTmpDirectory);
        System.out.println("Result file: " + finalResultFile);
        System.out.println("Reuse binary files: " + (reUseBinaryFiles ? "true" : "false"));
        System.out.println("Input files: " + peaklistFilenames.length);
        System.out.println("Using fast mode: " + (fastMode ? "yes" : "no"));

        System.out.println();
    }

    /**
     * Generates the actual thresholds to use based on the
     * start and end threshold and the number of iterations
     * to perform. The result is sorted from highest to
     * lowest threshold.
     *
     * @param startThreshold Starting threshold (highest threshold)
     * @param endThreshold Final threshold
     * @param rounds Number of rounds of clustering to perform
     * @return
     */
    private static List<Float> generateThreshold(float startThreshold, float endThreshold, int rounds) {
        List<Float> thresholds = new ArrayList<Float>(rounds);
        float stepSize = (startThreshold - endThreshold) / (rounds - 1);

        for (int i = 0; i < rounds; i++) {
            thresholds.add(startThreshold - (stepSize * i));
        }

        return thresholds;
    }

    private static void convertClusters(File combinedResultFile, File finalResultFile, float endThreshold) throws Exception {
        FileInputStream fileInputStream = new FileInputStream(combinedResultFile);
        FileWriter writer = new FileWriter(finalResultFile);

        int nClusterWritten = 0;
        System.out.print("Converting results to .clustering...");
        long start = System.currentTimeMillis();

        DotClusterClusterAppender.INSTANCE.appendStart(writer, "GreedyClustering_" + String.valueOf(endThreshold));
        CGFSpectrumIterable iterable = new CGFSpectrumIterable(fileInputStream);

        for (ICluster cluster : iterable) {
            DotClusterClusterAppender.INSTANCE.appendCluster(writer, cluster);
            nClusterWritten++;
        }

        writer.close();

        System.out.println("Done. (Took " + String.format("%.2f", (float) (System.currentTimeMillis() - start) / 1000 / 60) + " min. " + nClusterWritten + " final clusters)");

        // copy the final CGF file as well
        FileUtils.copyFile(combinedResultFile, new File(finalResultFile.getAbsolutePath() + ".cgf"));
        System.out.println("Copied CGF result to " + finalResultFile.getAbsolutePath() + ".cgf");
    }

    private static void mergeDuplicateClusters(File combinedResultFile, List<ClusterReference> clusterReferences, File finalResultFile, Map<String, SpectrumReference> spectrumReferencesPerId, String[] peaklistFilenames, float finalThreshold, List<List<IndexElement>> fileIndices) throws Exception {
        FileInputStream fileInputStream = new FileInputStream(combinedResultFile);

        FileWriter writer = new FileWriter(finalResultFile);
        final float WINDOW_SIZE = 1.0F;
        final double MIN_SHARED_SPECTRA = 0.4;

        int nClusterRead = 0;
        int nClusterWritten = 0;

        // process all clusters based on precursor m/z
        Collections.sort(clusterReferences);

        System.out.print("Merging duplicate clusters...");
        long start = System.currentTimeMillis();

        IIncrementalClusteringEngine clusteringEngine = new LoadingSimilarClusteringEngine(Defaults.getDefaultSpectrumComparator(), WINDOW_SIZE, MIN_SHARED_SPECTRA, spectrumReferencesPerId, peaklistFilenames, fileIndices);
        DotClusterClusterAppender.INSTANCE.appendStart(writer, "GreedyClustering_" + String.valueOf(finalThreshold)); // TODO: add better name

        for (ClusterReference clusterReference : clusterReferences) {
            nClusterRead++;

            // load the cluster
            fileInputStream.getChannel().position(clusterReference.getOffset());
            LineNumberReader reader = new LineNumberReader(new InputStreamReader(fileInputStream));
            ICluster cluster = ParserUtilities.readSpectralCluster(reader, null);

            if (!clusterReference.getId().equals(cluster.getId()))
                throw new IllegalStateException("Wrong cluster read");

            // if a cluster is merged, load its spectra from the binary file first
            Collection<ICluster> removedCluster = clusteringEngine.addClusterIncremental(cluster);

            if (!removedCluster.isEmpty()) {
                for (ICluster rc : removedCluster) {
                    nClusterWritten++;
                    DotClusterClusterAppender.INSTANCE.appendCluster(writer, rc);
                }
            }
        }

        Collection<ICluster> remainingClusters = clusteringEngine.getClusters();
        for (ICluster rc : remainingClusters) {
            nClusterWritten++;
            DotClusterClusterAppender.INSTANCE.appendCluster(writer, rc);
        }

        DotClusterClusterAppender.INSTANCE.appendEnd(writer);
        writer.close();
        System.out.println("Done. (Took " + String.format("%.2f", (float) (System.currentTimeMillis() - start) / 1000 / 60) + " min. Reduced " + nClusterRead + " to " + nClusterWritten + " final clusters)");
    }

    private static Map<String, SpectrumReference> getSpectrumReferencesPerId(List<SpectrumReference> spectrumReferences) {
        Map<String, SpectrumReference> spectrumReferencePerId = new HashMap<String, SpectrumReference>();

        // save the spectrum references per id
       for (SpectrumReference spectrumReference : spectrumReferences) {
            if (spectrumReferencePerId.containsKey(spectrumReference.getSpectrumId()))
                continue;

            spectrumReferencePerId.put(spectrumReference.getSpectrumId(), spectrumReference);
        }

        return spectrumReferencePerId;
    }

    private static void printDone(long start) {
        printDone(start, "Done");
    }

    private static void printDone(long start, String message) {
        long duration = System.currentTimeMillis() - start;
        System.out.println(message + " (" + String.format("%.2f", (double) duration / 1000 / 60) + " min)");
    }

    private static File createTemporaryDirectory(String prefix) throws Exception {
        File tmpFile = File.createTempFile(prefix, "");

        if (!tmpFile.delete())
            throw new Exception("Failed to delete temporary file");

        if (!tmpFile.mkdir())
            throw new Exception("Failed to create temporary directory");

        return tmpFile;
    }

    private static void printUsage() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Spectra Cluster - Clusterer",
                "Clusters the spectra found in an MGF file and writes the results in a text-based file.\n",
                CliOptions.getOptions(), "\n\n", true);
    }


}
