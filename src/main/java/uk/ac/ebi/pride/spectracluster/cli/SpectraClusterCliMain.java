package uk.ac.ebi.pride.spectracluster.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.io.FileUtils;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.ClusteringProcessLauncher;
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
    private static List<List<IndexElement>> fileIndices;

    public static void main(String[] args) {
        CommandLineParser parser = new GnuParser();

        try {
            CommandLine commandLine = parser.parse(CliOptions.getOptions(), args);

            // HELP
            if (commandLine.hasOption(CliOptions.OPTIONS.HELP.getValue())) {
                printUsage();
                return;
            }

            if (!commandLine.hasOption(CliOptions.OPTIONS.OUTPUT_PATH.getValue()))
                throw new Exception("Missing required option " + CliOptions.OPTIONS.OUTPUT_PATH.getValue());
            File finalResultFile = new File(commandLine.getOptionValue(CliOptions.OPTIONS.OUTPUT_PATH.getValue()));

            if (finalResultFile.exists())
                throw new Exception("Result file " + finalResultFile + " already exists");

            int nMajorPeakJobs = MAJOR_PEAK_CLUSTERING_JOBS;
            if (commandLine.hasOption(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue())) {
                nMajorPeakJobs = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.MAJOR_PEAK_JOBS.getValue()));
            }

            int rounds = 4;
            if (commandLine.hasOption(CliOptions.OPTIONS.ROUNDS.getValue()))
                rounds = Integer.parseInt(commandLine.getOptionValue(CliOptions.OPTIONS.ROUNDS.getValue()));

            float startThreshold = 0.999F;
            if (commandLine.hasOption(CliOptions.OPTIONS.START_THRESHOLD.getValue()))
                startThreshold = Float.parseFloat(commandLine.getOptionValue(CliOptions.OPTIONS.START_THRESHOLD.getValue()));

            float endThreshold = 0.99F;
            if (commandLine.hasOption(CliOptions.OPTIONS.END_THRESHOLD.getValue()))
                endThreshold = Float.parseFloat(commandLine.getOptionValue(CliOptions.OPTIONS.END_THRESHOLD.getValue()));

            List<Float> thresholds = new ArrayList<Float>(rounds);
            float stepSize = (startThreshold - endThreshold) / (rounds - 1);

            for (int i = 0; i < rounds; i++) {
                thresholds.add(startThreshold - (stepSize * i));
            }

            boolean mergeDuplicate = commandLine.hasOption(CliOptions.OPTIONS.MERGE_DUPLICATE.getValue());

            // 1.) Pre-process spectra and store in list per highest peak
            /**
             * pre-scan all files and create list of SpectrumReferences for all spectra
             * > already store each spectrum N times for each highest peak
             *   > NOTE: this already requires the precursor filter to be applied
             * sort each list according to precursor m/z
             * write spectra to files in this specific order (use ObjectOutputStream)
             */
            String[] peaklistFilenames = commandLine.getArgs();
            File tmpSpectraPerPeakDir = createTemporaryDirectory("spectra_per_peak");
            File tmpClusteredPeakDir = createTemporaryDirectory("clustering_results");

            List<SpectrumReference> spectrumReferences = prescanPeaklistFiles(peaklistFilenames);
            Map<String, SpectrumReference> spectrumReferencesPerId = getSpectrumReferencesPerId(spectrumReferences);

            // write out the spectra, one file per major peak and start the clustering job right away
            System.out.println("Converting spectra to binary format per major peak");

            // create the spectrum writer and add the executor service as listener,
            // thereby, as soon as a file is written, the clustering job is launched
            SpectrumWriter spectrumWriter = new SpectrumWriter(peaklistFilenames, fileIndices);
            ExecutorService executorService = Executors.newFixedThreadPool(nMajorPeakJobs);
            ClusteringProcessLauncher clusteringProcessLauncher = new ClusteringProcessLauncher(executorService, tmpClusteredPeakDir, thresholds);
            spectrumWriter.addListener(clusteringProcessLauncher);

            // write the major peak files
            /**
            for (int majorPeak : spectrumReferencesPerMajorPeak.keySet()) {
                File outputFile = getMajorPeakSourceFile(majorPeak, tmpSpectraPerPeakDir);

                spectrumWriter.writeSpectra(spectrumReferencesPerMajorPeak.get(majorPeak), outputFile);
            }
             */

            // group the spectrum references and write each group to a file
            ReferenceMzBinner binner = new ReferenceMzBinner();
            List<List<SpectrumReference>> groupedSpectrumReferences = binner.groupSpectrumReferences(spectrumReferences);
            System.out.println("Split " + spectrumReferences.size() + " spectra in " + groupedSpectrumReferences.size() + " bins.");
            int outputIndex = 0;
            for (List<SpectrumReference> spectrumReferenceList : groupedSpectrumReferences) {
                File outputFile = getMajorPeakSourceFile(outputIndex, tmpSpectraPerPeakDir);

                spectrumWriter.writeSpectra(spectrumReferenceList, outputFile);
                outputIndex++;
            }

            System.out.println("Completed writing binary files.");

            // wait until all clustering jobs are done - since all files were written, all
            // jobs have been submitted
            // start the termination process and merge the results into one file
            executorService.shutdown();
            boolean allDone = false;

            // write all clusters into on cgf and save each cluster's position
            File combinedResultFile = File.createTempFile("combined_clustering_results", ".cgf");
            List<ClusterReference> clusterReferences = new ArrayList<ClusterReference>();
            List<Future<File>> resultFileFutures = clusteringProcessLauncher.getFileFutures();

            // wait until all jobs are done
            Set<Integer> completedJobs = new HashSet<Integer>();

            while (!allDone) {
                allDone = true;

                for (int i = 0; i < resultFileFutures.size(); i++) {
                    if (completedJobs.contains(i))
                        continue;

                    Future<File> fileFuture = resultFileFutures.get(i);

                    if (fileFuture.isDone()) {
                        // scan the result file - this is done in the main thread to make sure that only one scanning process is running at a time
                        List<ClusterReference> completedClusterReferences = mergeClusteringResults(fileFuture.get(), combinedResultFile);
                        clusterReferences.addAll(completedClusterReferences);

                        // remove the result file
                        fileFuture.get().delete();
                        // remove major peak file
                        File majorPeakFile = new File(tmpSpectraPerPeakDir, fileFuture.get().getName());
                        majorPeakFile.delete();

                        completedJobs.add(i);
                    }
                    else {
                        allDone = false;
                    }
                }
            }

            executorService.awaitTermination(1, TimeUnit.MINUTES);

            if (!tmpSpectraPerPeakDir.delete())
                System.out.println("Warning: Failed to delete " + tmpSpectraPerPeakDir);

            if (!tmpClusteredPeakDir.delete())
                System.out.println("Warning: Failed to delete " + tmpClusteredPeakDir);

            // 3.) merge duplicated clusters
            /**
             * pre-scan all generated clusters
             *  > store in one big list, sort according to m/z
             * run incremental clustering directly on this list / big file
             */

            /**
             * Merge all clusters that share more than 40% of their spectra
             */
            if (mergeDuplicate)
                mergeDuplicateClusters(combinedResultFile, clusterReferences, finalResultFile, spectrumReferencesPerId, peaklistFilenames, endThreshold);
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

    private static List<ClusterReference> mergeClusteringResults(File binaryResultFile, File combinedResultFile) throws Exception {
        FileOutputStream outputStream = new FileOutputStream(combinedResultFile, true);

        ObjectInputStream objectInputStream = new ObjectInputStream(new FileInputStream(binaryResultFile));
        BinaryClusterIterable binaryClusterIterable = new BinaryClusterIterable(objectInputStream);

        List<ClusterReference> clusterReferences = new ArrayList<ClusterReference>();

        for (ICluster cluster : binaryClusterIterable) {
            long offset = outputStream.getChannel().position();
            OutputStreamWriter writer = new OutputStreamWriter(outputStream);
            CGFClusterAppender.INSTANCE.appendCluster(writer, cluster);
            writer.flush();

            // save the position of the cluster
            ClusterReference clusterReference = new ClusterReference(
                    0, offset, cluster.getPrecursorMz(), cluster.getClusteredSpectraCount(), cluster.getId());
            clusterReferences.add(clusterReference);
        }

        outputStream.close();

        return clusterReferences;
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

    private static void mergeDuplicateClusters(File combinedResultFile, List<ClusterReference> clusterReferences, File finalResultFile, Map<String, SpectrumReference> spectrumReferencesPerId, String[] peaklistFilenames, float finalThreshold) throws Exception {
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

    @Deprecated // random access doesn't work on a binary file
    private static void mergeClusteringResults(List<ClusterReference> clusterReferences, List<File> clusteringResultFiles, File combinedResultFile) throws Exception {
        // open the result file
        FileOutputStream outputStream = new FileOutputStream(combinedResultFile);
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(outputStream);

        // sort the cluster references, thereby clusters are written in ascending precursor m/z
        Collections.sort(clusterReferences);

        for (ClusterReference clusterReference : clusterReferences) {
            int fileId = clusterReference.getFileId();

            if (fileId < 0 || fileId >= clusteringResultFiles.size())
                throw new Exception("Invalid file id encountered: " + fileId);

            // open the result file
            FileInputStream inputStream = new FileInputStream(clusteringResultFiles.get(fileId));
            ObjectInputStream objectInputStream = new ObjectInputStream(inputStream);

            // move to the beginning of the current cluster
            inputStream.getChannel().position(clusterReference.getOffset());

            // read the cluster at that position
            ICluster loadedCluster = BinaryClusterParser.INSTANCE.parseNextCluster(objectInputStream, null);

            inputStream.close();

            // write the cluster to the output file
            BinaryClusterAppender.INSTANCE.appendCluster(objectOutputStream, loadedCluster);
        }

        BinaryClusterAppender.INSTANCE.appendEnd(objectOutputStream);
        objectOutputStream.close();
    }

    private static Map<Integer, List<SpectrumReference>> prescanPeaklistFilesPerMajorPeak(String[] peaklistFilenames) throws Exception {
        // pre-scan all files
        System.out.print("Pre-scanning " + peaklistFilenames.length + " input files...");
        long start = System.currentTimeMillis();

        PeakListFileScanner fileScanner = new PeakListFileScanner();
        Map<Integer, List<SpectrumReference>> loadedSpectrumReferenceMap = fileScanner.getSpectraPerMajorPeaks(peaklistFilenames, 5);
        fileIndices = fileScanner.getFileIndices();

        printDone(start);

        return loadedSpectrumReferenceMap;
    }

    private static List<SpectrumReference> prescanPeaklistFiles(String[] peaklistFilenames) throws Exception {
        // pre-scan all files
        System.out.print("Pre-scanning " + peaklistFilenames.length + " input files...");
        long start = System.currentTimeMillis();

        //IPeaklistScanner fileScanner = new PeakListFileScanner();
        IPeaklistScanner fileScanner = new ParsingMgfScanner();
        List<SpectrumReference> spectrumReferences = fileScanner.getSpectrumReferences(peaklistFilenames);
        fileIndices = fileScanner.getFileIndices();

        printDone(start);

        return spectrumReferences;
    }

    private static Map<String, SpectrumReference> getSpectrumReferencesPerId(Map<Integer, List<SpectrumReference>> spectrumReferencesPerMajorPeak) {
        Map<String, SpectrumReference> spectrumReferencePerId = new HashMap<String, SpectrumReference>();

        // save the spectrum references per id
        for (List<SpectrumReference> spectrumReferences : spectrumReferencesPerMajorPeak.values()) {
            for (SpectrumReference spectrumReference : spectrumReferences) {
                if (spectrumReferencePerId.containsKey(spectrumReference.getSpectrumId()))
                    continue;

                spectrumReferencePerId.put(spectrumReference.getSpectrumId(), spectrumReference);
            }
        }

        return spectrumReferencePerId;
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

    private static File getMajorPeakSourceFile(int majorPeak, File dir) {
        File outputFile = new File(dir, "bin_" + majorPeak + ".cls");

        return outputFile;
    }

    private static void printDone(long start) {
        long duration = System.currentTimeMillis() - start;
        System.out.println("Done (" + String.format("%.2f", (double) duration / 1000 / 60) + " min.)");
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
