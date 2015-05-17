package uk.ac.ebi.pride.spectracluster.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusteringCallable;
import uk.ac.ebi.pride.spectracluster.engine.IIncrementalClusteringEngine;
import uk.ac.ebi.pride.spectracluster.engine.SimilarClusterMergingEngine;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterIterable;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterParser;
import uk.ac.ebi.pride.spectracluster.io.DotClusterClusterAppender;
import uk.ac.ebi.pride.spectracluster.merging.LoadingSimilarClusteringEngine;
import uk.ac.ebi.pride.spectracluster.spectra_list.*;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

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
    public final static int MAJOR_PEAK_CLUSTERING_JOBS = 2;

    public static void main(String[] args) {
        CommandLineParser parser = new GnuParser();

        try {
            CommandLine commandLine = parser.parse(CliOptions.getOptions(), args);

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

            File finalResultFile = new File("/tmp/clustering_result.clustering");
            if (finalResultFile.exists())
                throw new Exception("Result file " + finalResultFile + " already exists");

            Set<Integer> majorPeaks = writeSpectraPerPeakFiles(tmpSpectraPerPeakDir, peaklistFilenames);


            // 2.) Cluster each list in one thread
            /**
             * For each binary file one thread is started
             * Cluster spectra using decreasing similarity thresholds
             * ?write to CGF file or return complete list > may cause memory issues
             */
            ExecutorService executorService = Executors.newFixedThreadPool(MAJOR_PEAK_CLUSTERING_JOBS);
            List<Future<File>> resultFileFutures = new ArrayList<Future<File>>();
            for (int majorPeak : majorPeaks) {
                // build the filenames
                File inputFile = getMajorPeakSourceFile(majorPeak, tmpSpectraPerPeakDir);
                File outputFile = getMajorPeakSourceFile(majorPeak, tmpClusteredPeakDir);

                BinaryFileClusteringCallable clusteringCallable = new BinaryFileClusteringCallable(outputFile, inputFile);
                Future<File> fileFuture = executorService.submit(clusteringCallable);
                resultFileFutures.add(fileFuture);
            }

            // start the termination process
            executorService.shutdown();
            boolean allDone = false;
            List<File> clusteringResultFiles = new ArrayList<File>();
            List<ClusterReference> clusterReferences = new ArrayList<ClusterReference>();
            IClusterScanner clusterScanner = new BinaryClusterFileScanner();
            Set<Integer> completedJobs = new HashSet<Integer>();

            while (!allDone) {
                allDone = true;

                for (int i = 0; i < resultFileFutures.size(); i++) {
                    if (completedJobs.contains(i))
                        continue;

                    Future<File> fileFuture = resultFileFutures.get(i);

                    if (fileFuture.isDone()) {
                        clusteringResultFiles.add(fileFuture.get());
                        // scan the result file - this is done in the main thread to make sure that only one scanning process is running at a time
                        List<ClusterReference> completedClusterReferences = clusterScanner.getClusterReferences(
                                fileFuture.get(), clusteringResultFiles.size() - 1);
                        clusterReferences.addAll(completedClusterReferences);

                        completedJobs.add(i);
                    }
                    else {
                        allDone = false;
                    }
                }
            }

            executorService.awaitTermination(1, TimeUnit.MINUTES);

            // delete the peak list files
            for (int majorPeak : majorPeaks) {
                File majorPeakFile = getMajorPeakSourceFile(majorPeak, tmpSpectraPerPeakDir);

                if (majorPeakFile.exists()) {
                    if (!majorPeakFile.delete())
                        throw new Exception("Failed to delete " + majorPeakFile);
                }
            }
            if (!tmpSpectraPerPeakDir.delete())
                System.out.println("Warning: Failed to delete " + tmpSpectraPerPeakDir);

            // 3.) merge duplicated clusters
            /**
             * pre-scan all generated clusters
             *  > store in one big list, sort according to m/z
             * run incremental clustering directly on this list / big file
             */
            File combinedResultFile = File.createTempFile("combined_clustering_results", ".cls");

            System.out.print("Merging clustering results...");
            long start = System.currentTimeMillis();

            mergeClusteringResults(clusterReferences, clusteringResultFiles, combinedResultFile);

            printDone(start);

            // delete the temporary files
            for (File file : clusteringResultFiles) {
                if (file.exists()) {
                    if (!file.delete())
                        throw new Exception("Failed to delete " + file.toString());
                }
            }
            if (!tmpClusteredPeakDir.delete())
                System.out.println("Warning: Failed to delete " + tmpClusteredPeakDir);

            /**
             * Merge all clusters that share more than 40% of their spectra
             */
            mergeDuplicateClusters(combinedResultFile, finalResultFile);

            if (!combinedResultFile.delete())
                System.out.println("Warning: Failed to delete " + combinedResultFile);

            System.out.println("Clustering completed. Results written to " + finalResultFile);

            // HELP
            if (commandLine.hasOption(CliOptions.OPTIONS.HELP.getValue())) {
                printUsage();
                return;
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Error: " + e.getMessage());

            System.exit(1);
        }
    }

    private static void mergeDuplicateClusters(File combinedResultFile, File finalResultFile) throws Exception {
        FileInputStream fileInputStream = new FileInputStream(combinedResultFile);
        ObjectInputStream objectInputStream = new ObjectInputStream(fileInputStream);
        BinaryClusterIterable clusterIterable = new BinaryClusterIterable(objectInputStream);

        FileWriter writer = new FileWriter(finalResultFile);
        final float WINDOW_SIZE = 1.0F;
        final double MIN_SHARED_SPECTRA = 0.4;

        int nClusterRead = 0;
        int nClusterWritten = 0;

        System.out.print("Merging duplicate clusters...");
        long start = System.currentTimeMillis();

        IIncrementalClusteringEngine clusteringEngine = new LoadingSimilarClusteringEngine(Defaults.getDefaultSpectrumComparator(), WINDOW_SIZE, MIN_SHARED_SPECTRA);
        DotClusterClusterAppender.INSTANCE.appendStart(writer, "GreedyClustering"); // TODO: add better name

        for (ICluster cluster : clusterIterable) {
            nClusterRead++;
            // TODO: Greedy clusters can't really be merged like that..
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
        objectInputStream.close();

        System.out.println("Done. (Took " + String.format("%.2f", (float) (System.currentTimeMillis() - start) / 1000 / 60) + " min. Reduced " + nClusterRead + " to " + nClusterWritten + " final clusters)");
    }

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

    private static Set<Integer> writeSpectraPerPeakFiles(File tmpSpectraPerPeakDir, String[] peaklistFilenames) throws Exception {
        // pre-scan all files
        System.out.print("Pre-scanning " + peaklistFilenames.length + " input files...");
        long start = System.currentTimeMillis();

        IPeaklistScanner fileScanner = new PeakListFileScanner();
        Map<Integer, List<SpectrumReference>> loadedSpectrumReferenceMap = fileScanner.getSpectraPerMajorPeaks(peaklistFilenames, 5);

        printDone(start);

        // write out the spectra, one file per major peak
        System.out.print("Converting spectra to binary format per major peak");
        start = System.currentTimeMillis();


        for (int majorPeak : loadedSpectrumReferenceMap.keySet()) {
            File outputFile = getMajorPeakSourceFile(majorPeak, tmpSpectraPerPeakDir);
            SpectrumWriter spectrumWriter = new SpectrumWriter(outputFile);
            spectrumWriter.writeSpectra(loadedSpectrumReferenceMap.get(majorPeak), peaklistFilenames);

            System.out.print(".");
        }

        printDone(start);

        return loadedSpectrumReferenceMap.keySet();
    }

    private static File getMajorPeakSourceFile(int majorPeak, File dir) {
        File outputFile = new File(dir, "peak_" + majorPeak + ".cls");

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
