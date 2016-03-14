package uk.ac.ebi.pride.spectracluster.clustering;

import org.apache.commons.io.FileUtils;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.engine.GreedyIncrementalClusteringEngine;
import uk.ac.ebi.pride.spectracluster.engine.IIncrementalClusteringEngine;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterIterable;
import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.util.ClusteringJobReference;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.FractionTICPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison.ClusterShareMajorPeakPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison.IsKnownComparisonsPredicate;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * Created by jg on 15.05.15.
 */
public class BinaryFileClusteringCallable implements Callable<ClusteringJobReference> {
    public static final ISimilarityChecker SIMILARITY_CHECKER = new CombinedFisherIntensityTest(Defaults.getFragmentIonTolerance());
    public static final int DEFAULT_MAJOR_PEAK_COUNT = 5;

    public final IFunction<List<IPeak>, List<IPeak>> peakFilterFunction;

    private final File outputFile;
    private final File inputFile;
    private final List<Float> thresholds;

    private final float minMz;
    private final float maxMz;

    public BinaryFileClusteringCallable(File outputFile, File inputFile, List<Float> thresholds, boolean fastMode) {
        this(outputFile, inputFile, thresholds, fastMode, -1, -1);
    }

    /**
     * Create a new BinaryFileClustering object.
     * @param outputFile File to write the output to.
     * @param inputFile The (binary) file to cluster.
     * @param thresholds The thresholds to use in the clustering rounds.
     * @param fastMode If set to true no peak filter will be applied. Otherwise, each spectrum will be filtered for the comparison.
     * @param minMz All clusters below the set m/z will be ignored and simply written to the output file.
     * @param maxMz All clusters above the set m/z will be ignored and simply written to the output file. If set to -1
     *              the value is ignored.
     */
    public BinaryFileClusteringCallable(File outputFile, File inputFile, List<Float> thresholds, boolean fastMode, float minMz, float maxMz) {
        this.minMz = minMz;
        this.maxMz = maxMz;
        this.outputFile = outputFile;
        this.inputFile = inputFile;
        this.thresholds = thresholds;

        if (fastMode) {
            peakFilterFunction = null;
        }
        else {
            peakFilterFunction = new FractionTICPeakFunction(0.5f, 20);
        }
    }

    @Override
    public ClusteringJobReference call() throws Exception {
        try {
            File currentInputFile = inputFile;
            long start = System.currentTimeMillis();
            int nSpectra = 0;
            float minMz = Float.MAX_VALUE, maxMz = 0;

            for (int nRound = 0; nRound < thresholds.size(); nRound++) {
                float threshold = thresholds.get(nRound);

                IComparisonPredicate<ICluster> comparisonPredicate;

                if (nRound == 0) {
                    // first round only compare spectra that share a major peak
                    comparisonPredicate = new ClusterShareMajorPeakPredicate(DEFAULT_MAJOR_PEAK_COUNT);
                } else {
                    // subsequent rounds only compare known matches
                    comparisonPredicate = new IsKnownComparisonsPredicate();
                }

                IIncrementalClusteringEngine incrementalClusteringEngine =
                        createIncrementalClusteringEngine(threshold, comparisonPredicate);

                // create the result file
                File tmpOutputfile = File.createTempFile("clustering_tmp", ".cls");
                ObjectOutputStream outputStream = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(tmpOutputfile)));

                // read the clusters
                ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(currentInputFile)));
                BinaryClusterIterable clusterIterable = new BinaryClusterIterable(inputStream);

                // do the actual clustering
                for (ICluster clusterToAdd : clusterIterable) {
                    if (nRound == 0) {
                        nSpectra++;
                        if (clusterToAdd.getPrecursorMz() < minMz)
                            minMz = clusterToAdd.getPrecursorMz();
                        if (clusterToAdd.getPrecursorMz() > maxMz)
                            maxMz = clusterToAdd.getPrecursorMz();
                    }

                    // write out clusters that are below of above the set m/z limit
                    if (clusterToAdd.getPrecursorMz() < minMz ||
                            (maxMz > -1 && clusterToAdd.getPrecursorMz() > maxMz)) {
                        List<ICluster> clusterList = new ArrayList<ICluster>(1);
                        clusterList.add(clusterToAdd);
                        writeOutClusters(clusterList, outputStream);
                        continue;
                    }

                    // do the clustering
                    Collection<ICluster> removedClusters = incrementalClusteringEngine.addClusterIncremental(clusterToAdd);

                    // write out the removed clusters
                    if (!removedClusters.isEmpty()) {
                        writeOutClusters(removedClusters, outputStream);
                    }
                }

                // write out the final clusters
                Collection<ICluster> clusters = incrementalClusteringEngine.getClusters();
                writeOutClusters(clusters, outputStream);

                // close the output file
                BinaryClusterAppender.INSTANCE.appendEnd(outputStream);
                outputStream.close();

                // close the input file
                inputStream.close();

                // copy the output file
                if (outputFile.exists()) {
                    if (!outputFile.delete())
                        throw new Exception("Failed to delete file " + outputFile.toString());
                }

                FileUtils.copyFile(tmpOutputfile, outputFile);

                if (!tmpOutputfile.delete()) {
                    throw new Exception("Failed to delete temporary file");
                }

                // read from the last output the next time
                currentInputFile = outputFile;
            }

            printCompletion(inputFile.getName(), start, nSpectra, minMz, maxMz);

            return new ClusteringJobReference(inputFile, new BinaryClusterFileReference(outputFile, minMz, maxMz));
        } catch (Exception e) {
            System.out.println("Error: " + e.getMessage());
            e.printStackTrace();
            throw (e);
        }
    }

    private void printCompletion(String filename, long start, int nSpectra, float minMz, float maxMz) {
        System.out.println(String.format("Processed %s in %.2f min (%d spectra, %.2f m/z - %.2f m/z)",
                filename, (float) (System.currentTimeMillis() - start) / 1000 / 60, nSpectra, minMz, maxMz));
    }

    private IIncrementalClusteringEngine createIncrementalClusteringEngine(double clusteringPrecision, IComparisonPredicate<ICluster> comparisonPredicate) {
        IIncrementalClusteringEngine clusteringEngine = new GreedyIncrementalClusteringEngine(
                SIMILARITY_CHECKER,
                Defaults.getDefaultSpectrumComparator(),
                Defaults.getDefaultPrecursorIonTolerance(),
                clusteringPrecision,
                peakFilterFunction,
                comparisonPredicate);

        return clusteringEngine;
    }

    protected void writeOutClusters(Collection<ICluster> clusters, ObjectOutputStream outputStream) {
        List<ICluster> sortedRemovedClusters = new ArrayList<ICluster>(clusters);
        Collections.sort(sortedRemovedClusters, ClusterMzComparator.INSTANCE);

        for (ICluster c : sortedRemovedClusters)
            BinaryClusterAppender.INSTANCE.appendCluster(outputStream, c);
    }
}
