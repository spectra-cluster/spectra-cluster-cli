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
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.util.ClusteringJobReference;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.predicate.ComparisonPredicates;
import uk.ac.ebi.pride.spectracluster.util.predicate.IComparisonPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.Predicates;
import uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison.ClusterPpmPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison.ClusterShareMajorPeakPredicate;
import uk.ac.ebi.pride.spectracluster.util.predicate.cluster_comparison.IsKnownComparisonsPredicate;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
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

    private final IPredicate<ICluster> clusterPredicate;

    private final File temporaryDirectory;

    public BinaryFileClusteringCallable(File outputFile, File inputFile, List<Float> thresholds, boolean fastMode, File temporaryDirectory) {
        this(outputFile, inputFile, thresholds, fastMode, -1, -1, temporaryDirectory, null);
    }

    public BinaryFileClusteringCallable(File outputFile, File inputFile, List<Float> thresholds, boolean fastMode, File temporaryDirectory, IPredicate<ICluster> clusterPredicate) {
        this(outputFile, inputFile, thresholds, fastMode, -1, -1, temporaryDirectory, clusterPredicate);
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
     * @param temporaryDirectory The directory where temporary clustering files should be stored
     * @param clusterPredicate If set, clusters that do not fulfill this predicate are being ignored.
     */
    public BinaryFileClusteringCallable(File outputFile, File inputFile, List<Float> thresholds, boolean fastMode, float minMz, float maxMz, File temporaryDirectory,
                                        IPredicate<ICluster> clusterPredicate) {
        this.minMz = minMz;
        this.maxMz = maxMz;
        this.outputFile = outputFile;
        this.inputFile = inputFile;
        this.thresholds = thresholds;
        this.temporaryDirectory = temporaryDirectory;
        this.clusterPredicate = clusterPredicate;

        if (fastMode) {
            peakFilterFunction = null;
        }
        else {
            peakFilterFunction = ClusteringSettings.getComparisonFilterFunction();
        }
    }

    @Override
    public ClusteringJobReference call() throws Exception {
        try {
            File currentInputFile = inputFile;
            long start = System.currentTimeMillis();
            int nSpectra = 0;
            float fileMinMz = Float.MAX_VALUE, fileMaxMz = 0;

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
                File tmpOutputfile = File.createTempFile("clustering_tmp", ".cls", temporaryDirectory);
                ObjectOutputStream outputStream = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(tmpOutputfile)));

                // read the clusters
                ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(currentInputFile)));
                BinaryClusterIterable clusterIterable = new BinaryClusterIterable(inputStream);

                // do the actual clustering
                for (ICluster clusterToAdd : clusterIterable) {
                    if (Thread.currentThread().isInterrupted()) {
                        outputStream.close();
                        inputStream.close();
                        throw new InterruptedException();
                    }

                    // ignore any cluster that does not fulfill the predicate
                    if (nRound == 0 && clusterPredicate != null && !clusterPredicate.apply(clusterToAdd)) {
                        continue;
                    }

                    if (nRound == 0) {
                        nSpectra++;
                    }

                    // update the file statistics in the last round
                    if (nRound == thresholds.size() - 1) {
                        if (clusterToAdd.getPrecursorMz() < fileMinMz)
                            fileMinMz = clusterToAdd.getPrecursorMz();
                        if (clusterToAdd.getPrecursorMz() > fileMaxMz)
                            fileMaxMz = clusterToAdd.getPrecursorMz();
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

                if (Thread.currentThread().isInterrupted()) {
                    outputStream.close();
                    inputStream.close();
                    throw new InterruptedException();
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

                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedException();
                }

                FileUtils.copyFile(tmpOutputfile, outputFile);

                if (!tmpOutputfile.delete()) {
                    throw new Exception("Failed to delete temporary file");
                }

                // read from the last output the next time
                currentInputFile = outputFile;
            }

            return new ClusteringJobReference(inputFile, new BinaryClusterFileReference(outputFile, fileMinMz, fileMaxMz, nSpectra));
        } catch (Exception e) {
            System.out.println("Error: " + e.getMessage());
            e.printStackTrace();
            throw (e);
        }
    }

    private IIncrementalClusteringEngine createIncrementalClusteringEngine(double clusteringPrecision, IComparisonPredicate<ICluster> comparisonPredicate) {
        // if the threshold is set in ppm, add a PPM predicate
        if (ClusteringSettings.ppmThreshold != null) {
            comparisonPredicate = ComparisonPredicates.and(comparisonPredicate,
                    new ClusterPpmPredicate(ClusteringSettings.ppmThreshold));
        }

        IIncrementalClusteringEngine clusteringEngine = new GreedyIncrementalClusteringEngine(
                SIMILARITY_CHECKER,
                Defaults.getDefaultSpectrumComparator(),
                Defaults.getDefaultPrecursorIonTolerance(),
                clusteringPrecision,
                peakFilterFunction,
                comparisonPredicate);

        return clusteringEngine;
    }

    protected void writeOutClusters(Collection<ICluster> clusters, ObjectOutputStream outputStream) throws InterruptedException {
        List<ICluster> sortedRemovedClusters = new ArrayList<ICluster>(clusters);
        Collections.sort(sortedRemovedClusters, ClusterMzComparator.INSTANCE);

        for (ICluster c : sortedRemovedClusters) {
            if (Thread.currentThread().isInterrupted())
                throw new InterruptedException();

            BinaryClusterAppender.INSTANCE.appendCluster(outputStream, c);
        }
    }
}
