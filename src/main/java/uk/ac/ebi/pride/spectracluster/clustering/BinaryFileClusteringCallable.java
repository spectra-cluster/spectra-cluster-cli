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
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.FractionTICPeakFunction;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * Created by jg on 15.05.15.
 */
public class BinaryFileClusteringCallable implements Callable<File> {
    public static final float FRAGMENT_TOLERANCE = 0.5F;
    public static final ISimilarityChecker SIMILARITY_CHECKER = new CombinedFisherIntensityTest(FRAGMENT_TOLERANCE);
    public static final float WINDOW_SIZE = 4.0F;
    public static final double[] CLUSTERING_PRECISION = {0.9999F, 0.9998F, 0.9996F, 0.9994F, 0.9992F, 0.999F};
    public static final IFunction<List<IPeak>, List<IPeak>> peakFilterFunction = new FractionTICPeakFunction(0.5f, 20);
    public static final boolean ONLY_COMPARE_KNOWN_MATCHES = true;

    private final File outputFile;
    private final File inputFile;

    public BinaryFileClusteringCallable(File outputFile, File inputFile) {
        this.outputFile = outputFile;
        this.inputFile = inputFile;
    }

    @Override
    public File call() throws Exception {
        File currentInputFile = inputFile;
        long start = System.currentTimeMillis();

        for (double clusteringPrecision : CLUSTERING_PRECISION) {
            IIncrementalClusteringEngine incrementalClusteringEngine =
                    createIncrementalClusteringEngine(clusteringPrecision);

            // create the result file
            File tmpOutputfile = File.createTempFile("clustering_tmp", ".cls");
            ObjectOutputStream outputStream = new ObjectOutputStream(new FileOutputStream(tmpOutputfile));

            // read the clusters
            ObjectInputStream inputStream = new ObjectInputStream(new FileInputStream(currentInputFile));
            BinaryClusterIterable clusterIterable = new BinaryClusterIterable(inputStream);

            // do the actual clustering
            for (ICluster clusterToAdd : clusterIterable) {
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

        System.out.println("Processed " + inputFile.getName() + " in " +
                String.format("%.2f", (double) (System.currentTimeMillis() - start) / 1000 / 60) + " min");

        return outputFile;
    }

    private IIncrementalClusteringEngine createIncrementalClusteringEngine(double clusteringPrecision) {
        IIncrementalClusteringEngine clusteringEngine = new GreedyIncrementalClusteringEngine(
                SIMILARITY_CHECKER,
                Defaults.getDefaultSpectrumComparator(),
                WINDOW_SIZE,
                clusteringPrecision,
                peakFilterFunction,
                ONLY_COMPARE_KNOWN_MATCHES);

        return clusteringEngine;
    }

    protected void writeOutClusters(Collection<ICluster> clusters, ObjectOutputStream outputStream) {
        List<ICluster> sortedRemovedClusters = new ArrayList<ICluster>(clusters);
        Collections.sort(sortedRemovedClusters, ClusterMzComparator.INSTANCE);

        for (ICluster c : sortedRemovedClusters)
            BinaryClusterAppender.INSTANCE.appendCluster(outputStream, c);
    }
}
