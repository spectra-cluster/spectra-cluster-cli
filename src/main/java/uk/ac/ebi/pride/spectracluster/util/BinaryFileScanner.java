package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.cdf.SpectraPerBinNumberComparisonAssessor;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterIterable;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.util.predicate.IPredicate;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * This class holds methods to process a list of
 * binary clustering files and report the metadata
 * as ClusteringResult objects.
 *
 * Created by jg on 13.03.16.
 */
public class BinaryFileScanner {
    /**
     * Scans the passed binary clustering files and returns
     * the matching CLusteringResultS with the associated
     * metadata. Since this function is primarily dependent
     * on the disk I/O speed it should not run in parallel.
     * @param inputFiles
     * @return
     */
    public static List<BinaryClusterFileReference> scanBinaryFiles(File... inputFiles) throws IOException {
        return scanBinaryFiles(null, null, inputFiles);
    }

    /**
     * Scans the passed binary clustering files and returns
     * the matching CLusteringResultS with the associated
     * metadata. Since this function is primarily dependent
     * on the disk I/O speed it should not run in parallel.
     * @param spectraPerBinNumberComparisonAssessor If set, the found cluster are counted as spectra per bin.
     * @param clusterAddingPredicate Predicate which needs to be fulfilled for clusters to be considered. If set
     *                               to NULL, all clusteres are processed.
     * @param inputFiles
     * @return
     */
    public static List<BinaryClusterFileReference> scanBinaryFiles(SpectraPerBinNumberComparisonAssessor spectraPerBinNumberComparisonAssessor,
                                                                   IPredicate<ICluster> clusterAddingPredicate,
                                                                   File... inputFiles) throws IOException {
        List<BinaryClusterFileReference> binaryClusterFileReferences = new ArrayList<BinaryClusterFileReference>(inputFiles.length);

        for (File currentInputFile : inputFiles) {
            ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(currentInputFile)));
            BinaryClusterIterable clusterIterable = new BinaryClusterIterable(inputStream);

            // scan the file to get the min and max m/z
            double minMz = Double.MAX_VALUE, maxMz = 0;
            int nCluster = 0;

            for (ICluster cluster : clusterIterable) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedIOException();
                }

                if (clusterAddingPredicate != null && !clusterAddingPredicate.apply(cluster)) {
                    continue;
                }

                if (cluster.getPrecursorMz() < minMz) {
                    minMz = cluster.getPrecursorMz();
                }
                if (cluster.getPrecursorMz() > maxMz) {
                    maxMz = cluster.getPrecursorMz();
                }

                if (spectraPerBinNumberComparisonAssessor != null) {
                    spectraPerBinNumberComparisonAssessor.countSpectrum(cluster.getPrecursorMz());
                }

                nCluster++;
            }

            // save the file reference as a CLusteringResult object
            BinaryClusterFileReference binaryClusterFileReference = new BinaryClusterFileReference(
                    currentInputFile, minMz, maxMz, nCluster);
            binaryClusterFileReferences.add(binaryClusterFileReference);
        }

        return binaryClusterFileReferences;
    }
}
