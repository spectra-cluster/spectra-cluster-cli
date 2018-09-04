package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterIterable;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Writes the binary spectra / cluster from the input
 * files into new files where cluster from "neighbouring"
 * result files within the defined window size are are
 * written to the same file.
 *
 * Created by jg on 14.03.16.
 */
public class BinaryFileRebinner {
    /**
     * Writes the spectra / cluster from the input files into
     * new files where spectra within the windowSize of each
     * other from neighbouring files are written into the same
     * file. Each input file must span at least 2 x the
     * window Size.
     * @param inputFiles
     * @param outputDirectory
     * @param windowSize Size of the overlapping region IN EACH FILE.
     */
    public static List<BinaryClusterFileReference> rebinBinaryFiles(List<BinaryClusterFileReference> inputFiles,
                                                                    File outputDirectory, double windowSize)
            throws Exception {
        // test whether anything needs to be re-binned at all
        if (inputFiles.size() < 2) {
            return inputFiles;
        }

        // sort the input files according to m/z
        Collections.sort(inputFiles);

        // set the current max m/z
        double maxMz = inputFiles.get(0).getMaxMz() - windowSize;
        File outputFile = getResultFile(outputDirectory, inputFiles.get(0).getMinMz(), maxMz);
        ObjectOutputStream outputStream = new ObjectOutputStream(new BufferedOutputStream(
                new FileOutputStream(outputFile)));

        List<BinaryClusterFileReference> outputFiles = new ArrayList<BinaryClusterFileReference>();
        double outputMinMz = Double.MAX_VALUE, outputMaxMz = 0;
        int nCluster = 0;

        for (int i = 0; i < inputFiles.size(); i++) {
            BinaryClusterFileReference clusterFileReference = inputFiles.get(i);

            // open the file
            ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(
                    new FileInputStream(clusterFileReference.getResultFile())));
            BinaryClusterIterable clusterIterable = new BinaryClusterIterable(inputStream);

            // make sure the clusters are sorted according to m/z
            double lastClusterMz = 0;

            for (ICluster cluster : clusterIterable) {
                // if the cluster's precursor m/z is larger than the current maximum, create the next output file
                if (cluster.getPrecursorMz() > maxMz) {
                    BinaryClusterAppender.INSTANCE.appendEnd(outputStream);
                    outputStream.close();

                    // save the file reference
                    outputFiles.add(new BinaryClusterFileReference(outputFile, outputMinMz, outputMaxMz, nCluster));

                    // set the next max m/z
                    if (i < inputFiles.size() - 1) {
                        // use the next file's maximum
                        maxMz = inputFiles.get(i + 1).getMaxMz() - windowSize;
                    } else {
                        // simply use this file's
                        maxMz = clusterFileReference.getMaxMz();
                    }

                    // create the new file
                    outputFile = getResultFile(outputDirectory, outputMaxMz, maxMz);
                    outputStream = new ObjectOutputStream(new BufferedOutputStream(
                            new FileOutputStream(outputFile)));

                    // reset the lastClusterMz
                    lastClusterMz = 0;
                    // reset the statistics for the output file
                    outputMinMz = Double.MAX_VALUE;
                    outputMaxMz = 0;
                    nCluster = 0;
                }

                // make sure the clusters are sorted according to m/z
                if (lastClusterMz > cluster.getPrecursorMz() + 0.1 ) {
                    throw new Exception("Clusters are not sorted according to precursor m/z in " +
                            clusterFileReference.getResultFile().getName());
                }

                lastClusterMz = cluster.getPrecursorMz();

                // write the cluster
                BinaryClusterAppender.INSTANCE.appendCluster(outputStream, cluster);

                // update the file statistics
                if (cluster.getPrecursorMz() < outputMinMz) {
                    outputMinMz = cluster.getPrecursorMz();
                }
                if (cluster.getPrecursorMz() > outputMaxMz) {
                    outputMaxMz = cluster.getPrecursorMz();
                }
                nCluster++;
            }

            inputStream.close();
        }

        // close the currently written file and save it
        BinaryClusterAppender.INSTANCE.appendEnd(outputStream);
        outputStream.close();
        outputFiles.add(new BinaryClusterFileReference(outputFile, outputMinMz, outputMaxMz, nCluster));

        return outputFiles;
    }

    /**
     * Create a result file object for the current bin. Basically just
     * generates the output filename.
     * @param outputDirectory
     * @param minMz
     * @param maxMz
     * @return
     */
    private static File getResultFile(File outputDirectory, double minMz, double maxMz) {
        String filename = String.format("rebinnedSpectra_%04d_%04d.cls",
                (int) Math.floor(minMz),
                (int) Math.floor(maxMz));

        return new File(outputDirectory, filename);
    }
}
