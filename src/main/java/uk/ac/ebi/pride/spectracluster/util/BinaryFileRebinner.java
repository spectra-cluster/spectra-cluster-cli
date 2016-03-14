package uk.ac.ebi.pride.spectracluster.util;

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
        // make sure the file window sizes are large enough
        checkFileWindows(inputFiles, windowSize);

        // sort the input files according to m/z
        Collections.sort(inputFiles);

        int currentBin = 1;
        File outputFile = getResultFile(outputDirectory, currentBin);
        ObjectOutputStream outputStream = new ObjectOutputStream(new BufferedOutputStream(
                new FileOutputStream(outputFile)));

        List<BinaryClusterFileReference> outputFiles = new ArrayList<BinaryClusterFileReference>();
        double outputMinMz = Double.MAX_VALUE, outputMaxMz = 0;

        for (BinaryClusterFileReference clusterFileReference : inputFiles) {
            // set the current max m/z
            double maxMz = clusterFileReference.getMaxMz() - windowSize;

            // open the file
            ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(
                    new FileInputStream(clusterFileReference.getResultFile())));
            BinaryClusterIterable clusterIterable = new BinaryClusterIterable(inputStream);

            // make sure the clusters are sorted according to m/z
            double lastClusterMz = 0;

            for (ICluster cluster : clusterIterable) {
                // if the cluster's precursor m/z is larger than the current maximum, create the next output file
                if (cluster.getPrecursorMz() > maxMz) {
                    outputStream.close();
                    currentBin++;

                    // save the file reference
                    outputFiles.add(new BinaryClusterFileReference(outputFile, outputMinMz, outputMaxMz));

                    // create the new file
                    outputFile = getResultFile(outputDirectory, currentBin);
                    outputStream = new ObjectOutputStream(new BufferedOutputStream(
                            new FileOutputStream(outputFile)));

                    // reset the lastClusterMz
                    lastClusterMz = 0;
                    // reset the statistics for the output file
                    outputMinMz = Double.MAX_VALUE;
                    outputMaxMz = 0;
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
            }

            inputStream.close();
        }

        // close the currently written file and save it
        outputStream.close();
        outputFiles.add(new BinaryClusterFileReference(outputFile, outputMinMz, outputMaxMz));

        return outputFiles;
    }

    private static void checkFileWindows(List<BinaryClusterFileReference> inputFiles, double windowSize) throws Exception {
        // make sure the input files are "large enough"
        for (BinaryClusterFileReference clusterFileReference : inputFiles) {
            double fileWindowSize = clusterFileReference.getMaxMz() - clusterFileReference.getMinMz();

            if (fileWindowSize < windowSize * 2) {
                throw new Exception(clusterFileReference.getResultFile().getName() + " only spans " +
                        fileWindowSize + " m/z. Window size of " + windowSize + " is too large for re-binning.");
            }
        }
    }

    /**
     * Create a result file object for the current bin. Basically just
     * generates the output filename.
     * @param outputDirectory
     * @param currentBin
     * @return
     */
    private static File getResultFile(File outputDirectory, int currentBin) {
        String stringBin;

        if (currentBin >= 1000) {
            stringBin = String.valueOf(currentBin);
        }
        else if (currentBin >= 100) {
            stringBin = "0" + String.valueOf(currentBin);
        }
        else if (currentBin >= 10) {
            stringBin = "00" + String.valueOf(currentBin);
        }
        else {
            stringBin = "000" + String.valueOf(currentBin);
        }

        return new File(outputDirectory, "rebinned_bin_" + stringBin + ".cls");
    }
}
