package uk.ac.ebi.pride.spectracluster.conversion;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterIterable;
import uk.ac.ebi.pride.spectracluster.io.CGFClusterAppender;
import uk.ac.ebi.pride.spectracluster.spectra_list.ClusterReference;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by jg on 30.05.15.
 */
public class MergingCGFConverter implements IBinaryClusteringResultListener {
    private final boolean deleteTemporaryClusteringFiles;
    private final File resultFile;

    private List<ClusterReference> clusterReferences;

    /**
     * Creates a new MergingCGFConverter instance
     * @param resultFile The result file to write to.
     * @param deleteTemporaryClusteringFiles whether to delete clustering result files.
     */
    public MergingCGFConverter(File resultFile, boolean deleteTemporaryClusteringFiles) {
        this.resultFile = resultFile;
        this.deleteTemporaryClusteringFiles = deleteTemporaryClusteringFiles;

        this.clusterReferences = new ArrayList<ClusterReference>();
    }

    @Override
    public void onNewResultFile(BinaryClusterFileReference binaryClusteringResultFile) {
        try {
            // open the result file to write to
            FileOutputStream outputStream = new FileOutputStream(resultFile, true);

            // read from the clustering result file
            ObjectInputStream objectInputStream = new ObjectInputStream(new FileInputStream(
                    binaryClusteringResultFile.getResultFile()));
            BinaryClusterIterable binaryClusterIterable = new BinaryClusterIterable(objectInputStream);

            for (ICluster cluster : binaryClusterIterable) {
                long offset = outputStream.getChannel().position();
                BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream));
                CGFClusterAppender.INSTANCE.appendCluster(writer, cluster);
                writer.flush();

                // save the position of the cluster
                ClusterReference clusterReference = new ClusterReference(
                        0, offset, cluster.getPrecursorMz(), cluster.getClusteredSpectraCount(), cluster.getId());
                clusterReferences.add(clusterReference);
            }

            // close the output file
            outputStream.close();

            // delete the files
            deleteTemporaryFiles(binaryClusteringResultFile.getResultFile());
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private void deleteTemporaryFiles(File binaryClusteringResultFile) {
        if (deleteTemporaryClusteringFiles) {
            binaryClusteringResultFile.delete();
        }
    }

    public List<ClusterReference> getClusterReferences() {
        return Collections.unmodifiableList(clusterReferences);
    }
}
