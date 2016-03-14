package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;

/**
 * Created by jg on 14.03.16.
 */
public class ClusteringJobReference {
    private final BinaryClusterFileReference inputFile;
    private final BinaryClusterFileReference outputFile;

    public ClusteringJobReference(BinaryClusterFileReference inputFile, BinaryClusterFileReference outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
    }

    public BinaryClusterFileReference getInputFile() {
        return inputFile;
    }

    public BinaryClusterFileReference getOutputFile() {
        return outputFile;
    }
}
