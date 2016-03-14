package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;

import java.io.File;

/**
 * Created by jg on 14.03.16.
 */
public class ClusteringJobReference {
    private final File inputFile;
    private final BinaryClusterFileReference outputFile;

    public ClusteringJobReference(File inputFile, BinaryClusterFileReference outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
    }

    public File getInputFile() {
        return inputFile;
    }

    public BinaryClusterFileReference getOutputFile() {
        return outputFile;
    }
}
