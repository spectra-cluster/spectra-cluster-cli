package uk.ac.ebi.pride.spectracluster.clustering;

import java.io.File;

/**
 * Created by jg on 18.05.15.
 */
public interface IBinaryClusteringResultListener {
    public void onNewResultFile(File binaryClusteringResultFile);
}
