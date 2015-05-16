package uk.ac.ebi.pride.spectracluster.spectra_list;

import java.io.File;
import java.util.List;

/**
 * Created by jg on 16.05.15.
 */
public interface IClusterScanner {
    public List<ClusterReference> getClusterReferences(File inputFile, int fileId) throws Exception;
}
