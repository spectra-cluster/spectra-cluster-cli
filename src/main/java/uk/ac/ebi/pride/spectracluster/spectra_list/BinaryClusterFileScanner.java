package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterParser;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 16.05.15.
 */
public class BinaryClusterFileScanner implements IClusterScanner {
    @Override
    public List<ClusterReference> getClusterReferences(File inputFile, int fileId) throws Exception {
        // read all clusters
        FileInputStream inputStream = new FileInputStream(inputFile);
        ObjectInputStream objectInputStream = new ObjectInputStream(inputStream);

        List<ClusterReference> clusterReferences = new ArrayList<ClusterReference>();

        long offset = 0;
        String lastLine = (String) objectInputStream.readObject();

        while (!lastLine.equals("END")) {
            ICluster currentCluster = BinaryClusterParser.INSTANCE.parseNextCluster(objectInputStream, lastLine);
            ClusterReference clusterReference =
                    new ClusterReference(fileId, offset, currentCluster.getPrecursorMz(),
                            currentCluster.getClusteredSpectraCount());

            clusterReferences.add(clusterReference);

            // read the next line
            offset = inputStream.getChannel().position();
            lastLine = (String) objectInputStream.readObject();
        }

        return clusterReferences;
    }
}
