package uk.ac.ebi.pride.spectracluster.clustering;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;

import java.util.Comparator;

/**
 * Created by jg on 15.05.15.
 */
public class ClusterMzComparator implements Comparator<ICluster> {

    public final static ClusterMzComparator INSTANCE = new ClusterMzComparator();

    protected ClusterMzComparator() {

    }

    @Override
    public int compare(ICluster o1, ICluster o2) {
        return Float.compare(o1.getPrecursorMz(), o2.getPrecursorMz());
    }
}
