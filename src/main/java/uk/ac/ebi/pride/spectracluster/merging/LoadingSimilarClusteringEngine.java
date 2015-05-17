package uk.ac.ebi.pride.spectracluster.merging;

import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.engine.SimilarClusterMergingEngine;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;

import java.util.Comparator;
import java.util.Set;

/**
 * This clustering engine also merges clusters if their share
 * a certain proportion of spectra. But if this is the case,
 * the peak lists are loaded again from the source file.
 *
 * Created by jg on 17.05.15.
 */
public class LoadingSimilarClusteringEngine extends SimilarClusterMergingEngine {
    public LoadingSimilarClusteringEngine(Comparator<ICluster> scm, float windowSize, double requiredSharedSpectra) {
        super(scm, windowSize, requiredSharedSpectra);
    }

    /**
     * this method is called by guaranteeClean to place any added clusters in play
     * for further clustering
     */
    protected void addToClusters(final ICluster clusterToAdd) {
        Set<String> spectraIdsToAdd = clusterToAdd.getSpectralIds();

        for (ICluster existingCluster : clusters) {
            Set<String> existingSpectraIds = existingCluster.getSpectralIds();
            double sharedSpectra = calculateSharedSpectra(spectraIdsToAdd, existingSpectraIds);

            if (sharedSpectra == 1) {
                // cluster is already stored
                return;
            }

            if (sharedSpectra >= requiredSharedSpectra) {
                if (clusterToAdd.storesPeakLists()) {
                    ISpectrum[] buffer = new ISpectrum[clusterToAdd.getClusteredSpectra().size()];
                    existingCluster.addSpectra(clusterToAdd.getClusteredSpectra().toArray(buffer));
                }
                else {
                    // it must be a greedy cluster
                    if (!GreedySpectralCluster.class.isInstance(existingCluster)) {
                        throw new IllegalStateException("Cannot add greedy cluster to non-greedy cluster");
                    }

                    // TODO: throw an exeception for now
                    if (true)
                        throw new UnsupportedOperationException("Loading peak lists not implemented yet.");
                    GreedySpectralCluster greedySpectralCluster = (GreedySpectralCluster) existingCluster;
                    greedySpectralCluster.addCluster(clusterToAdd);
                }
            }
        }

        // since the cluster wasn't merged, add it as new
        clusters.add(clusterToAdd);
    }
}
