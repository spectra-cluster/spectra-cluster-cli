package uk.ac.ebi.pride.spectracluster.merging;

import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.engine.SimilarClusterMergingEngine;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumWriter;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.File;
import java.util.*;

/**
 * This clustering engine also merges clusters if their share
 * a certain proportion of spectra. But if this is the case,
 * the peak lists are loaded again from the source file.
 * <p/>
 * Created by jg on 17.05.15.
 */
public class LoadingSimilarClusteringEngine extends SimilarClusterMergingEngine {
    private final Map<String, SpectrumReference> spectrumReferencesPerId;
    private final String[] peakListFiles;
    private Map<Integer, JMzReader> readersPerFileid = new HashMap<Integer, JMzReader>();
    private final List<List<IndexElement>> peakListIndices;

    public LoadingSimilarClusteringEngine(Comparator<ICluster> scm, float windowSize, double requiredSharedSpectra, Map<String, SpectrumReference> spectrumReferencesPerId, String[] peakListFiles, List<List<IndexElement>> peakListIndices) {
        super(scm, windowSize, requiredSharedSpectra);
        this.spectrumReferencesPerId = spectrumReferencesPerId;
        this.peakListFiles = peakListFiles;
        this.peakListIndices = peakListIndices;
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
                ISpectrum[] spectraToAdd = null;

                if (clusterToAdd.storesPeakLists()) {
                    spectraToAdd = new ISpectrum[clusterToAdd.getClusteredSpectra().size()];
                    spectraToAdd = clusterToAdd.getClusteredSpectra().toArray(spectraToAdd);
                } else {
                    // identify the actual spectra to load
                    Set<String> spectraToLoad = new HashSet<String>(clusterToAdd.getSpectralIds());
                    spectraToLoad.removeAll(existingSpectraIds);
                    // load the spectra from file
                    spectraToAdd = loadSpectraFromFile(spectraToLoad);
                }

                existingCluster.addSpectra(spectraToAdd);

                // keep the larger id
                if (clusterToAdd.getClusteredSpectraCount() > existingCluster.getClusteredSpectraCount())
                    existingCluster.setId(clusterToAdd.getId());

                return;
            }
        }

        // since the cluster wasn't merged, add it as new
        clusters.add(clusterToAdd);
    }

    private ISpectrum[] loadSpectraFromFile(Set<String> spectraToLoad) {
        try {
            List<ISpectrum> loadedSpectra = new ArrayList<ISpectrum>();

            for (String specId : spectraToLoad) {
                // get the spectrum reference
                SpectrumReference reference = spectrumReferencesPerId.get(specId);

                if (reference == null)
                    throw new IllegalStateException("Missing SpecrumReference for spectrum " + specId);

                // read the spectrum
                JMzReader reader;
                if (readersPerFileid.containsKey(reference.getFileId())) {
                    reader = readersPerFileid.get(reference.getFileId());
                } else {
                    // TODO: support multiple file types
                    reader = new MgfFile(new File(peakListFiles[reference.getFileId()]), peakListIndices.get(reference.getFileId()));
                    readersPerFileid.put(reference.getFileId(), reader);
                }

                ISpectrum loadedSpectrum = SpectrumConverter.convertJmzReaderSpectrum(reader.getSpectrumByIndex(reference.getSpectrumIndex()), reference.getSpectrumId());
                ISpectrum processedSpectrum = SpectrumWriter.filterFunction.apply(loadedSpectrum);
                // normalize the spectrum
                processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(
                        processedSpectrum, Defaults.getDefaultIntensityNormalizer().normalizePeaks(processedSpectrum.getPeaks()));
                loadedSpectra.add(processedSpectrum);
            }

            // convert to array and return
            ISpectrum[] returnValue = new ISpectrum[loadedSpectra.size()];
            return loadedSpectra.toArray(returnValue);
        } catch (Exception e) {
            throw new IllegalStateException(e);
        }
    }
}
