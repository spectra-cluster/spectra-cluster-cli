package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.clusteringfilereader.indexing.ClusteringFileIndex;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.indexing.ClusteringFileIndexer;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.indexing.ClusteringIndexElement;
import uk.ac.ebi.pride.tools.braf.BufferedRandomAccessFile;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.IndexElementImpl;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Scans a .clustering file and returns the clusters as a list of
 * SpectrumReferenceS.
 *
 * Created by jg on 28.05.15.
 */
public class ParsingClusteringScanner implements IPeaklistScanner {
    private List<ClusteringFileIndex> fileIndices;

    @Override
    public Map<Integer, List<SpectrumReference>> getSpectraPerMajorPeaks(String[] filenames, int nMajorPeaks) throws Exception {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<SpectrumReference> getSpectrumReferences(String[] filenames) throws Exception {
        fileIndices = new ArrayList<>();

        List<SpectrumReference> spectrumReferences = new ArrayList<SpectrumReference>();
        ClusteringFileIndexer indexer = new ClusteringFileIndexer();

        for (int i = 0; i < filenames.length; i++) {
            // index the .clustering file
            ClusteringFileIndex fileIndex = indexer.indexFile(new File(filenames[i]));

            // convert to a list of spectrum references
            for (ClusteringIndexElement indexElement : fileIndex.getIndex().values()) {
                SpectrumReference spectrumReference = new SpectrumReference(i, SpectrumReference.IS_CLUSTER, indexElement.getPrecursorMz(), indexElement.getId());
                spectrumReferences.add(spectrumReference);
            }

            fileIndices.add(fileIndex);
        }

        return spectrumReferences;
    }

    public List<List<IndexElement>> getFileIndices() {
        throw new UnsupportedOperationException();
    }

    public List<ClusteringFileIndex> getClusteringFileIndices() {
        return fileIndices;
    }
}
