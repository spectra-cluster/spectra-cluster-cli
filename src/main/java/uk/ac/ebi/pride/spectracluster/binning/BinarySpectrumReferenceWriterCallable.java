package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.indexing.ClusteringFileIndex;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

/**
 * Created by jg on 28.05.15.
 */
public class BinarySpectrumReferenceWriterCallable implements Callable<BinaryClusterFileReference> {
    private final List<String> peaklistFiles;
    private final List<String> clusteringFiles;
    private final List<List<IndexElement>> fileIndices;
    private final List<ClusteringFileIndex> clusteringFileIndices;
    private final List<SpectrumReference> spectrumReferencesToWrite;
    private final File outputFile;
    private final boolean fastMode;

    public BinarySpectrumReferenceWriterCallable(List<String> peaklistFiles, List<String> clusteringFiles, List<List<IndexElement>> fileIndices, List<ClusteringFileIndex> clusteringFileIndices, List<SpectrumReference> spectrumReferencesToWrite, File outputFile, boolean fastMode) {
        this.peaklistFiles = peaklistFiles;
        this.clusteringFiles = clusteringFiles;
        this.fileIndices = fileIndices;
        this.clusteringFileIndices = clusteringFileIndices;
        this.spectrumReferencesToWrite = spectrumReferencesToWrite;
        this.outputFile = outputFile;
        this.fastMode = fastMode;
    }

    @Override
    public BinaryClusterFileReference call() throws Exception {
        ISpectrumReferenceWriter writer = new BinarySpectrumReferenceWriter(peaklistFiles, clusteringFiles,
                fileIndices, clusteringFileIndices, fastMode);

        writer.writeSpectra(spectrumReferencesToWrite, outputFile, peaklistFiles);

        // get the min and max m/z
        double minMz = Double.MAX_VALUE, maxMz = 0;

        for (SpectrumReference spectrumReference : spectrumReferencesToWrite) {
            if (spectrumReference.getPrecursorMz() < minMz) {
                minMz = spectrumReference.getPrecursorMz();
            }
            if (spectrumReference.getPrecursorMz() > maxMz) {
                maxMz = spectrumReference.getPrecursorMz();
            }
        }

        return new BinaryClusterFileReference(outputFile, minMz, maxMz, spectrumReferencesToWrite.size());
    }
}
