package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

/**
 * Created by jg on 28.05.15.
 */
public class BinarySpectrumReferenceWriterCallable implements Callable<BinaryClusterFileReference> {
    private final String[] peaklistFiles;
    private final List<List<IndexElement>> fileIndices;
    private final List<SpectrumReference> spectrumReferencesToWrite;
    private final File outputFile;
    private final boolean fastMode;

    public BinarySpectrumReferenceWriterCallable(String[] peaklistFiles, List<List<IndexElement>> fileIndices, List<SpectrumReference> spectrumReferencesToWrite, File outputFile, boolean fastMode) {
        this.peaklistFiles = peaklistFiles;
        this.fileIndices = fileIndices;
        this.spectrumReferencesToWrite = spectrumReferencesToWrite;
        this.outputFile = outputFile;
        this.fastMode = fastMode;
    }

    @Override
    public BinaryClusterFileReference call() throws Exception {
        ISpectrumReferenceWriter writer = new BinarySpectrumReferenceWriter(peaklistFiles, fileIndices, fastMode);

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
