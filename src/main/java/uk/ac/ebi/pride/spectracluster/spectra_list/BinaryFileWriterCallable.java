package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

/**
 * Created by jg on 28.05.15.
 */
public class BinaryFileWriterCallable implements Callable<File> {
    private final String[] peaklistFiles;
    private final List<List<IndexElement>> fileIndices;
    private final List<SpectrumReference> spectrumReferencesToWrite;
    private final File outputFile;

    public BinaryFileWriterCallable(String[] peaklistFiles, List<List<IndexElement>> fileIndices, List<SpectrumReference> spectrumReferencesToWrite, File outputFile) {
        this.peaklistFiles = peaklistFiles;
        this.fileIndices = fileIndices;
        this.spectrumReferencesToWrite = spectrumReferencesToWrite;
        this.outputFile = outputFile;
    }

    @Override
    public File call() throws Exception {
        SpectrumWriter writer = new SpectrumWriter(peaklistFiles, fileIndices);

        writer.writeSpectra(spectrumReferencesToWrite, outputFile);

        return outputFile;
    }
}
