package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.io.File;
import java.io.FileDescriptor;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by jg on 13.05.15.
 */
public class SpectrumWriter {
    private final File outputFile;

    public SpectrumWriter(File outputFile) {
        this.outputFile = outputFile;
    }

    public void writeSpectra(List<SpectrumReference> spectrumReferences, String[] peakListFilenames) throws Exception {
        // initialize storage for various readers
        Map<Integer, JMzReader> readerPerFileIndex = new HashMap<Integer, JMzReader>();

        // sort the references
        Collections.sort(spectrumReferences);

        // open the file and write to it
        FileOutputStream outputStream = new FileOutputStream(outputFile);


    }
}
