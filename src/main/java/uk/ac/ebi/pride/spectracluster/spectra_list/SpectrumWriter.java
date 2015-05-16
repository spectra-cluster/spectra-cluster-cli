package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.spectracluster.util.SpectrumUtilities;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.*;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Writes spectra stored as SpectrumReferences to a file.
 *
 * Created by jg on 13.05.15.
 */
public class SpectrumWriter {
    private final File outputFile;
    protected static final IFunction<ISpectrum, ISpectrum> filterFunction = Defaults.getDefaultPeakFilter();

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
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(outputStream);

        for (SpectrumReference spectrumReference : spectrumReferences) {
            int fileIndex = spectrumReference.getFileId();

            if (fileIndex >= peakListFilenames.length)
                throw new Exception("Invalid file id for spectrum reference");

            if (!readerPerFileIndex.containsKey(fileIndex))
                readerPerFileIndex.put(fileIndex, openFile(peakListFilenames[fileIndex]));

            JMzReader fileReader = readerPerFileIndex.get(fileIndex);

            // load the spectrum
            Spectrum spectrum = fileReader.getSpectrumByIndex(spectrumReference.getSpectrumIndex());

            // pre-process the spectrum
            ISpectrum processedSpectrum = filterFunction.apply(SpectrumConverter.convertJmzReaderSpectrum(spectrum));

            // write it to the file
            BinaryClusterAppender.INSTANCE.appendCluster(objectOutputStream, ClusterUtilities.asCluster(processedSpectrum));
        }

        // close the file
        BinaryClusterAppender.INSTANCE.appendEnd(objectOutputStream);
        objectOutputStream.close();
        outputStream.close();
    }

    private JMzReader openFile(String peakListFilename) throws Exception {
        if (peakListFilename.toLowerCase().endsWith(".mgf"))
            return new MgfFile(new File(peakListFilename));

        throw new Exception("Unknown file extension encountered: " + peakListFilename);
    }
}
