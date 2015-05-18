package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.File;
import java.util.*;

/**
 * Created by jg on 13.05.15.
 */
public class PeakListFileScanner implements IPeaklistScanner {
    private List<List<IndexElement>> fileIndices = new ArrayList<List<IndexElement>>();

    @Override
    public Map<Integer, List<SpectrumReference>> getSpectraPerMajorPeaks(String[] filenames, int nMajorPeaks) throws Exception {
        // map each file
        Map<Integer, List<SpectrumReference>> spectraPerMajorPeak = new HashMap<Integer, List<SpectrumReference>>();

        for (int i = 0; i < filenames.length; i++) {
            List<IndexElement> fileIndex = scanFile(filenames[i], i, spectraPerMajorPeak, nMajorPeaks);
            fileIndices.add(fileIndex);
        }

        return spectraPerMajorPeak;
    }

    private List<IndexElement> scanFile(String filename, int fileIndex, Map<Integer, List<SpectrumReference>> spectraPerMajorPeak, int nMajorPeaks) throws Exception {
        MgfFile reader = openFile(filename);

        Iterator<Spectrum> spectrumIterator = reader.getSpectrumIterator();

        int spectrumIndex = 0;

        while(spectrumIterator.hasNext()) {
            Spectrum readSpectrum = spectrumIterator.next();
            ISpectrum spectrum = SpectrumConverter.convertJmzReaderSpectrum(readSpectrum, readSpectrum.getId()); // id is not used

            // spectrum index is 1-based
            SpectrumReference spectrumReference = new SpectrumReference(fileIndex, spectrumIndex + 1, spectrum.getPrecursorMz());

            int[] majorPeaks = spectrum.asMajorPeakMZs(nMajorPeaks);

            for (int majorPeak : majorPeaks) {
                if (!spectraPerMajorPeak.containsKey(majorPeak))
                    spectraPerMajorPeak.put(majorPeak, new ArrayList<SpectrumReference>());

                spectraPerMajorPeak.get(majorPeak).add(spectrumReference);
            }

            spectrumIndex++;
        }

        return reader.getIndex();
    }

    private MgfFile openFile(String filename) throws Exception {
        if (filename.toLowerCase().endsWith(".mgf"))
            return new MgfFile(new File(filename));

        throw new Exception("Unsupported filetype encountered: " + filename);
    }

    public List<List<IndexElement>> getFileIndices() {
        return Collections.unmodifiableList(fileIndices);
    }
}
