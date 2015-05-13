package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.File;
import java.util.*;

/**
 * Created by jg on 13.05.15.
 */
public class PeakListFileScanner implements IFileScanner {
    private final int nMajorPeaks;

    public PeakListFileScanner(int nMajorPeaks) {
        this.nMajorPeaks = nMajorPeaks;
    }

    @Override
    public Map<Integer, List<SpectrumReference>> getSpectraPerMajorPeaks(String[] filenames, int nMajorPeaks) throws Exception {
        // map each file
        Map<Integer, List<SpectrumReference>> spectraPerMajorPeak = new HashMap<Integer, List<SpectrumReference>>();

        for (int i = 0; i < filenames.length; i++)
            scanFile(filenames[i], i, spectraPerMajorPeak);

        return spectraPerMajorPeak;
    }

    private void scanFile(String filename, int fileIndex, Map<Integer, List<SpectrumReference>> spectraPerMajorPeak) throws Exception {
        JMzReader reader = openFile(filename);

        Iterator<Spectrum> spectrumIterator = reader.getSpectrumIterator();

        int spectrumIndex = 0;

        while(spectrumIterator.hasNext()) {
            Spectrum readSpectrum = spectrumIterator.next();
            ISpectrum spectrum = SpectrumConverter.convertJmzReaderSpectrum(readSpectrum);

            SpectrumReference spectrumReference = new SpectrumReference(fileIndex, spectrumIndex, spectrum.getPrecursorMz());

            int[] majorPeaks = spectrum.asMajorPeakMZs(nMajorPeaks);

            for (int majorPeak : majorPeaks) {
                if (!spectraPerMajorPeak.containsKey(majorPeak))
                    spectraPerMajorPeak.put(majorPeak, new ArrayList<SpectrumReference>());

                spectraPerMajorPeak.get(majorPeak).add(spectrumReference);
            }
        }
    }

    private JMzReader openFile(String filename) throws Exception {
        if (filename.toLowerCase().endsWith(".mgf"))
            return new MgfFile(new File(filename));

        throw new Exception("Unsupported filetype encountered: " + filename);
    }
}
