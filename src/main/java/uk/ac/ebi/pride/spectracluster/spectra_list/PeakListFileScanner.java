package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.binning.BinarySpectrumReferenceWriter;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.File;
import java.util.*;

/**
 * Created by jg on 13.05.15.
 */
@Deprecated // this class was replaced by ParsingMgfScanner
public class PeakListFileScanner implements IPeaklistScanner {
    private List<List<IndexElement>> fileIndices = new ArrayList<List<IndexElement>>();

    @Override
    public List<SpectrumReference> getSpectrumReferences(String[] filenames) throws Exception {
        List<SpectrumReference> spectrumReferences = new ArrayList<SpectrumReference>();


        for (int i = 0; i < filenames.length; i++) {
            List<SpectrumReference> loadedSpectrumReferences = scanFile(filenames[i], i);
            spectrumReferences.addAll(loadedSpectrumReferences);
        }

        return spectrumReferences;
    }

    private List<SpectrumReference> scanFile(String filename, int fileIndex) throws Exception {
        MgfFile reader = openFile(filename);
        List<SpectrumReference> spectrumReferences = new ArrayList<SpectrumReference>();

        Iterator<Spectrum> spectrumIterator = reader.getSpectrumIterator();

        int spectrumIndex = 0;

        while(spectrumIterator.hasNext()) {
            Spectrum readSpectrum = spectrumIterator.next();
            ISpectrum spectrum = SpectrumConverter.convertJmzReaderSpectrum(readSpectrum, readSpectrum.getId()); // id is not used
            ISpectrum processedSpectrum = BinarySpectrumReferenceWriter.filterFunction.apply(spectrum);

            // spectrum index is 1-based
            SpectrumReference spectrumReference = new SpectrumReference(fileIndex, spectrumIndex + 1, processedSpectrum.getPrecursorMz());
            spectrumReferences.add(spectrumReference);

            spectrumIndex++;
        }

        if (fileIndex >= fileIndices.size())
            fileIndices.add(reader.getIndex());

        return spectrumReferences;
    }

    @Override
    public Map<Integer, List<SpectrumReference>> getSpectraPerMajorPeaks(String[] filenames, int nMajorPeaks) throws Exception {
        // map each file
        Map<Integer, List<SpectrumReference>> spectraPerMajorPeak = new HashMap<Integer, List<SpectrumReference>>();

        for (int i = 0; i < filenames.length; i++) {
            List<IndexElement> fileIndex = scanFilePerMajorPeak(filenames[i], i, spectraPerMajorPeak, nMajorPeaks);
            fileIndices.add(fileIndex);
        }

        return spectraPerMajorPeak;
    }

    private List<IndexElement> scanFilePerMajorPeak(String filename, int fileIndex, Map<Integer, List<SpectrumReference>> spectraPerMajorPeak, int nMajorPeaks) throws Exception {
        MgfFile reader = openFile(filename);

        Iterator<Spectrum> spectrumIterator = reader.getSpectrumIterator();

        int spectrumIndex = 0;

        while(spectrumIterator.hasNext()) {
            Spectrum readSpectrum = spectrumIterator.next();
            ISpectrum spectrum = SpectrumConverter.convertJmzReaderSpectrum(readSpectrum, readSpectrum.getId()); // id is not used
            ISpectrum processedSpectrum = BinarySpectrumReferenceWriter.filterFunction.apply(spectrum);

            // spectrum index is 1-based
            SpectrumReference spectrumReference = new SpectrumReference(fileIndex, spectrumIndex + 1, processedSpectrum.getPrecursorMz());

            int[] majorPeaks = processedSpectrum.asMajorPeakMZs(nMajorPeaks);

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
