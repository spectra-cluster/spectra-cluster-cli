package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.util.List;
import java.util.Map;

/**
 * Created by jg on 13.05.15.
 */
public interface IPeaklistScanner {
    public Map<Integer, List<SpectrumReference>> getSpectraPerMajorPeaks(String[] filenames, int nMajorPeaks) throws Exception;

    public List<SpectrumReference> getSpectrumReferences(String[] filenames) throws Exception;

    public List<List<IndexElement>> getFileIndices();
}
