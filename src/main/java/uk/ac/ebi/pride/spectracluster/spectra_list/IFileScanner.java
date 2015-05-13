package uk.ac.ebi.pride.spectracluster.spectra_list;

import java.util.List;
import java.util.Map;

/**
 * Created by jg on 13.05.15.
 */
public interface IFileScanner {
    public Map<Integer, List<SpectrumReference>> getSpectraPerMajorPeaks(String[] filenames, int nMajorPeaks) throws Exception;
}
