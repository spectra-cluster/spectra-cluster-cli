package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;

import java.io.File;
import java.util.List;

/**
 * Created by jg on 30.05.15.
 */
public interface ISpectrumReferenceWriter {
    public void writeSpectra(List<SpectrumReference> spectrumReferences, File outputFile, List<String> peakListFilenames) throws Exception;
}
