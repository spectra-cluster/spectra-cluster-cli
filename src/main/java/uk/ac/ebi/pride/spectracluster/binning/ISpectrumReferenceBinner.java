package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;

import java.util.List;

/**
 * Created by jg on 30.05.15.
 */
public interface ISpectrumReferenceBinner {
    public List<List<SpectrumReference>> binSpectrumReferences(List<SpectrumReference> spectrumReferences);
}
