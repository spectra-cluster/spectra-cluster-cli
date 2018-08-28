package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;

import java.util.ArrayList;
import java.util.List;

/**
 * This ReferenceMzBinner makes sure that spectra are
 * always assigned to the same bin. It is only intended
 * for cases where spectra should be binned equally irrespective
 * of whether they are processed together or not.
 *
 * Created by jg on 23.12.17.
 */
public class FixedReferenceMzBinner implements ISpectrumReferenceBinner {
    private final float windowSize;

    /**
     * Create a new FixedReferenceMzBinner.
     * @param windowSize The window size in m/z
     */
    public FixedReferenceMzBinner(int windowSize) {
        this.windowSize = windowSize;
    }

    @Override
    public List<List<SpectrumReference>> binSpectrumReferences(List<SpectrumReference> spectrumReferences) {
        List<List<SpectrumReference>> binnedSpectrumReferences = new ArrayList<>();

        for (SpectrumReference specRef : spectrumReferences) {
            int bin = (int) Math.ceil(specRef.getPrecursorMz() / windowSize);

            // make sure the list for this bin already exists
            for (int i = binnedSpectrumReferences.size(); i <= bin; i++) {
                binnedSpectrumReferences.add(new ArrayList<>());
            }

            // add the spectrum reference
            binnedSpectrumReferences.get(bin).add(specRef);
        }

        return binnedSpectrumReferences;
    }
}
