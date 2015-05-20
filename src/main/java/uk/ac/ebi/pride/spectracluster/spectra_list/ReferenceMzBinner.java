package uk.ac.ebi.pride.spectracluster.spectra_list;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Groups spectrum references in bins
 *
 * Created by jg on 20.05.15.
 */
public class ReferenceMzBinner {
    public static final int[] WINDOW_SIZES = {25, 10, 4};
    public static final int MAXIMAL_SPECTRA = 500000;

    public List<List<SpectrumReference>> groupSpectrumReferences(List<SpectrumReference> spectrumReferences) {
        // sort according to m/z
        Collections.sort(spectrumReferences);

        // start with single group
        List<List<SpectrumReference>> groupedSpectrumReferences = new ArrayList<List<SpectrumReference>>();
        groupedSpectrumReferences.add(spectrumReferences);

        for (int i = 0; i < WINDOW_SIZES.length; i++) {
            int windowSize = WINDOW_SIZES[i];
            int maximalSpectra = MAXIMAL_SPECTRA;

            // force a re-distribution in the first round
            if (i == 0)
                maximalSpectra = 1;

            groupedSpectrumReferences = groupReferences(groupedSpectrumReferences, windowSize, maximalSpectra);
        }

        return groupedSpectrumReferences;
    }

    private List<List<SpectrumReference>> groupReferences(List<List<SpectrumReference>> spectrumReferences, int windowSize, int maximalSpectra) {
        List<List<SpectrumReference>> regroupedSpectrumReferences = new ArrayList<List<SpectrumReference>>();

        for (List<SpectrumReference> spectrumReferenceGroup : spectrumReferences) {
            // if the group is small enough, just leave it
            if (spectrumReferenceGroup.size() <= maximalSpectra) {
                regroupedSpectrumReferences.add(spectrumReferenceGroup);
                continue;
            }

            List<List<SpectrumReference>> splitGroup = splitReferenceGroup(spectrumReferenceGroup, windowSize);

            regroupedSpectrumReferences.addAll(splitGroup);
        }

        return regroupedSpectrumReferences;
    }

    private List<List<SpectrumReference>> splitReferenceGroup(List<SpectrumReference> spectrumReferenceGroup, int windowSize) {
        int startMz = (int) Math.floor(spectrumReferenceGroup.get(0).getPrecursorMz());
        int maxMz = (int) Math.ceil(spectrumReferenceGroup.get(spectrumReferenceGroup.size() - 1).getPrecursorMz());

        List<List<SpectrumReference>> splitGroups = new ArrayList<List<SpectrumReference>>();
        int lastIndex = 0;

        for (int currentMinMz = startMz; currentMinMz < maxMz; currentMinMz = currentMinMz + windowSize) {
            int currentMaxMz = currentMinMz + windowSize;

            List<SpectrumReference> currentGroup = new ArrayList<SpectrumReference>();

            for (int i = lastIndex; i < spectrumReferenceGroup.size(); i++) {
                SpectrumReference spectrumReference = spectrumReferenceGroup.get(i);

                if (spectrumReference.getPrecursorMz() < currentMinMz) {
                    lastIndex = i;
                    continue;
                }

                if (spectrumReference.getPrecursorMz() >= currentMaxMz) {
                    break;
                }

                currentGroup.add(spectrumReference);
                lastIndex = i + 1; // start with next item
            }

            splitGroups.add(currentGroup);
        }

        return splitGroups;
    }
}
