package uk.ac.ebi.pride.spectracluster.spectra_list;

import java.util.Comparator;

/**
 * Compares SpectrumReferenceS based on the files
 * they refer to
 * Created by jg on 07.04.16.
 */
public class SpectrumReferenceFileComparator implements Comparator<SpectrumReference> {
    public final static SpectrumReferenceFileComparator INSTANCE = new SpectrumReferenceFileComparator();

    /**
     * Make sure only the single static instance is used
     */
    private SpectrumReferenceFileComparator() {

    }

    @Override
    public int compare(SpectrumReference o1, SpectrumReference o2) {
        // for spectra in the same file, compare based on index in the file
        if (o1.getFileId() == o2.getFileId()) {
            return Integer.compare(o1.getSpectrumIndex(), o2.getSpectrumIndex());
        }

        // otherwise sort based on the input file
        return Integer.compare(o1.getFileId(), o2.getFileId());
    }
}
