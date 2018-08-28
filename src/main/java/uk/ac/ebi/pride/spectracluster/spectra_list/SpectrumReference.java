package uk.ac.ebi.pride.spectracluster.spectra_list;

import java.util.UUID;

/**
 * References a spectrum in an external file
 * Created by jg on 13.05.15.
 */
public final class SpectrumReference implements Comparable<SpectrumReference> {
    private final int fileId;
    private final int spectrumIndex;
    private final float precursorMz;
    private final String spectrumId;
    /**
     * Indicates that the current spectrum is referencing a cluster.
     */
    public static final int IS_CLUSTER = -1;

    public SpectrumReference(int fileId, int spectrumIndex, float precursorMz) {
        this(fileId, spectrumIndex, precursorMz, UUID.randomUUID().toString());
    }

    public SpectrumReference(int fileId, int spectrumIndex, float precursorMz, String spectrumId) {
        this.fileId = fileId;
        this.spectrumIndex = spectrumIndex;
        this.precursorMz = precursorMz;
        this.spectrumId = spectrumId;
    }

    public int getFileId() {
        return fileId;
    }

    public int getSpectrumIndex() {
        return spectrumIndex;
    }

    public float getPrecursorMz() {
        return precursorMz;
    }

    public String getSpectrumId() {
        return spectrumId;
    }

    @Override
    public int compareTo(SpectrumReference o) {
        return Float.compare(this.precursorMz, o.precursorMz);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SpectrumReference that = (SpectrumReference) o;

        if (fileId != that.fileId) return false;
        if (spectrumIndex != that.spectrumIndex) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = fileId;
        result = 31 * result + spectrumIndex;
        return result;
    }
}
