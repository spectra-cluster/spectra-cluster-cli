package uk.ac.ebi.pride.spectracluster.spectra_list;

/**
 * References a spectrum in an external file
 * Created by jg on 13.05.15.
 */
public final class SpectrumReference implements Comparable<SpectrumReference> {
    private final int fileId;
    private final int spectrumIndex;
    private final float precursorMz;

    public SpectrumReference(int fileId, int spectrumIndex, float precursorMz) {
        this.fileId = fileId;
        this.spectrumIndex = spectrumIndex;
        this.precursorMz = precursorMz;
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
