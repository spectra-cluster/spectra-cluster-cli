package uk.ac.ebi.pride.spectracluster.spectra_list;

/**
 * Created by jg on 16.05.15.
 */
public final class ClusterReference implements Comparable<ClusterReference> {
    private final int fileId;
    private final long offset;
    private final float precursorMz;
    private final int size;

    public ClusterReference(int fileId, long offset, float precursorMz, int size) {
        this.fileId = fileId;
        this.offset = offset;
        this.precursorMz = precursorMz;
        this.size = size;
    }

    public int getFileId() {
        return fileId;
    }

    public long getOffset() {
        return offset;
    }

    public float getPrecursorMz() {
        return precursorMz;
    }

    public int getSize() {
        return size;
    }

    /**
     * Sort according to precursor m/z in ascending order and size in
     * descending order.
     * @param o
     * @return
     */
    @Override
    public int compareTo(ClusterReference o) {
        if (this.getPrecursorMz() == o.getPrecursorMz())
            return Integer.compare(o.getSize(), this.getSize());

        return Float.compare(this.getPrecursorMz(), o.getPrecursorMz());
    }
}
