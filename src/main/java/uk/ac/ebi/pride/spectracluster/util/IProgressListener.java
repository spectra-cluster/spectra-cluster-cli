package uk.ac.ebi.pride.spectracluster.util;

/**
 * Defines an interface for classes listening to progress
 * updates.
 *
 * Created by jg on 10.06.16.
 */
public interface IProgressListener {
    public void onProgressUpdate(ProgressUpdate progressUpdate);
}
