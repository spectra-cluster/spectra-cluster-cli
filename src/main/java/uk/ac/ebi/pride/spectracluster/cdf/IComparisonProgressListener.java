package uk.ac.ebi.pride.spectracluster.cdf;

/**
 * Created by jg on 07.04.16.
 */
public interface IComparisonProgressListener {
    public void progress(int completedCalculations, int totalCalculations);
}
