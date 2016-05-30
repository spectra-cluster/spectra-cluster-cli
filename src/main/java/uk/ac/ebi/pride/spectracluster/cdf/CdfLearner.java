package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.spectra_list.IPeaklistScanner;
import uk.ac.ebi.pride.spectracluster.spectra_list.ParsingMgfScanner;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * This class is used to learn the CDF based on
 * the passed spectra.
 * Created by jg on 07.04.16.
 */
public class CdfLearner {
    /**
     * Default number of comparisons to perform
     */
    public final static int DEFAULT_NUMBER_OF_COMARPSISONS = 1000000;
    /**
     * Minimum difference in m/z for spectra to certainly originate
     * from different peptides.
     */
    public final double MIN_MZ_DIFFERENCE = 4.0;

    /**
     * Number of comparisons to perform per thread
     */
    private final static int COMPARISONS_PER_THREAD = 10000;

    private IPeaklistScanner peaklistScanner = new ParsingMgfScanner();
    private int numberOfComparisons = DEFAULT_NUMBER_OF_COMARPSISONS;
    private List<IComparisonProgressListener> listeners = new ArrayList<IComparisonProgressListener>();

    /**
     * Derives the cumulative distribution function based on the
     * passed peak list files.
     * @param peaklistFilenames Array of peaklist filenames.
     * @return The derived CumulativeDistributionFunction object.
     * @throws Exception
     */
    public CdfResult learnCumulativeDistribution(String[] peaklistFilenames, int nJobs) throws Exception {
        // prescan the peaklists
        List<SpectrumReference> spectrumReferences = peaklistScanner.getSpectrumReferences(peaklistFilenames);

        // make sure a sufficient number of spectra are available
        // to calculate the CDF using a 1.5 safety margin
        if (Math.sqrt(numberOfComparisons) > spectrumReferences.size() * 1.5) {
            throw new Exception("Insufficient number of spectra passed to derive a cumlative distribution function");
        }

        // create the defined number of random matches
        List<SpectrumMatch> matches = createRandomMatches(spectrumReferences);
        // sort the matches based on the lower spectrum index
        // Thereby it's more likely that spectra from the same
        // file are used for the comparison within one thread
        Collections.sort(matches);

        // launch the comparison jobs, 10000 per thread
        List<SpectrumMatch> currentMatches = new ArrayList<SpectrumMatch>();
        ExecutorService executorService = Executors.newFixedThreadPool(nJobs);
        List<Future<CdfResult>> cdfResultFutures = new ArrayList<Future<CdfResult>>();

        for (SpectrumMatch spectrumMatch : matches) {
            currentMatches.add(spectrumMatch);

            if (currentMatches.size() >= COMPARISONS_PER_THREAD) {
                CdfComparisonCallable cdfComparisonCallable = new CdfComparisonCallable(
                        currentMatches, spectrumReferences, peaklistFilenames, peaklistScanner.getFileIndices());
                Future<CdfResult> cdfResultFuture = executorService.submit(cdfComparisonCallable);
                cdfResultFutures.add(cdfResultFuture);

                currentMatches = new ArrayList<SpectrumMatch>();
            }
        }

        if (currentMatches.size() > 0) {
            CdfComparisonCallable cdfComparisonCallable = new CdfComparisonCallable(
                    currentMatches, spectrumReferences, peaklistFilenames, peaklistScanner.getFileIndices());
            Future<CdfResult> cdfResultFuture = executorService.submit(cdfComparisonCallable);
            cdfResultFutures.add(cdfResultFuture);
        }

        // show that we are done
        executorService.shutdown();

        // wait for everything to complete
        CdfResult mergedCdfResult = null;
        boolean threadsRunning = true;
        Set<Integer> completedJobIds = new HashSet<Integer>();
        int completed = 0;

        while (threadsRunning) {
            threadsRunning = false;
            for (int i = 0; i < cdfResultFutures.size(); i++) {
                // ignore all threads that have completed before
                if (completedJobIds.contains(i)) {
                    continue;
                }

                Future<CdfResult> cdfResultFuture = cdfResultFutures.get(i);

                // save new results
                if (cdfResultFuture.isDone()) {
                    // save the result
                    CdfResult cdfResult = cdfResultFuture.get();
                    completed += cdfResult.getTotalComparisons();
                    if (mergedCdfResult == null) {
                        mergedCdfResult = cdfResult;
                    } else {
                        mergedCdfResult.addCdfResult(cdfResult);
                    }

                    notifyListeners(completed, numberOfComparisons);

                    completedJobIds.add(i);
                }
                else {
                    // indicate that some threads are still running
                    threadsRunning = true;
                }
            }

            // wait a second for new threads to complete
            Thread.sleep(1000);
        }

        // only wait 1 second since everything should be done already
        executorService.awaitTermination(1, TimeUnit.SECONDS);

        // create the cumulative distribution function from the results
        return mergedCdfResult;
    }

    /**
     * Derives the cumulative distribution function based on the
     * passed peak list files.
     * @param peaklistFilenames Array of peaklist filenames.
     * @return The derived CumulativeDistributionFunction object.
     * @throws Exception
     */
    public CumulativeDistributionFunction learnCumulativeDistributionFunction(String[] peaklistFilenames, int nJobs) throws Exception {
        CdfResult cdfResult = learnCumulativeDistribution(peaklistFilenames, nJobs);

        return CumulativeDistributionFunction.fromString(cdfResult.toString());
    }

    /**
     * Create the required number of random spectrum pairs
     * to derive the cumulative distribution function from.
     * @param spectrumReferences
     * @return
     */
    private List<SpectrumMatch> createRandomMatches(List<SpectrumReference> spectrumReferences) {
        List<SpectrumMatch> matches = new ArrayList<SpectrumMatch>(numberOfComparisons);
        Random random = new Random();
        int max = spectrumReferences.size() - 1;

        while (matches.size() < numberOfComparisons) {
            int specIndex1 = random.nextInt(max + 1);
            int specIndex2 = random.nextInt(max + 1);

            double mz1 = spectrumReferences.get(specIndex1).getPrecursorMz();
            double mz2 = spectrumReferences.get(specIndex2).getPrecursorMz();

            double mzDifference = Math.abs(mz1 - mz2);

            if (mzDifference < MIN_MZ_DIFFERENCE) {
                continue;
            }

            matches.add(new SpectrumMatch(specIndex1, specIndex2));
        }

        return matches;
    }

    private void notifyListeners(int completed, int total) {
        for (IComparisonProgressListener listener : listeners) {
            listener.progress(completed, total);
        }
    }

    public void addListener(IComparisonProgressListener listener) {
        listeners.add(listener);
    }

    /**
     * The number of comparisons to perform when building
     * the cumulative distribution function.
     */
    public int getNumberOfComparisons() {
        return numberOfComparisons;
    }

    /**
     * Set the number of comparisons to perform when building
     * the cumulative distribution function.
     * @param numberOfComparisons
     */
    public void setNumberOfComparisons(int numberOfComparisons) {
        this.numberOfComparisons = numberOfComparisons;
    }

    /**
     * Stores a comparison to perform between to spectra
     */
    public class SpectrumMatch implements Comparable<SpectrumMatch> {
        private final int specIndex1;
        private final int specIndex2;

        public SpectrumMatch(int specIndex1, int specIndex2) {
            // make sure specIndex1 is always smaller
            if (specIndex1 < specIndex2) {
                this.specIndex1 = specIndex1;
                this.specIndex2 = specIndex2;
            }
            else {
                this.specIndex1 = specIndex2;
                this.specIndex2 = specIndex1;
            }
        }

        public int getSpecIndex1() {
            return specIndex1;
        }

        public int getSpecIndex2() {
            return specIndex2;
        }

        @Override
        public int compareTo(SpectrumMatch o) {
            return Integer.compare(this.getSpecIndex1(), o.getSpecIndex1());
        }
    }
}
