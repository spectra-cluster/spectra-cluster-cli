package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.normalizer.IIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.function.Functions;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.FractionTICPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.HighestNPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveImpossiblyHighPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemovePrecursorPeaksFunction;

import java.util.List;

/**
 * This is a (static) class that holds certain user
 * defined settings.
 * Created by jg on 08.04.16.
 */
public class CliSettings {
    /**
     * This filter is applied to all spectra as soon as they are loaded from file.
     */
    public static final IFunction<ISpectrum, ISpectrum> DEFAULT_INITIAL_SPECTRUM_FILTER =  Functions.join(
            new RemoveImpossiblyHighPeaksFunction(),
            new RemovePrecursorPeaksFunction(Defaults.getFragmentIonTolerance()));

    private static IFunction<ISpectrum, ISpectrum> initialSpectrumFilter = DEFAULT_INITIAL_SPECTRUM_FILTER;

    /**
     * This filter is applied to all spectra when loaded from the file - after the initial
     * spectrum filter.
     */
    public static final IFunction<List<IPeak>, List<IPeak>> DEFAULT_LOADING_SPECTRUM_FILTER = new HighestNPeakFunction(70);
    private static IFunction<List<IPeak>, List<IPeak>> loadingSpectrumFilter = DEFAULT_LOADING_SPECTRUM_FILTER;

    /**
     * The intensity normalizer used
     */
    public static final IIntensityNormalizer DEFAULT_INTENTISTY_NORMALIZER = Defaults.getDefaultIntensityNormalizer();
    private static IIntensityNormalizer intensityNormalizer = DEFAULT_INTENTISTY_NORMALIZER;

    /**
     * The filter applied to every spectrum when performing the actual comparison. This
     * filter does not affect the spectra from which the consensus spectrum is build.
     */
    public static final IFunction<List<IPeak>, List<IPeak>> DEFAULT_COMPARISON_FILTER_FUNCTION = new FractionTICPeakFunction(0.5f, 20);
    private static IFunction<List<IPeak>, List<IPeak>> comparisonFilterFunction = DEFAULT_COMPARISON_FILTER_FUNCTION;

    /**
     * The initial filter applied to every spectrum when loaded from file.
     * @return
     */
    public static IFunction<ISpectrum, ISpectrum> getInitialSpectrumFilter() {
        return initialSpectrumFilter;
    }

    /**
     * Set the filter that is applied to every spectrum as soon as its loaded from
     * the input file.
     * @param initialSpectrumFilter
     */
    public static void setInitialSpectrumFilter(IFunction<ISpectrum, ISpectrum> initialSpectrumFilter) {
        CliSettings.initialSpectrumFilter = initialSpectrumFilter;
    }

    /**
     * The filter applied to every spectrum when performing the actual comparison. This
     * filter does not affect the spectra from which the consensus spectrum is build.
     * @return
     */
    public static IFunction<List<IPeak>, List<IPeak>> getComparisonFilterFunction() {
        return comparisonFilterFunction;
    }

    /**
     * Set the filter applied to every spectrum when performing the actual comparison. This
     * filter does not affect the spectra from which the consensus spectrum is build.
     * @param comparisonFilterFunction
     */
    public static void setComparisonFilterFunction(IFunction<List<IPeak>, List<IPeak>> comparisonFilterFunction) {
        CliSettings.comparisonFilterFunction = comparisonFilterFunction;
    }

    /**
     * The initial filter applied to every spectrum when loaded from file.
     * @return
     */
    public static IIntensityNormalizer getIntensityNormalizer() {
        return intensityNormalizer;
    }

    /**
     * Set the initial filter applied to every spectrum when loaded from file.
     * @param intensityNormalizer
     */
    public static void setIntensityNormalizer(IIntensityNormalizer intensityNormalizer) {
        CliSettings.intensityNormalizer = intensityNormalizer;
    }

    /**
     * This filter is applied to all spectra when loaded from the file - after the initial
     * spectrum filter.
     * @return
     */
    public static IFunction<List<IPeak>, List<IPeak>> getLoadingSpectrumFilter() {
        return loadingSpectrumFilter;
    }

    /**
     * Set the filter which is applied to all spectra when loaded from the file - after the initial
     * spectrum filter.
     * @param loadingSpectrumFilter
     */
    public static void setLoadingSpectrumFilter(IFunction<List<IPeak>, List<IPeak>> loadingSpectrumFilter) {
        CliSettings.loadingSpectrumFilter = loadingSpectrumFilter;
    }
}
