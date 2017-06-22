package uk.ac.ebi.pride.spectracluster.implementation;

import uk.ac.ebi.pride.spectracluster.normalizer.IIntensityNormalizer;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.function.Functions;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.FractionTICPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.HighestNPeakFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveImpossiblyHighPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveIonContaminantsPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemovePrecursorPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveWindowPeaksFunction;

import java.util.List;
import java.util.function.Function;

/**
 * This is a (static) class that holds certain user
 * defined settings.
 * Created by jg on 08.04.16.
 */
public class ClusteringSettings {
    /**
     * Available filters to apply to the spectrum on loading
     * (added to the initial spectrum filter)
     */
    public static enum SPECTRUM_FILTER {
        IMMONIUM_IONS("immonium_ions", new RemoveIonContaminantsPeaksFunction(0.1F)),
        MZ_150("mz_150", new RemoveWindowPeaksFunction(150.0F, Float.MAX_VALUE)),
        MZ_200("mz_200", new RemoveWindowPeaksFunction(200.0F, Float.MAX_VALUE));

        public final String name;
        public final IFunction<ISpectrum, ISpectrum> filter;

        private SPECTRUM_FILTER(String name, IFunction<ISpectrum, ISpectrum> filter) {
            this.name = name;
            this.filter = filter;
        }

        public static SPECTRUM_FILTER getFilterForName(String name) {
            for (SPECTRUM_FILTER spectrumFilter : values()) {
                if (spectrumFilter.name.equals(name)) {
                    return spectrumFilter;
                }
            }

            return null;
        }
    }

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
        ClusteringSettings.initialSpectrumFilter = initialSpectrumFilter;
    }

    /**
     * Adds an additional filter to the current initial spectrum
     * filter.
     * @param filter The filter to be added
     */
    public static void addIntitalSpectrumFilter(IFunction<ISpectrum, ISpectrum> filter) {
        ClusteringSettings.initialSpectrumFilter = Functions.join(ClusteringSettings.initialSpectrumFilter, filter);
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
        ClusteringSettings.comparisonFilterFunction = comparisonFilterFunction;
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
        ClusteringSettings.intensityNormalizer = intensityNormalizer;
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
        ClusteringSettings.loadingSpectrumFilter = loadingSpectrumFilter;
    }

    /**
     * If this parameter is set to true, support for comment strings
     * in MGF files is disabled. This increases performance but only
     * works for MGF files that do not contain any comment strings.
     *
     * Comment support has been disabled since its strict interpretation
     * also treats ';' as comment characters which breaks some of the
     * PRIDE special fields.
     */
    public static boolean disableMGFCommentSupport = true;
}
