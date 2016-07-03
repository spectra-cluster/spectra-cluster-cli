package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReferenceFileComparator;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.File;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * Created by jg on 07.04.16.
 */
public class CdfComparisonCallable implements Callable<CdfResult> {
    private final List<CdfLearner.SpectrumMatch> spectrumMatches;
    private final List<SpectrumReference> spectrumReferences;
    private final List<SpectrumReference> allSpectrumReferences;
    private final String[] filenames;
    private final List<List<IndexElement>> fileIndices;
    private final ISimilarityChecker similarityChecker = new CombinedFisherIntensityTest(Defaults.getFragmentIonTolerance());

    /**
     * Score increments to use when creating the cumulative
     * distribution function
     */
    public final double SCORE_INCREMENTS = 0.5;

    /**
     * Creates a new CdfComparisonCallable which will compare the passed spectra
     * defined in the spectrumMatches list and return the result of theses
     * comparisons as a CumulativeDistributionFunction object.
     * @param spectrumMatches The spectrum matches to compare.
     * @param spectrumReferences Spectrum references to which the SpectrumMatch indices relate to
     * @param filenames List of peaklist filenames.
     * @param fileIndices File indices of these peaklist files used to speedup random access to the files.
     */
    public CdfComparisonCallable(List<CdfLearner.SpectrumMatch> spectrumMatches, List<SpectrumReference> spectrumReferences, String[] filenames, List<List<IndexElement>> fileIndices) {
        this.spectrumMatches = spectrumMatches;
        this.filenames = filenames;
        this.fileIndices = fileIndices;
        this.allSpectrumReferences = spectrumReferences;

        // only store the spectrum references which are actually needed
        Set<Integer> spectrumReferenceIndices = new HashSet<Integer>();
        for (CdfLearner.SpectrumMatch match : spectrumMatches) {
            spectrumReferenceIndices.add(match.getSpecIndex1());
            spectrumReferenceIndices.add(match.getSpecIndex2());
        }

        this.spectrumReferences = new ArrayList<SpectrumReference>(spectrumReferenceIndices.size());

        for (Integer spectrumReferenceIndex : spectrumReferenceIndices) {
            this.spectrumReferences.add(spectrumReferences.get(spectrumReferenceIndex));
        }

        // sort based on the passed input file
        Collections.sort(this.spectrumReferences, SpectrumReferenceFileComparator.INSTANCE);
    }

    @Override
    public CdfResult call() throws Exception {
        // load all spectra
        Map<SpectrumReference, ISpectrum> loadedSpectra = loadSpectra();

        // start the comparisons
        CdfResult cdfResult = new CdfResult(0.5);

        for (CdfLearner.SpectrumMatch spectrumMatch : spectrumMatches) {
            ISpectrum spectrum1 = loadedSpectra.get(allSpectrumReferences.get(spectrumMatch.getSpecIndex1()));
            ISpectrum spectrum2 = loadedSpectra.get(allSpectrumReferences.get(spectrumMatch.getSpecIndex2()));

            double similarity = similarityChecker.assessSimilarity(spectrum1, spectrum2);

            cdfResult.saveRandomMatchResult(similarity);
        }

        return cdfResult;
    }

    private Map<SpectrumReference, ISpectrum> loadSpectra() throws Exception {
        Map<SpectrumReference, ISpectrum> loadedSpectra = new HashMap<SpectrumReference, ISpectrum>();

        // spectrum references are sorted by input file. Therefore a file
        // only has to be opened once
        int currentFileIndex = -1;
        JMzReader currentFileReader = null;

        for (SpectrumReference spectrumReference : spectrumReferences) {
            // open the file if a new file id is encountered
            if (spectrumReference.getFileId() != currentFileIndex) {
                currentFileIndex = spectrumReference.getFileId();
                currentFileReader = openFile(filenames[currentFileIndex], fileIndices.get(currentFileIndex));
            }

            // load the spectrum
            uk.ac.ebi.pride.tools.jmzreader.model.Spectrum spectrum = currentFileReader.getSpectrumByIndex(spectrumReference.getSpectrumIndex());

            // pre-process the spectrum
            ISpectrum convertedSpectrum = SpectrumConverter.convertJmzReaderSpectrum(spectrum, spectrumReference.getSpectrumId(), "testfile");
            ISpectrum processedSpectrum = ClusteringSettings.getInitialSpectrumFilter().apply(convertedSpectrum);
            // normalize the spectrum
            processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(
                    processedSpectrum, ClusteringSettings.getIntensityNormalizer().normalizePeaks(processedSpectrum.getPeaks()));

            // only retain the N-highest peaks
            processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(processedSpectrum,
                    ClusteringSettings.getLoadingSpectrumFilter().apply(processedSpectrum.getPeaks()));

            // perform the actual "comparison" filter
            processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(processedSpectrum,
                    ClusteringSettings.getComparisonFilterFunction().apply(processedSpectrum.getPeaks()));

            loadedSpectra.put(spectrumReference, processedSpectrum);
        }

        return loadedSpectra;
    }

    private JMzReader openFile(String peakListFilename, List<IndexElement> fileIndex) throws Exception {
        if (peakListFilename.toLowerCase().endsWith(".mgf")) {
            MgfFile mgfFile = new MgfFile(new File(peakListFilename), fileIndex, true);
            mgfFile.setDisableCommentSupport(true);

            return mgfFile;
        }

        throw new Exception("Unknown file extension encountered: " + peakListFilename);
    }
}
