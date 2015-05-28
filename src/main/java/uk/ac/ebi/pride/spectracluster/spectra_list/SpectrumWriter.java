package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.spectracluster.util.SpectrumUtilities;
import uk.ac.ebi.pride.spectracluster.util.function.IFunction;
import uk.ac.ebi.pride.spectracluster.util.function.peak.FractionTICPeakFunction;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.*;
import java.util.*;

/**
 * Writes spectra stored as SpectrumReferences to a file.
 *
 * Created by jg on 13.05.15.
 */
public class SpectrumWriter {
    public static final IFunction<ISpectrum, ISpectrum> filterFunction = Defaults.getDefaultPeakFilter();
    public static final IFunction<List<IPeak>, List<IPeak>> comparisonFilterFunction = new FractionTICPeakFunction(0.5f, 20);
    private Map<Integer, JMzReader> readerPerFileIndex = new HashMap<Integer, JMzReader>();
    private final String[] peakListFilenames;
    private final List<List<IndexElement>> fileIndices;
    private List<IBinaryClusteringResultListener> listeners = new ArrayList<IBinaryClusteringResultListener>();


    public SpectrumWriter(String[] peakListFilenames, List<List<IndexElement>> fileIndices) {
        this.peakListFilenames = peakListFilenames;
        this.fileIndices = fileIndices;
    }

    public void writeSpectra(List<SpectrumReference> spectrumReferences, File outputFile) throws Exception {
        // sort the references
        Collections.sort(spectrumReferences);

        // open the file and write to it
        FileOutputStream outputStream = new FileOutputStream(outputFile);
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(outputStream);

        for (SpectrumReference spectrumReference : spectrumReferences) {
            int fileIndex = spectrumReference.getFileId();

            if (fileIndex >= peakListFilenames.length)
                throw new Exception("Invalid file id for spectrum reference");

            if (!readerPerFileIndex.containsKey(fileIndex))
                readerPerFileIndex.put(fileIndex, openFile(peakListFilenames[fileIndex], fileIndices.get(fileIndex)));

            JMzReader fileReader = readerPerFileIndex.get(fileIndex);

            // load the spectrum
            Spectrum spectrum = fileReader.getSpectrumByIndex(spectrumReference.getSpectrumIndex());

            // pre-process the spectrum
            ISpectrum processedSpectrum = filterFunction.apply(SpectrumConverter.convertJmzReaderSpectrum(spectrum, spectrumReference.getSpectrumId()));
            // normalize the spectrum
            processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(
                    processedSpectrum, Defaults.getDefaultIntensityNormalizer().normalizePeaks(processedSpectrum.getPeaks()));

            // do the radical peak filtering already now
            processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(processedSpectrum, comparisonFilterFunction.apply(processedSpectrum.getPeaks()));

            // write it to the file
            BinaryClusterAppender.INSTANCE.appendCluster(objectOutputStream, ClusterUtilities.asCluster(processedSpectrum));
        }

        // close the file
        BinaryClusterAppender.INSTANCE.appendEnd(objectOutputStream);
        objectOutputStream.close();
        outputStream.close();

        // notify the listeners
        for (IBinaryClusteringResultListener listener : listeners)
            listener.onNewResultFile(outputFile);
    }

    private JMzReader openFile(String peakListFilename, List<IndexElement> fileIndex) throws Exception {
        if (peakListFilename.toLowerCase().endsWith(".mgf"))
            return new MgfFile(new File(peakListFilename), fileIndex);

        throw new Exception("Unknown file extension encountered: " + peakListFilename);
    }

    public void addListener(IBinaryClusteringResultListener listener) {
        listeners.add(listener);
    }
}
