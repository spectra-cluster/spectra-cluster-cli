package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.util.ClusterUtilities;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.*;

/**
 * Writes spectra stored as SpectrumReferences to a file and performs the complete
 * spectrum pre-processing (normalization, peak filtering etc.).
 *
 * Created by jg on 13.05.15.
 */
public class BinarySpectrumReferenceWriter implements ISpectrumReferenceWriter {
    private Map<Integer, JMzReader> readerPerFileIndex = new HashMap<Integer, JMzReader>();
    private final String[] peakListFilenames;
    private final List<List<IndexElement>> fileIndices;
    /**
     * If this is set, the comparison peak filtering is applied
     * to every input spectrum during loading.
     */
    private final boolean fastMode;
    private List<IBinaryClusteringResultListener> listeners = new ArrayList<IBinaryClusteringResultListener>();


    public BinarySpectrumReferenceWriter(String[] peakListFilenames, List<List<IndexElement>> fileIndices, boolean fastMode) {
        this.peakListFilenames = peakListFilenames;
        this.fileIndices = fileIndices;
        this.fastMode = fastMode;
    }

    @Override
    public void writeSpectra(List<SpectrumReference> spectrumReferences, File outputFile, String[] peakListFilenames) throws Exception {
        // sort the references
        Collections.sort(spectrumReferences);

        // open the file and write to it
        FileOutputStream outputStream = new FileOutputStream(outputFile);
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(outputStream);

        double minMz = Double.MAX_VALUE, maxMz = 0;

        for (SpectrumReference spectrumReference : spectrumReferences) {
            try {
                checkInterrupt(outputStream);

                int fileIndex = spectrumReference.getFileId();

                if (fileIndex >= peakListFilenames.length)
                    throw new Exception("Invalid file id for spectrum reference");

                String peakListFilename = peakListFilenames[fileIndex];

                if (!readerPerFileIndex.containsKey(fileIndex))
                    readerPerFileIndex.put(fileIndex, openFile(peakListFilenames[fileIndex], fileIndices.get(fileIndex)));

                JMzReader fileReader = readerPerFileIndex.get(fileIndex);

                // load the spectrum
                Spectrum spectrum = fileReader.getSpectrumByIndex(spectrumReference.getSpectrumIndex());

                // pre-process the spectrum
                ISpectrum convertedSpectrum = SpectrumConverter.convertJmzReaderSpectrum(spectrum, spectrumReference.getSpectrumId(), peakListFilename);
                ISpectrum processedSpectrum = ClusteringSettings.getInitialSpectrumFilter().apply(convertedSpectrum);
                // normalize the spectrum
                processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(
                        processedSpectrum, ClusteringSettings.getIntensityNormalizer().normalizePeaks(processedSpectrum.getPeaks()));

                // apply the comparison filter function right when loading the specturm in fast mode.
                if (fastMode) {
                    processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(processedSpectrum, ClusteringSettings.getComparisonFilterFunction().apply(processedSpectrum.getPeaks()));
                }
                else {
                    processedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(processedSpectrum, ClusteringSettings.getLoadingSpectrumFilter().apply(processedSpectrum.getPeaks()));
                }

                // update the statistics
                if (processedSpectrum.getPrecursorMz() < minMz) {
                    minMz = processedSpectrum.getPrecursorMz();
                }
                if (processedSpectrum.getPrecursorMz() > maxMz) {
                    maxMz = processedSpectrum.getPrecursorMz();
                }

                checkInterrupt(outputStream);

                // write it to the file
                ICluster spectrumAsCluster = ClusterUtilities.asCluster(processedSpectrum);
                BinaryClusterAppender.INSTANCE.appendCluster(objectOutputStream, spectrumAsCluster);
            }
            catch (Exception e) {
                throw new Exception("Error while processing spectrum reference (file id = " + spectrumReference.getFileId() + " > " + peakListFilenames[spectrumReference.getFileId()] + ", spec index = " + spectrumReference.getSpectrumIndex() + ", m/z = " + spectrumReference.getPrecursorMz() + ")", e);
            }
        }

        // close the file
        BinaryClusterAppender.INSTANCE.appendEnd(objectOutputStream);
        objectOutputStream.close();
        outputStream.close();

        // notify the listeners
        for (IBinaryClusteringResultListener listener : listeners)
            listener.onNewResultFile(new BinaryClusterFileReference(outputFile, minMz, maxMz, spectrumReferences.size()));
    }

    private void checkInterrupt(OutputStream outputStream) throws Exception {
        if (Thread.currentThread().isInterrupted()) {
            outputStream.close();
            throw new InterruptedException();
        }
    }

    private JMzReader openFile(String peakListFilename, List<IndexElement> fileIndex) throws Exception {
        if (peakListFilename.toLowerCase().endsWith(".mgf")) {
            MgfFile mgfFile = new MgfFile(new File(peakListFilename), fileIndex, true);
            mgfFile.setDisableCommentSupport(ClusteringSettings.disableMGFCommentSupport);

            return mgfFile;
        }

        throw new Exception("Unknown file extension encountered: " + peakListFilename);
    }

    public void addListener(IBinaryClusteringResultListener listener) {
        listeners.add(listener);
    }
}
