package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.indexing.ClusteringFileIndex;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.io.ClusteringFileReader;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterAppender;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
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
    private final List<String> peakListFilenames;
    private final List<String> clusteringFilenames;
    private final List<List<IndexElement>> fileIndices;
    private final List<ClusteringFileIndex> clusteringFileIndices;
    /**
     * If this is set, the comparison peak filtering is applied
     * to every input spectrum during loading.
     */
    private final boolean fastMode;
    private List<IBinaryClusteringResultListener> listeners = new ArrayList<IBinaryClusteringResultListener>();


    public BinarySpectrumReferenceWriter(List<String> peakListFilenames, List<String> clusteringFilenames, List<List<IndexElement>> fileIndices, List<ClusteringFileIndex> clusteringFileIndices, boolean fastMode) {
        this.peakListFilenames = peakListFilenames;
        this.clusteringFilenames = clusteringFilenames;
        this.fileIndices = fileIndices;
        this.clusteringFileIndices = clusteringFileIndices;
        this.fastMode = fastMode;
    }

    @Override
    public void writeSpectra(List<SpectrumReference> spectrumReferences, File outputFile, List<String> peakListFilenames) throws Exception {
        // sort the references
        Collections.sort(spectrumReferences);

        // open the file and write to it
        FileOutputStream outputStream = new FileOutputStream(outputFile);
        ObjectOutputStream objectOutputStream = new ObjectOutputStream(outputStream);

        double minMz = Double.MAX_VALUE, maxMz = 0;

        for (SpectrumReference spectrumReference : spectrumReferences) {
            try {
                checkInterrupt(outputStream);
                ICluster cluster = null;

                if (spectrumReference.getSpectrumIndex() == SpectrumReference.IS_CLUSTER) {
                    cluster = processCluster(spectrumReference);
                } else {
                    cluster = processSpectrum(spectrumReference);
                }

                // ignore empty spectra
                if (cluster == null) {
                    continue;
                }

                // update the statistics
                if (cluster.getPrecursorMz() < minMz) {
                    minMz = cluster.getPrecursorMz();
                }
                if (cluster.getPrecursorMz() > maxMz) {
                    maxMz = cluster.getPrecursorMz();
                }

                checkInterrupt(outputStream);

                // write it to the file
                BinaryClusterAppender.INSTANCE.appendCluster(objectOutputStream, cluster);
            }
            catch (Exception e) {
                String filename;
                if (spectrumReference.getSpectrumIndex() == SpectrumReference.IS_CLUSTER) {
                    filename = clusteringFilenames.get(spectrumReference.getFileId());
                } else {
                    filename = peakListFilenames.get(spectrumReference.getFileId());
                }

                throw new Exception("Error while processing spectrum reference (" + filename +
                        ", spec id = " + spectrumReference.getSpectrumId() + ", m/z = " +
                        spectrumReference.getPrecursorMz() + ")", e);
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

    private ICluster processCluster(SpectrumReference spectrumReference) throws Exception {
        int fileIndex = spectrumReference.getFileId();

        if (fileIndex >= clusteringFilenames.size())
            throw new Exception("Invalid file id for spectrum reference");

        String peakListFilename = clusteringFilenames.get(fileIndex);

        ClusteringFileReader reader = new ClusteringFileReader(new File(peakListFilename), clusteringFileIndices.get(fileIndex));

        // load the cluster
        uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster readerCluster =
                reader.readCluster(spectrumReference.getSpectrumId());

        ICluster cluster = SpectrumConverter.convertClusteringFileReaderCluster(readerCluster);

        // TODO: in fast mode, the consensus spectrum should be filtered

        return cluster;
    }

    private ICluster processSpectrum(SpectrumReference spectrumReference) throws Exception {
        int fileIndex = spectrumReference.getFileId();

        if (fileIndex >= peakListFilenames.size())
            throw new Exception("Invalid file id for spectrum reference");

        String peakListFilename = peakListFilenames.get(fileIndex);

        if (!readerPerFileIndex.containsKey(fileIndex))
            readerPerFileIndex.put(fileIndex, openFile(peakListFilenames.get(fileIndex), fileIndices.get(fileIndex)));

        JMzReader fileReader = readerPerFileIndex.get(fileIndex);

        // load the spectrum
        Spectrum spectrum = fileReader.getSpectrumByIndex(spectrumReference.getSpectrumIndex());

        // ignore empty spectra
        if (spectrum.getPeakList() == null || spectrum.getPeakList().size() < 1) {
            return null;
        }

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

        // convert to a cluster
        ICluster spectrumAsCluster = ClusterUtilities.asCluster(processedSpectrum);

        return spectrumAsCluster;
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
