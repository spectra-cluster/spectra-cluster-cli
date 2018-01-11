package uk.ac.ebi.pride.spectracluster.binning;

import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.implementation.SpectraClusterStandalone;
import uk.ac.ebi.pride.spectracluster.spectra_list.IPeaklistScanner;
import uk.ac.ebi.pride.spectracluster.spectra_list.ParsingClusteringScanner;
import uk.ac.ebi.pride.spectracluster.spectra_list.ParsingMgfScanner;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.util.IProgressListener;
import uk.ac.ebi.pride.spectracluster.util.ProgressUpdate;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Processes the input files and bins the spectra
 * based on the used binner. Each bin is then
 * written into one binary file.
 *
 * Created by jg on 30.05.15.
 */
public class BinningSpectrumConverter {
    private ISpectrumReferenceBinner spectrumReferenceBinner = new ReferenceMzBinner();
    private IPeaklistScanner peaklistScanner = new ParsingMgfScanner();
    private ParsingClusteringScanner clusteringScanner = new ParsingClusteringScanner();

    private List<IBinaryClusteringResultListener> listeners = new ArrayList<IBinaryClusteringResultListener>();
    private List<IProgressListener> progressListeners = new ArrayList<IProgressListener>();
    private ExecutorService writingJobsExecutorService;
    private List<Future<BinaryClusterFileReference>> writtenBinaryFileFutures;

    private final File outputDirectory;
    private final int nJobs;
    private final boolean fastMode;
    private List<BinaryClusterFileReference> writtenFiles;
    private List<SpectrumReference> spectrumReferences;

    private ClusteringSettings.LOADING_MODE loadingMode = ClusteringSettings.DEFAULT_LOADING_MODE;

    public BinningSpectrumConverter(File outputDirectory, int nJobs, boolean fastMode) {
        this.outputDirectory = outputDirectory;
        this.nJobs = nJobs;
        this.fastMode = fastMode;
    }

    public BinningSpectrumConverter(File outputDirectory, int nJobs, boolean fastMode, ISpectrumReferenceBinner spectrumReferenceBinner) {
        this(outputDirectory, nJobs, fastMode);
        this.spectrumReferenceBinner = spectrumReferenceBinner;
    }

    public void processPeaklistFiles(String[] filenames) throws Exception {
        // process MGF and .clustering files separately
        List<String> mgfFilenames = new ArrayList<>(filenames.length);
        List<String> clusteringFilenames = new ArrayList<>(filenames.length);

        ((ParsingMgfScanner) peaklistScanner).setLoadingMode(loadingMode);

        for (String filename : filenames) {
            if (filename.toLowerCase().endsWith(".mgf")) {
                mgfFilenames.add(filename);
            }
            else if (filename.toLowerCase().endsWith(".clustering")) {
                clusteringFilenames.add(filename);
            }
            else {
                throw new Exception("Unsupported filetype in " + filename);
            }
        }

        // pre-scan the spectrum references
        spectrumReferences = new ArrayList<>();
        if (mgfFilenames.size() > 0) {
            spectrumReferences.addAll(
                    peaklistScanner.getSpectrumReferences(mgfFilenames.toArray(new String[mgfFilenames.size()])));
        }
        if (clusteringFilenames.size() > 0) {
            spectrumReferences.addAll(
                    clusteringScanner.getSpectrumReferences(clusteringFilenames.toArray(new String[clusteringFilenames.size()])));
        }

        // bin the references
        List<List<SpectrumReference>> binnedSpectrumReferences = spectrumReferenceBinner.binSpectrumReferences(spectrumReferences);

        // create the files
        launchFileWritingJobs(binnedSpectrumReferences, mgfFilenames, clusteringFilenames);

        // wait for the completed jobs
        waitForCompletedJobs();
    }

    private void waitForCompletedJobs() throws Exception {
        try {
            boolean allDone = false;
            writtenFiles = new ArrayList<BinaryClusterFileReference>();
            Set<Integer> completedWritingJobs = new HashSet<Integer>();

            while (!allDone) {
                allDone = true;

                for (int i = 0; i < writtenBinaryFileFutures.size(); i++) {
                    if (Thread.currentThread().isInterrupted()) {
                        throw new InterruptedException();
                    }

                    if (completedWritingJobs.contains(i))
                        continue;

                    Future<BinaryClusterFileReference> fileFuture = writtenBinaryFileFutures.get(i);

                    if (!fileFuture.isDone()) {
                        allDone = false;
                    } else {
                        // save the written file
                        BinaryClusterFileReference writtenFile = fileFuture.get();
                        writtenFiles.add(writtenFile);
                        // notify all listeners
                        notifyListeners(writtenFile);
                        notifyProgressListeners(writtenFiles.size(), writtenBinaryFileFutures.size());

                        completedWritingJobs.add(i);
                    }
                }
            }

            // terminate the executor service
            writingJobsExecutorService.awaitTermination(1, TimeUnit.SECONDS);
        }
        catch (InterruptedException e) {
            writingJobsExecutorService.shutdownNow();
            throw(e);
        }
    }

    private void notifyListeners(BinaryClusterFileReference writtenFile) {
        for (IBinaryClusteringResultListener listener : listeners) {
            listener.onNewResultFile(writtenFile);
        }
    }

    private void notifyProgressListeners(int writtenFiles, int totalFiles) {
        ProgressUpdate progressUpdate = new ProgressUpdate(
                String.format("Converting %d / %d binary files.", writtenFiles, totalFiles),
                ProgressUpdate.CLUSTERING_STAGE.CONVERSION, writtenFiles, totalFiles);

        for (IProgressListener progressListener : progressListeners) {
            progressListener.onProgressUpdate(progressUpdate);
        }
    }

    private void launchFileWritingJobs(List<List<SpectrumReference>> binnedSpectrumReferences, List<String> peaklistFilenames,
                                       List<String> clusteringFilenames) throws Exception {
        if (outputDirectory == null || !outputDirectory.exists() || !outputDirectory.isDirectory())
            throw new Exception("Invalid output directory for converted spectra set");

        writingJobsExecutorService = Executors.newFixedThreadPool(nJobs);
        writtenBinaryFileFutures = new ArrayList<Future<BinaryClusterFileReference>>(binnedSpectrumReferences.size());

        // launch the jobs
        for (int i = 0; i < binnedSpectrumReferences.size(); i++) {
            List<SpectrumReference> spectrumReferences = binnedSpectrumReferences.get(i);

            if (spectrumReferences.size() < 1) {
                continue;
            }

            File outputFile = generateOutputfile(spectrumReferences);

            BinarySpectrumReferenceWriterCallable writerCallable =
                    new BinarySpectrumReferenceWriterCallable(
                            peaklistFilenames,
                            clusteringFilenames,
                            peaklistScanner.getFileIndices(),
                            clusteringScanner.getClusteringFileIndices(),
                            spectrumReferences,
                            outputFile,
                            fastMode);

            Future<BinaryClusterFileReference> fileFuture = writingJobsExecutorService.submit(writerCallable);
            writtenBinaryFileFutures.add(fileFuture);
        }

        writingJobsExecutorService.shutdown();
    }

    public File generateOutputfile(List<SpectrumReference> spectrumReferences) {
        float minMz = (float) spectrumReferences.stream()
                .mapToDouble(SpectrumReference::getPrecursorMz)
                .min()
                .orElse(0.0);
        float maxMz = (float) spectrumReferences.stream()
                .mapToDouble(SpectrumReference::getPrecursorMz)
                .max()
                .orElse(9999.0);

        String filename = String.format("convertedSpectra_%04d_%04d.cls",
                (int) Math.floor(minMz),
                (int) Math.ceil(maxMz));

        return new File(outputDirectory, filename);
    }

    public void addWrittenFileListener(IBinaryClusteringResultListener listener) {
        listeners.add(listener);
    }

    public List<BinaryClusterFileReference> getWrittenFiles() {
        return Collections.unmodifiableList(writtenFiles);
    }

    public List<List<IndexElement>> getFileIndices() {
        return Collections.unmodifiableList(peaklistScanner.getFileIndices());
    }

    public List<SpectrumReference> getSpectrumReferences() {
        return Collections.unmodifiableList(spectrumReferences);
    }

    public void addProgressListener(IProgressListener listener) {
        progressListeners.add(listener);
    }

    public ClusteringSettings.LOADING_MODE getLoadingMode() {
        return loadingMode;
    }

    public void setLoadingMode(ClusteringSettings.LOADING_MODE loadingMode) {
        this.loadingMode = loadingMode;
    }
}
