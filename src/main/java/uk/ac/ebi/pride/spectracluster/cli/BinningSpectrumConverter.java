package uk.ac.ebi.pride.spectracluster.cli;

import uk.ac.ebi.pride.spectracluster.binning.*;
import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;
import uk.ac.ebi.pride.spectracluster.spectra_list.IPeaklistScanner;
import uk.ac.ebi.pride.spectracluster.spectra_list.ParsingMgfScanner;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
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

    private List<IBinaryClusteringResultListener> listeners = new ArrayList<IBinaryClusteringResultListener>();
    private ExecutorService writingJobsExecutorService;
    private List<Future<File>> writtenBinaryFileFutures;

    private final File outputDirectory;
    private final int nJobs;
    private final boolean fastMode;
    private List<File> writtenFiles;
    private List<SpectrumReference> spectrumReferences;

    public BinningSpectrumConverter(File outputDirectory, int nJobs, boolean fastMode) {
        this.outputDirectory = outputDirectory;
        this.nJobs = nJobs;
        this.fastMode = fastMode;
    }

    public void processPeaklistFiles(String[] filenames) throws Exception {
        // pre-scan the spectrum references
        spectrumReferences = peaklistScanner.getSpectrumReferences(filenames);

        // bin the references
        List<List<SpectrumReference>> binnedSpectrumReferences = spectrumReferenceBinner.binSpectrumReferences(spectrumReferences);

        // create the files
        launchFileWritingJobs(binnedSpectrumReferences, filenames);

        // wait for the completed jobs
        waitForCompletedJobs();
    }

    private void waitForCompletedJobs() throws Exception {
        boolean allDone = false;
        writtenFiles = new ArrayList<File>();
        Set<Integer> completedWritingJobs = new HashSet<Integer>();

        while (!allDone) {
            allDone = true;

            for (int i = 0; i < writtenBinaryFileFutures.size(); i++) {
                if (completedWritingJobs.contains(i))
                    continue;

                Future<File> fileFuture = writtenBinaryFileFutures.get(i);

                if (!fileFuture.isDone()) {
                    allDone = false;
                }
                else {
                    // save the written file
                    File writtenFile = fileFuture.get();
                    writtenFiles.add(writtenFile);
                    // notify all listeners
                    notifyListeners(writtenFile);

                    completedWritingJobs.add(i);
                }
            }
        }

        // terminate the executor service
        writingJobsExecutorService.awaitTermination(1, TimeUnit.MINUTES);
    }

    private void notifyListeners(File writtenFile) {
        for (IBinaryClusteringResultListener listener : listeners) {
            listener.onNewResultFile(writtenFile);
        }
    }

    private void launchFileWritingJobs(List<List<SpectrumReference>> binnedSpectrumReferences, String[] filenames) throws Exception {
        if (outputDirectory == null || !outputDirectory.exists() || !outputDirectory.isDirectory())
            throw new Exception("Invalid output directory for converted spectra set");

        writingJobsExecutorService = Executors.newFixedThreadPool(nJobs);
        writtenBinaryFileFutures = new ArrayList<Future<File>>(binnedSpectrumReferences.size());

        // launch the jobs
        for (int i = 0; i < binnedSpectrumReferences.size(); i++) {
            List<SpectrumReference> spectrumReferences = binnedSpectrumReferences.get(i);
            File outputFile = generateOutputfile(i);

            BinarySpectrumReferenceWriterCallable writerCallable =
                    new BinarySpectrumReferenceWriterCallable(
                            filenames,
                            peaklistScanner.getFileIndices(),
                            spectrumReferences,
                            outputFile,
                            fastMode);

            Future<File> fileFuture = writingJobsExecutorService.submit(writerCallable);
            writtenBinaryFileFutures.add(fileFuture);
        }

        writingJobsExecutorService.shutdown();
    }

    public File generateOutputfile(int bin) {
        return new File(outputDirectory, "convertedSpectra_" + bin + ".cls");
    }

    public void addWrittenFileListener(IBinaryClusteringResultListener listener) {
        listeners.add(listener);
    }

    public List<File> getWrittenFiles() {
        return Collections.unmodifiableList(writtenFiles);
    }

    public List<List<IndexElement>> getFileIndices() {
        return Collections.unmodifiableList(peaklistScanner.getFileIndices());
    }

    public List<SpectrumReference> getSpectrumReferences() {
        return Collections.unmodifiableList(spectrumReferences);
    }
}
