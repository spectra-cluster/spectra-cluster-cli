package uk.ac.ebi.pride.spectracluster.merging;

import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusteringCallable;
import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;
import uk.ac.ebi.pride.spectracluster.util.ClusteringJobReference;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * This class is used to cluster the re-binned
 * files. These files should contain the overlapping
 * region of the previous file. Therefore, the first
 * file will not be clustered and in subsequent files
 * only the spectra within the set overlapping window
 * (ie in a window of 2 m/z the spectra within the file's
 * min m/z - [min m/z + 2 m/z])
 *
 * Created by jg on 30.05.15.
 */
public class BinaryFileMergingClusterer {
    private final int nJobs;
    private final File outputDirectory;
    private final List<Float> thresholds;
    private final boolean fastMode;
    private final double windowSize;
    private final boolean deleteInputFiles;
    private final File temporaryDirectory;

    private ExecutorService clusteringExecuteService;
    private List<IBinaryClusteringResultListener> listeners = new ArrayList<IBinaryClusteringResultListener>();
    private List<Future<ClusteringJobReference>> clusteringFutures = new ArrayList<Future<ClusteringJobReference>>();

    private List<BinaryClusterFileReference> resultFiles;

    public BinaryFileMergingClusterer(int nJobs, File outputDirectory, List<Float> thresholds, boolean fastMode,
                                      double windowSize, File temporaryDirectory) {
        this(nJobs, outputDirectory, thresholds, fastMode, windowSize, false, temporaryDirectory);
    }

    public BinaryFileMergingClusterer(int nJobs, File outputDirectory, List<Float> thresholds, boolean fastMode,
                                      double windowSize, boolean deleteInputFiles, File temporaryDirectory) {
        this.nJobs = nJobs;
        this.outputDirectory = outputDirectory;
        this.thresholds = thresholds;
        this.fastMode = fastMode;
        this.windowSize = windowSize;
        this.deleteInputFiles = deleteInputFiles;
        this.temporaryDirectory = temporaryDirectory;
    }

    public void clusterFiles(List<BinaryClusterFileReference> binaryFiles) throws Exception {
        launchClusteringJobs(binaryFiles);

        waitForCompletedJobs();
    }

    private void waitForCompletedJobs() throws Exception {
        boolean allDone = false;
        resultFiles = new ArrayList<BinaryClusterFileReference>();
        Set<Integer> completedJobs = new HashSet<Integer>();

        while (!allDone) {
            allDone = true;

            for (int i = 0; i < clusteringFutures.size(); i++) {
                if (Thread.currentThread().isInterrupted()) {
                    clusteringExecuteService.shutdownNow();
                    throw new InterruptedException();
                }

                if (completedJobs.contains(i))
                    continue;

                Future<ClusteringJobReference> fileFuture = clusteringFutures.get(i);

                if (!fileFuture.isDone()) {
                    allDone = false;
                }
                else {
                    // save the written file
                    BinaryClusterFileReference resultFile = fileFuture.get().getOutputFile();
                    resultFiles.add(resultFile);
                    // notify all listeners
                    notifyListeners(resultFile);

                    // delete the input file if set
                    if (deleteInputFiles) {
                        // ignore if the deletion fails
                        fileFuture.get().getInputFile().delete();
                    }

                    completedJobs.add(i);
                }
            }

            Thread.sleep(1000); // only check every second
        }

        // terminate the executor service
        clusteringExecuteService.awaitTermination(2, TimeUnit.SECONDS);
    }

    private void notifyListeners(BinaryClusterFileReference writtenFile) {
        for (IBinaryClusteringResultListener listener : listeners) {
            listener.onNewResultFile(writtenFile);
        }
    }

    private void launchClusteringJobs(List<BinaryClusterFileReference> binaryFiles) {
        clusteringExecuteService = Executors.newFixedThreadPool(nJobs);

        // make sure the files are sorted according to m/z
        Collections.sort(binaryFiles);

        for (int i = 0; i < binaryFiles.size(); i++) {
            if (Thread.currentThread().isInterrupted()) {
                clusteringExecuteService.shutdownNow();
                return;
            }

            BinaryClusterFileReference binaryClusterFileReference = binaryFiles.get(i);

            // ignore the first file since it will not be clustered
            if (i == 0) {
                notifyListeners(binaryClusterFileReference);
                continue;
            }

            // launch the clustering process
            double maxMz = binaryClusterFileReference.getMinMz() + windowSize;
            File outputFile = new File(outputDirectory, binaryClusterFileReference.getResultFile().getName());
            BinaryFileClusteringCallable clusteringCallable =
                    new BinaryFileClusteringCallable(outputFile, binaryClusterFileReference.getResultFile(),
                            thresholds, fastMode, 0, (float) maxMz, temporaryDirectory);
            Future<ClusteringJobReference> resultFileFuture = clusteringExecuteService.submit(clusteringCallable);
            clusteringFutures.add(resultFileFuture);
        }

        clusteringExecuteService.shutdown();
    }

    public void addListener(IBinaryClusteringResultListener listener) {
        listeners.add(listener);
    }

    public List<BinaryClusterFileReference> getResultFiles() {
        return Collections.unmodifiableList(resultFiles);
    }
}
