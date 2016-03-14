package uk.ac.ebi.pride.spectracluster.clustering;

import uk.ac.ebi.pride.spectracluster.clustering.ClusteringProcessLauncher;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.IBinaryClusteringResultListener;

import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Created by jg on 30.05.15.
 */
public class BinaryFileClusterer {
    private final int nJobs;
    private final File outputDirectory;
    private final List<Float> thresholds;
    private final boolean fastMode;

    private ExecutorService clusteringExecuteService;
    private ClusteringProcessLauncher clusteringProcessLauncher;
    private List<IBinaryClusteringResultListener> listeners = new ArrayList<IBinaryClusteringResultListener>();

    private List<BinaryClusterFileReference> resultFiles;

    public BinaryFileClusterer(int nJobs, File outputDirectory, List<Float> thresholds, boolean fastMode) {
        this.nJobs = nJobs;
        this.outputDirectory = outputDirectory;
        this.thresholds = thresholds;
        this.fastMode = fastMode;
    }

    public void clusterFiles(List<BinaryClusterFileReference> binaryFiles) throws Exception {
        launchClusteringJobs(binaryFiles);

        waitForCompletedJobs();
    }

    private void waitForCompletedJobs() throws Exception {
        List<Future<BinaryClusterFileReference>> fileFutures = clusteringProcessLauncher.getResultFileFutures();

        boolean allDone = false;
        resultFiles = new ArrayList<BinaryClusterFileReference>();
        Set<Integer> completedJobs = new HashSet<Integer>();

        while (!allDone) {
            allDone = true;

            for (int i = 0; i < fileFutures.size(); i++) {
                if (completedJobs.contains(i))
                    continue;

                Future<BinaryClusterFileReference> fileFuture = fileFutures.get(i);

                if (!fileFuture.isDone()) {
                    allDone = false;
                }
                else {
                    // save the written file
                    BinaryClusterFileReference resultFile = fileFuture.get();
                    resultFiles.add(resultFile);
                    // notify all listeners
                    notifyListeners(resultFile);

                    completedJobs.add(i);
                }

                Thread.sleep(1000); // only check every second
            }
        }

        // terminate the executor service
        clusteringExecuteService.awaitTermination(1, TimeUnit.MINUTES);
    }

    private void notifyListeners(BinaryClusterFileReference writtenFile) {
        for (IBinaryClusteringResultListener listener : listeners) {
            listener.onNewResultFile(writtenFile);
        }
    }

    private void launchClusteringJobs(List<BinaryClusterFileReference> binaryFiles) {
        clusteringExecuteService = Executors.newFixedThreadPool(nJobs);
        clusteringProcessLauncher = new ClusteringProcessLauncher(clusteringExecuteService, outputDirectory, thresholds, fastMode);

        for (BinaryClusterFileReference binaryFile : binaryFiles) {
            clusteringProcessLauncher.onNewResultFile(binaryFile);
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
