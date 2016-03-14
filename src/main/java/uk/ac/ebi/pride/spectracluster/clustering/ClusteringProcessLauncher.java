package uk.ac.ebi.pride.spectracluster.clustering;

import uk.ac.ebi.pride.spectracluster.util.ClusteringJobReference;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/**
 * Created by jg on 18.05.15.
 */
public class ClusteringProcessLauncher implements  IBinaryClusteringResultListener {
    private final ExecutorService executorService;
    private final File outputDirectory;
    private List<Future<ClusteringJobReference>> fileFutures = new ArrayList<Future<ClusteringJobReference>>();
    private final List<Float> thresholds;
    private final boolean fastMode;


    public ClusteringProcessLauncher(ExecutorService executorService, File outputDirectory, List<Float> thresholds, boolean fastMode) {
        this.executorService = executorService;
        this.outputDirectory = outputDirectory;
        this.thresholds = thresholds;
        this.fastMode = fastMode;
    }

    @Override
    public void onNewResultFile(BinaryClusterFileReference binaryBinaryClusterFileReferenceFile) {
        File outputFile = new File(outputDirectory, binaryBinaryClusterFileReferenceFile.getResultFile().getName());
        BinaryFileClusteringCallable clusteringCallable = new BinaryFileClusteringCallable(outputFile, binaryBinaryClusterFileReferenceFile.getResultFile(), thresholds, fastMode);
        Future<ClusteringJobReference> fileFuture = executorService.submit(clusteringCallable);

        fileFutures.add(fileFuture);
    }

    public List<Future<ClusteringJobReference>> getResultFileFutures() {
        return Collections.unmodifiableList(fileFutures);
    }
}
