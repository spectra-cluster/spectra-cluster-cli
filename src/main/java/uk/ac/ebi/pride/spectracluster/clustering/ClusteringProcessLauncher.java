package uk.ac.ebi.pride.spectracluster.clustering;

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
    private List<Future<File>> fileFutures = new ArrayList<Future<File>>();


    public ClusteringProcessLauncher(ExecutorService executorService, File outputDirectory) {
        this.executorService = executorService;
        this.outputDirectory = outputDirectory;
    }

    @Override
    public void onNewResultFile(File binaryClusteringResultFile) {
        File outputFile = new File(outputDirectory, binaryClusteringResultFile.getName());
        BinaryFileClusteringCallable clusteringCallable = new BinaryFileClusteringCallable(outputFile, binaryClusteringResultFile);
        Future<File> fileFuture = executorService.submit(clusteringCallable);

        fileFutures.add(fileFuture);
    }

    public List<Future<File>> getFileFutures() {
        return Collections.unmodifiableList(fileFutures);
    }
}
