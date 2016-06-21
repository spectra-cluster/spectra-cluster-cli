package uk.ac.ebi.pride.spectracluster.util;

/**
 * Used to inform listeners on the progress of
 * a certain task (ie. number of files processed, etc.)
 *
 * Created by jg on 10.06.16.
 */
public class ProgressUpdate {
    public enum CLUSTERING_STAGE {
        CDF_LEARNING,
        CONVERSION,
        CLUSTERING,
        MERGING,
        OUTPUT
    };

    private String message;
    private Integer tasksCompleted;
    private Integer totalTasks;
    private Float percentCompleted;
    private CLUSTERING_STAGE clusteringStage;

    public ProgressUpdate(String message, CLUSTERING_STAGE clusteringStage) {
        this.message = message;
        this.clusteringStage = clusteringStage;
    }

    public ProgressUpdate(String message, CLUSTERING_STAGE clusteringStage, float percentCompleted) {
        this.message = message;
        this.clusteringStage = clusteringStage;
        this.percentCompleted = percentCompleted;
    }

    public ProgressUpdate(float percentCompleted) {
        this.percentCompleted = percentCompleted;
    }

    public ProgressUpdate(String message, CLUSTERING_STAGE clusteringStage, int tasksCompleted, int totalTasks) {
        this.message = message;
        this.clusteringStage = clusteringStage;
        this.tasksCompleted = tasksCompleted;
        this.totalTasks = totalTasks;
        this.percentCompleted = (float) tasksCompleted / totalTasks;
    }

    public boolean hasMessage() {
        return this.message != null;
    }

    public boolean hasTasksCompleted() {
        return this.tasksCompleted != null;
    }

    public boolean hasPercentCompleted() {
        return this.percentCompleted != null;
    }

    public String getMessage() {
        return message;
    }

    public Integer getTasksCompleted() {
        return tasksCompleted;
    }

    public Integer getTotalTasks() {
        return totalTasks;
    }

    public Float getPercentCompleted() {
        return percentCompleted;
    }

    public CLUSTERING_STAGE getClusteringStage() {
        return clusteringStage;
    }
}
