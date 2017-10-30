package uk.ac.ebi.pride.spectracluster.implementation;

import org.apache.commons.io.filefilter.FileFilterUtils;
import uk.ac.ebi.pride.spectracluster.binning.BinningSpectrumConverter;
import uk.ac.ebi.pride.spectracluster.cdf.SpectraPerBinNumberComparisonAssessor;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryFileClusterer;
import uk.ac.ebi.pride.spectracluster.conversion.MergingCGFConverter;
import uk.ac.ebi.pride.spectracluster.io.CGFSpectrumIterable;
import uk.ac.ebi.pride.spectracluster.io.DotClusterClusterAppender;
import uk.ac.ebi.pride.spectracluster.merging.BinaryFileMergingClusterer;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.util.BinaryFileScanner;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.IProgressListener;
import uk.ac.ebi.pride.spectracluster.util.ProgressUpdate;
import uk.ac.ebi.pride.spectracluster.util.function.Functions;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.HighestNSpectrumPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveReporterIonPeaksFunction;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * This class takes care of the actual clustering
 * logic. Settings used during the clustering process
 * are either taken from ClusteringSettings or from
 * the Defaults class which is part of the spectra-cluster
 * API.
 *
 * Settings taken from the spectra-cluster Defaults class:
 * - precursor tolerance
 * - fragment tolerance
 * - custom cumulative distribution function
 *
 * Settings taken from ClusteringSettings:
 * - peak filtering functions
 *
 * Created by jg on 23.06.16.
 */
public class SpectraClusterStandalone {
    private List<IProgressListener> progressListeners = new ArrayList<IProgressListener>();
    private boolean verbose = false;
    private boolean keepBinaryFiles = false;
    private File temporaryDirectory;
    private RemoveReporterIonPeaksFunction.REPORTER_TYPE reporterType = null;
    private int parallelJobs;
    private boolean useFastMode = false;

    /**
     * This variable only exists for debugging purposes. If set to false, temporary files
     * (ie. intermediate files created during merging) are kept.
     */
    private boolean deleteTemporaryFiles = true;

    /**
     * Create a new SpectraClusterStandalone object.
     * @throws Exception In case the temporary directory cannot be initialized (ie. created)
     */
    public SpectraClusterStandalone() throws Exception {
        // initialize the temporary directory
        this.temporaryDirectory = createTemporaryDirectory("spectra_cluster_cli");
        // default number of parallel jobs is the number of cores
        this.parallelJobs = Runtime.getRuntime().availableProcessors();
    }

    /**
     * Clusters the defined peak list files.
     *
     * @param peaklistFiles Peak list files to cluster. Currently, only MGF files are supported.
     * @param clusteringThresholds The clustering thresholds to use. The thresholds are processed according
     *                             to the order of the list. For most similarity metrics a higher threshold
     *                             means a higher accuracy. In these cases the list should be sorted from
     *                             highest to lowest (clustering should start with the highest accuracy and
     *                             decrease with each subsequent clustering round).
     * @param resultFile File object were the final result file (.clustering) format will be written to.
     */
    public void clusterPeaklistFiles(List<File> peaklistFiles, List<Float> clusteringThresholds, File resultFile)
        throws Exception {
        File binarySpectraDirectory = new File(temporaryDirectory, "spectra");

        if (reporterType != null) {
            ClusteringSettings.addIntitalSpectrumFilter(
                    new RemoveReporterIonPeaksFunction(Defaults.getFragmentIonTolerance(), reporterType)
            );
        }

        // convert the binary files
        List<BinaryClusterFileReference> binaryFiles = convertInputFiles(peaklistFiles, binarySpectraDirectory);

        clusterBinaryFiles(binaryFiles, clusteringThresholds, resultFile);
    }

    /**
     * Cluster the converted or existing binary files.
     * @param binaryFiles BinaryClusterFileRefernece representing the binary files to cluster.
     * @param clusteringThresholds The clustering thresholds to use. The thresholds are processed according
     *                             to the order of the list. For most similarity metrics a higher threshold
     *                             means a higher accuracy. In these cases the list should be sorted from
     *                             highest to lowest (clustering should start with the highest accuracy and
     *                             decrease with each subsequent clustering round).
     * @param resultFile
     */
    private void clusterBinaryFiles(List<BinaryClusterFileReference> binaryFiles, List<Float> clusteringThresholds,
                                    File resultFile) throws Exception {
        // the actual clustering
        List<BinaryClusterFileReference> clusteredFiles = clusterFiles(binaryFiles, clusteringThresholds);

        // merge the results
        File combinedResultFile = mergeClusteringResults(clusteredFiles, clusteringThresholds);

        // create the output file
        convertCgfToClustering(combinedResultFile, resultFile, clusteringThresholds.get(clusteringThresholds.size() - 1));

        if (!combinedResultFile.delete()) {
            // TODO: add notification at later stage
        }
    }

    /**
     * Clusters existing binary files in the defined directory. This directory must have been used as
     * "temporary" directory in a previous clustering run. Several parameters do not take effect if this
     * approach is used, most importantly the initial peak filter is not applied in this approach.
     * @param clusteringDirectory The directory previously set as "temporary" directory.
     * @param clusteringThresholds The clustering thresholds to use. The thresholds are processed according
     *                             to the order of the list. For most similarity metrics a higher threshold
     *                             means a higher accuracy. In these cases the list should be sorted from
     *                             highest to lowest (clustering should start with the highest accuracy and
     *                             decrease with each subsequent clustering round).
     * @param resultFile File object were the final result file (.clustering) format will be written to.
     */
    public void clusterExistingBinaryFiles(File clusteringDirectory, List<Float> clusteringThresholds, File resultFile)
            throws Exception {
        File binarySpectraDirectory = new File(clusteringDirectory, "spectra");

        // pre-scan the binary files
        List<BinaryClusterFileReference> binaryFiles = findExistingBinaryFiles(binarySpectraDirectory);

        if (binaryFiles.size() < 1) {
            throw new Exception("No existing binary files found.");
        }

        clusterBinaryFiles(binaryFiles, clusteringThresholds, resultFile);
    }

    /**
     * Tests whether the passed directory contains the expected structure found in a directory previously
     * used as "temporary" directory during a clustering run.
     * @param directory Directory to check.
     * @return
     */
    public boolean isValidClusteringDirectory(File directory) {
        if (!directory.exists())
            return false;

        // make sure the spectra directory exists
        File spectraDirectory = new File(directory, "spectra");

        return spectraDirectory.exists();
    }

    /**
     * Adds a new IProgressListener to receive status updates.
     * @param progressListener
     */
    public void addProgressListener(IProgressListener progressListener) {
        progressListeners.add(progressListener);
    }

    /**
     * Remove a progressListener to receive updates.
     * @param progressListener
     */
    public void removeProgressListener(IProgressListener progressListener) {
        progressListeners.remove(progressListener);
    }

    /**
     * If set more status messages are sent.
     * @return
     */
    public boolean isVerbose() {
        return verbose;
    }

    /**
     * If set to true more status messages are sent to the
     * added IProgressListenerS
     * @param verbose
     */
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    /**
     * Defines whether the generated binary files are deleted
     * after the respective clustering step is completed or not.
     * @return
     */
    public boolean isKeepBinaryFiles() {
        return keepBinaryFiles;
    }

    /**
     * Defines whether the generated binary files are deleted
     * after the respective clustering step is completed or not.
     * @param keepBinaryFiles
     */
    public void setKeepBinaryFiles(boolean keepBinaryFiles) {
        this.keepBinaryFiles = keepBinaryFiles;
    }

    /**
     * Temporary directory were all temporary clustering files will
     * be created. These include the converted binary files (input files
     * are converted before clustering) as well as intermediate files during
     * round-based clustering.
     * The default behaviour is to create a temporary directory in the OS'
     * standard location.
     * @return
     */
    public File getTemporaryDirectory() {
        return temporaryDirectory;
    }

    /**
     * Temporary directory were all temporary clustering files will
     * be created. These include the converted binary files (input files
     * are converted before clustering) as well as intermediate files during
     * round-based clustering.
     * The default behaviour is to create a temporary directory in the OS'
     * standard location.
     * @param temporaryDirectory
     */
    public void setTemporaryDirectory(File temporaryDirectory) {
        // try to delete the currently set temporary directory but ignore any problems
        if (this.temporaryDirectory != null && this.temporaryDirectory.exists()) {
            this.temporaryDirectory.delete();
        }

        this.temporaryDirectory = temporaryDirectory;
    }

    /**
     * If this parameter is set to not null the corresponding reporter peaks
     * are removed from the spectra. Set to null to disable reporter peak
     * filtering.
     * @return
     */
    public RemoveReporterIonPeaksFunction.REPORTER_TYPE getReporterType() {
        return reporterType;
    }

    /**
     * If this parameter is set to not null the corresponding reporter peaks
     * are removed from the spectra. Set to null to disable reporter peak
     * filtering.
     * @param reporterType
     */
    public void setReporterType(RemoveReporterIonPeaksFunction.REPORTER_TYPE reporterType) {
        this.reporterType = reporterType;
    }

    /**
     * The number of parallel jobs to launch during the clustering process.
     * @return
     */
    public int getParallelJobs() {
        return parallelJobs;
    }

    /**
     * The number of parallel jobs to launch during the clustering process.
     * @param parallelJobs
     */
    public void setParallelJobs(int parallelJobs) {
        this.parallelJobs = parallelJobs;
    }

    /**
     * If set radical peak filtering is already performed during spectrum loading
     * and not for the actual comparisons. This increases clustering performance but
     * may decrease accuracy.
     * @return
     */
    public boolean isUseFastMode() {
        return useFastMode;
    }

    /**
     * If set radical peak filtering is already performed during spectrum loading
     * and not for the actual comparisons. This increases clustering performance but
     * may decrease accuracy.
     * @param useFastMode
     */
    public void setUseFastMode(boolean useFastMode) {
        this.useFastMode = useFastMode;
    }

    /**
     * This variable only exists for debugging purposes. If set to false, temporary files
     * (ie. intermediate files created during merging) are kept.
     * @return
     */
    public boolean isDeleteTemporaryFiles() {
        return deleteTemporaryFiles;
    }

    /**
     * This variable only exists for debugging purposes. If set to false, temporary files
     * (ie. intermediate files created during merging) are kept.
     * @param deleteTemporaryFiles
     */
    public void setDeleteTemporaryFiles(boolean deleteTemporaryFiles) {
        this.deleteTemporaryFiles = deleteTemporaryFiles;
    }

    /**
     * Create a temporary directory in the OS' default location.
     * @param prefix The prefix to use for the temporary directory.
     * @return A File object pointing to the newly created temporary directory.
     * @throws Exception
     */
    private File createTemporaryDirectory(String prefix) throws Exception {
        return createTemporaryDirectory(prefix, null);
    }

    /**
     * Creates a temporary directory with the specified prefix.
     * @param prefix The prefix to use for the temporary directory.
     * @param tmpDirectory Parent directory to use instead of the OS' default temporary directory. Set
     *                     to null to ignore.
     * @return A File object pointing to the newly created temporary directory.
     * @throws Exception
     */
    private File createTemporaryDirectory(String prefix, File tmpDirectory) throws Exception {
        // first a temporary file is created. This is then deleted and the generated name re-used
        // to create the directory.
        File tmpFile;

        if (tmpDirectory != null && tmpDirectory.isDirectory()) {
            tmpFile = File.createTempFile(prefix, "", tmpDirectory);
        }
        else {
            tmpFile = File.createTempFile(prefix, "");
        }

        if (!tmpFile.delete())
            throw new Exception("Failed to delete temporary file");

        if (!tmpFile.mkdir())
            throw new Exception("Failed to create temporary directory");

        return tmpFile;
    }

    /**
     * This function either converts the input files and writes them to the
     * binary format or simply re-loads already existing binary files.
     * @param peaklistFilenames Array of input peaklist filenames.
     * @param binarySpectraDirectory Directory into which the binary files should be written to (or from where they are loaded).
     * @return A list of BinaryClusterFileReferenceS representing the generated / existing binary files.
     */
    private List<BinaryClusterFileReference> convertInputFiles(List<File> peaklistFilenames,
                                                               File binarySpectraDirectory)
            throws Exception {
        // if the binary spectra directory doesn't exist, create it
        if (!binarySpectraDirectory.exists()) {
            if (!binarySpectraDirectory.mkdir()) {
                throw new Exception("Failed to create temporary binary directory: " + binarySpectraDirectory);
            }
        }

        BinningSpectrumConverter binningSpectrumConverter = new BinningSpectrumConverter(binarySpectraDirectory,
                parallelJobs, useFastMode);

        // if verbose mode is enabled, send all progress updates to the progress listeners
        if (verbose) {
            for (IProgressListener progressListener : progressListeners) {
                binningSpectrumConverter.addProgressListener(progressListener);
            }
        }

        // tell the listeners that the conversion has started
        notifyProgressListeners(new ProgressUpdate(
                String.format("Converting %d input files...", peaklistFilenames.size()),
                ProgressUpdate.CLUSTERING_STAGE.CONVERSION));

        // convert to a String[] array - not a very nice way at the moment
        String[] peakListFilenameArray = new String[peaklistFilenames.size()];
        for (int i = 0; i < peaklistFilenames.size(); i++) {
            peakListFilenameArray[i] = peaklistFilenames.get(i).toString();
        }

        binningSpectrumConverter.processPeaklistFiles(peakListFilenameArray);

        // count the spectra per bin
        if (Defaults.getNumberOfComparisonAssessor().getClass() == SpectraPerBinNumberComparisonAssessor.class) {
            SpectraPerBinNumberComparisonAssessor assessor = (SpectraPerBinNumberComparisonAssessor) Defaults.getNumberOfComparisonAssessor();
            List<SpectrumReference> specRefs = binningSpectrumConverter.getSpectrumReferences();

            for (SpectrumReference specRef : specRefs) {
                assessor.countSpectrum(specRef.getPrecursorMz());
            }
        }

        return binningSpectrumConverter.getWrittenFiles();
    }

    /**
     * Find and index existing binary files. Binary files are identified by their file ending ".cls".
     * @param binarySpectraDirectory The directory in which the binary files should be searched for.
     * @return
     * @throws Exception
     */
    private List<BinaryClusterFileReference> findExistingBinaryFiles(File binarySpectraDirectory) throws Exception {
        // get the list of files
        File[] existingBinaryFiles = binarySpectraDirectory.listFiles(
                (FilenameFilter) FileFilterUtils.suffixFileFilter(".cls"));

        if (existingBinaryFiles == null) {
            throw new Exception("Failed to find binary files in " + binarySpectraDirectory.toString());
        }

        // index the existing files
        // tell the listeners that the scanning has started
        notifyProgressListeners(new ProgressUpdate(
                String.format("Scanning %d binary files...", existingBinaryFiles.length),
                ProgressUpdate.CLUSTERING_STAGE.CONVERSION));

        // count the spectra per bin while scanning the files
        SpectraPerBinNumberComparisonAssessor spectraPerBinNumberComparisonAssessor = null;
        if (Defaults.getNumberOfComparisonAssessor().getClass() == SpectraPerBinNumberComparisonAssessor.class) {
            spectraPerBinNumberComparisonAssessor = (SpectraPerBinNumberComparisonAssessor) Defaults.getNumberOfComparisonAssessor();
        }

        return BinaryFileScanner.scanBinaryFiles(spectraPerBinNumberComparisonAssessor, existingBinaryFiles);
    }

    /**
     * Clusters the passed binary files.
     * @param binaryFiles The input files to cluster in binary format.
     * @param clusteringThresholds The clustering thresholds to use. The thresholds are processed according
     *                             to the order of the list. For most similarity metrics a higher threshold
     *                             means a higher accuracy. In these cases the list should be sorted from
     *                             highest to lowest (clustering should start with the highest accuracy and
     *                             decrease with each subsequent clustering round).
     * @return A list of BinaryClusterFileReferenceS that represent the result files.
     * @throws Exception
     */
    private List<BinaryClusterFileReference> clusterFiles(List<BinaryClusterFileReference> binaryFiles,
                                                          List<Float> clusteringThresholds)
            throws Exception {
        // create the needed temporary directories
        File clusteringResultDirectory = createTemporaryDirectory("clustering_results", temporaryDirectory);
        File tmpClusteringResultDirectory = createTemporaryDirectory("clustering_results_tmp", temporaryDirectory);

        // cluster the files
        BinaryFileClusterer binaryFileClusterer = new BinaryFileClusterer(parallelJobs, clusteringResultDirectory,
                clusteringThresholds, useFastMode, tmpClusteringResultDirectory);

        // if verbose mode is enabled add the progress listeners to receive all updates
        if (verbose) {
            for (IProgressListener progressListener : progressListeners) {
                binaryFileClusterer.addProgressListener(progressListener);
            }
        }

        notifyProgressListeners(new ProgressUpdate(
                String.format("Clustering %d binary files...", binaryFiles.size()),
                ProgressUpdate.CLUSTERING_STAGE.CLUSTERING
        ));

        // start the clustering
        binaryFileClusterer.clusterFiles(binaryFiles);

        // delete the temporary directory
        if (!tmpClusteringResultDirectory.delete()) {
            // TODO: a notification could be added at a later stage.
        }
        // delete the binary input files if set
        if (!keepBinaryFiles) {
            for (BinaryClusterFileReference binaryFile : binaryFiles) {
                binaryFile.getResultFile().delete();
            }

            // delete the directory
            if (binaryFiles.size() > 0) {
                File parentDirectory = new File(binaryFiles.get(0).getResultFile().getParent());
                parentDirectory.delete();
            }
        }

        // sort result files by min m/z
        List<BinaryClusterFileReference> clusteredFiles = binaryFileClusterer.getResultFiles();
        clusteredFiles = new ArrayList<BinaryClusterFileReference>(clusteredFiles);
        Collections.sort(clusteredFiles);

        return clusteredFiles;
    }

    /**
     * Merge existing binary files that were created through a previous clustering
     * process.
     * @param binaryFilenames
     * @param thresholds
     * @param finalResultFile
     */
    public void mergeBinaryFiles(String[] binaryFilenames, List<Float> thresholds, File finalResultFile) throws Exception{
        // convert into an array of FileS
        File[] binaryFiles = new File[binaryFilenames.length];
        for (int i = 0; i < binaryFilenames.length; i++) {
            binaryFiles[i] = new File(binaryFilenames[i]);
        }

        // scan binary files
        notifyProgressListeners(new ProgressUpdate(
                String.format("Scanning %d binary files...", binaryFilenames.length),
                ProgressUpdate.CLUSTERING_STAGE.MERGING
        ));

        List<BinaryClusterFileReference> binaryFileReferences = BinaryFileScanner.scanBinaryFiles(binaryFiles);

        // merge the files
        File combinedResultFile = mergeClusteringResults(binaryFileReferences, thresholds);

        // create the output file
        convertCgfToClustering(combinedResultFile, finalResultFile, thresholds.get(thresholds.size() - 1));
    }

    /**
     * Merges the result files (clusters at the borders of two result files are re-processed based on the precursor
     * tolerance). The resulting clustering are automatically merged into a single CGF file.
     * @param clusteredFiles Clustered (binary) files to merge.
     * @param clusteringThresholds The clustering thresholds to use. The thresholds are processed according
     *                             to the order of the list. For most similarity metrics a higher threshold                             means a higher accuracy. In these cases the list should be sorted from                             highest to lowest (clustering should start with the highest accuracy and                             decrease with each subsequent clustering round).
     * @return The resulting CGF file as a File object.
     * @throws Exception
     */
    private File mergeClusteringResults(List<BinaryClusterFileReference> clusteredFiles, List<Float> clusteringThresholds)
            throws Exception {
        if (clusteredFiles.size() < 1) {
            throw new Exception("No clustering result files found for merging.");
        }

        // create the required temporary directories
        File mergedResultsDirectory = createTemporaryDirectory("merged_results", temporaryDirectory);
        File mergedResultsDirectoryTmp = createTemporaryDirectory("merged_results_tmp", temporaryDirectory);

        // create the merger
        BinaryFileMergingClusterer mergingClusterer = new BinaryFileMergingClusterer(parallelJobs, mergedResultsDirectory,
                clusteringThresholds, useFastMode, Defaults.getDefaultPrecursorIonTolerance(), deleteTemporaryFiles,
                mergedResultsDirectoryTmp);

        // if verbose mode is enabled add the progress listeners to receive all updates
        if (verbose) {
            for (IProgressListener progressListener : progressListeners) {
                mergingClusterer.addProgressListener(progressListener);
            }
        }

        // create the combined output file as soon as a job is done
        File combinedResultFile = File.createTempFile("combined_clustering_results", ".cgf", temporaryDirectory);
        MergingCGFConverter mergingCGFConverter = new MergingCGFConverter(combinedResultFile, deleteTemporaryFiles);
        mergingClusterer.addListener(mergingCGFConverter);

        // launch the merging job
        notifyProgressListeners(new ProgressUpdate(
                String.format("Merging %d binary files...", clusteredFiles.size()),
                ProgressUpdate.CLUSTERING_STAGE.MERGING
        ));

        mergingClusterer.clusterFiles(clusteredFiles);

        // delete the temporary directory after merging - this directory should be empty
        if (!mergedResultsDirectoryTmp.delete()) {
            // TODO: add notification at a later stage
        }

        // delete the temporary directories if set
        if (deleteTemporaryFiles) {
            File clusteredFilesDirectory = clusteredFiles.get(0).getResultFile().getParentFile();

            if (!clusteredFilesDirectory.delete()) {
                // TODO: add notification at a later stage
            }

            if (!mergedResultsDirectory.delete()) {
                // TODO: add notification at a later stage
            }
        }

        return combinedResultFile;
    }

    /**
     * Convert a cgf file to a .clustering file.
     * @param cgfResultFile The .cgf input file to convert.
     * @param clusteringOutputFile File to write the .clustering file to. This file will be overwritten if it exists.
     * @param targetThreshold The final threshold used during the clustering process. This will be written as
     *                        as metadata in the header of the .clustering file.
     * @throws Exception
     */
    public void convertCgfToClustering(File cgfResultFile, File clusteringOutputFile, float targetThreshold) throws Exception {
        FileInputStream fileInputStream = new FileInputStream(cgfResultFile);
        FileWriter writer = new FileWriter(clusteringOutputFile);

        int nClusterWritten = 0;

        notifyProgressListeners(new ProgressUpdate("Converting results to .clustering format...",
                ProgressUpdate.CLUSTERING_STAGE.OUTPUT));

        DotClusterClusterAppender.INSTANCE.appendStart(writer, "GreedyClustering_" + String.valueOf(targetThreshold));
        CGFSpectrumIterable iterable = new CGFSpectrumIterable(fileInputStream);

        for (ICluster cluster : iterable) {
            DotClusterClusterAppender.INSTANCE.appendCluster(writer, cluster);
            nClusterWritten++;

            if (nClusterWritten % 10000 == 0 && isVerbose()) {
                notifyProgressListeners(new ProgressUpdate(
                        String.format("%d clusters successfully written to %s", nClusterWritten, clusteringOutputFile.toString()),
                        ProgressUpdate.CLUSTERING_STAGE.OUTPUT));
            }
        }

        writer.close();

        notifyProgressListeners(new ProgressUpdate(
                String.format("%d clusters successfully written to %s", nClusterWritten, clusteringOutputFile.toString()),
                ProgressUpdate.CLUSTERING_STAGE.OUTPUT));
    }

    /**
     * Helper function to create the List of clustering thresholds based on the initial starting threshold
     * (highest precision), the target threshold, and the number of clustering rounds to perform.
     * @param startThreshold Highest starting precision to use.
     * @param endThreshold Lowest target precision for the clustering run.
     * @param clusteringRounds Number of round of clustering to perform.
     * @return A List of thresholds used as parameter for the clustering runs.
     */
    public List<Float> generateClusteringThresholds(Float startThreshold, Float endThreshold, int clusteringRounds) {
        List<Float> thresholds = new ArrayList<Float>(clusteringRounds);
        float stepSize = (startThreshold - endThreshold) / (clusteringRounds - 1);

        for (int i = 0; i < clusteringRounds; i++) {
            thresholds.add(startThreshold - (stepSize * i));
        }

        return thresholds;
    }

    /**
     * Send a ProgressUpdate to all progress listeners.
     * @param progressUpdate
     */
    private void notifyProgressListeners(ProgressUpdate progressUpdate) {
        for (IProgressListener progressListener : progressListeners) {
            progressListener.onProgressUpdate(progressUpdate);
        }
    }
}
