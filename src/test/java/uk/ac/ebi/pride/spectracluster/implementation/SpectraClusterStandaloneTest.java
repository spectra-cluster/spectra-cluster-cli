package uk.ac.ebi.pride.spectracluster.implementation;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cdf.SpectraPerBinNumberComparisonAssessor;
import uk.ac.ebi.pride.spectracluster.util.Defaults;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class SpectraClusterStandaloneTest {
    private List<File> testMgfFiles;

    @Before
    public void setUp() throws Exception {
        File testMgf = new File(SpectraClusterStandaloneTest.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        testMgfFiles = new ArrayList<>(1);
        testMgfFiles.add(testMgf);
    }

    @Test
    public void testClusteringKeepBinary() throws Exception {
        SpectraClusterStandalone spectraClusterStandalone = new SpectraClusterStandalone();

        // set all clustering parameters
        Defaults.setDefaultPrecursorIonTolerance(0.5f);
        Defaults.setFragmentIonTolerance(0.5f);
        List<Float> thresholds = spectraClusterStandalone.generateClusteringThresholds(1.0f, 0.99f, 3);

        spectraClusterStandalone.setKeepBinaryFiles(true);
        spectraClusterStandalone.setDeleteTemporaryFiles(false);
        spectraClusterStandalone.setParallelJobs(1);
        ClusteringSettings.addIntitalSpectrumFilter(ClusteringSettings.SPECTRUM_FILTER.MZ_150.filter);

        Defaults.setNumberOfComparisonAssessor(new SpectraPerBinNumberComparisonAssessor(Defaults.getDefaultPrecursorIonTolerance()));

        File resultFile = File.createTempFile("clustering_test_", ".clustering");
        resultFile.deleteOnExit();
        spectraClusterStandalone.clusterPeaklistFiles(testMgfFiles, thresholds, resultFile);

        Assert.assertTrue(resultFile.exists());
    }

    @Test
    public void testClusteringDeleteBinary() throws Exception {
        SpectraClusterStandalone spectraClusterStandalone = new SpectraClusterStandalone();

        // set all clustering parameters
        Defaults.setDefaultPrecursorIonTolerance(0.5f);
        Defaults.setFragmentIonTolerance(0.5f);
        List<Float> thresholds = spectraClusterStandalone.generateClusteringThresholds(1.0f, 0.99f, 3);

        spectraClusterStandalone.setKeepBinaryFiles(false);
        spectraClusterStandalone.setDeleteTemporaryFiles(false);
        spectraClusterStandalone.setParallelJobs(1);
        ClusteringSettings.addIntitalSpectrumFilter(ClusteringSettings.SPECTRUM_FILTER.MZ_150.filter);

        Defaults.setNumberOfComparisonAssessor(new SpectraPerBinNumberComparisonAssessor(Defaults.getDefaultPrecursorIonTolerance()));

        File resultFile = File.createTempFile("clustering_test_", ".clustering");
        resultFile.deleteOnExit();
        spectraClusterStandalone.clusterPeaklistFiles(testMgfFiles, thresholds, resultFile);

        Assert.assertTrue(resultFile.exists());
    }

    @Test
    public void testClusteringDeleteAll() throws Exception {
        SpectraClusterStandalone spectraClusterStandalone = new SpectraClusterStandalone();

        // set all clustering parameters
        Defaults.setDefaultPrecursorIonTolerance(0.5f);
        Defaults.setFragmentIonTolerance(0.5f);
        List<Float> thresholds = spectraClusterStandalone.generateClusteringThresholds(1.0f, 0.99f, 3);

        spectraClusterStandalone.setKeepBinaryFiles(false);
        spectraClusterStandalone.setDeleteTemporaryFiles(true);
        spectraClusterStandalone.setParallelJobs(1);
        ClusteringSettings.addIntitalSpectrumFilter(ClusteringSettings.SPECTRUM_FILTER.MZ_150.filter);

        Defaults.setNumberOfComparisonAssessor(new SpectraPerBinNumberComparisonAssessor(Defaults.getDefaultPrecursorIonTolerance()));

        File resultFile = File.createTempFile("clustering_test_", ".clustering");
        resultFile.deleteOnExit();
        spectraClusterStandalone.clusterPeaklistFiles(testMgfFiles, thresholds, resultFile);

        Assert.assertTrue(resultFile.exists());
    }
}
