package uk.ac.ebi.pride.spectracluster.binning;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.implementation.ClusteringSettings;
import uk.ac.ebi.pride.spectracluster.io.BinaryClusterIterable;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Created by jg on 12.10.17.
 */
public class BinningSpectrumConverterTests {
    private File mgfTestFile;
    private File clusteringTestFile;

    @Before
    public void setUp() throws Exception {
        mgfTestFile = new File(BinningSpectrumConverterTests.class.getClassLoader().getResource("spectra_400.0_4.0.mgf").toURI());
        clusteringTestFile = new File(BinningSpectrumConverterTests.class.getClassLoader().getResource("imp_hela_test_API11.clustering").toURI());
    }

    @Test
    public void testReferenceMzBinner() throws Exception {
        String[] filenames = {mgfTestFile.getAbsolutePath(), clusteringTestFile.getAbsolutePath()};

        // create temporary directory
        Path tmpDir = Files.createTempDirectory("spectra_cluster_cli_test");

        BinningSpectrumConverter spectrumConverter = new BinningSpectrumConverter(tmpDir.toFile(), 1, false);
        spectrumConverter.processPeaklistFiles(filenames);

        // read all created files
        int nClusters = 0;
        int nFiles = 0;
        File[] binaryResultFiles = tmpDir.toFile().listFiles();

        for (File binaryResultFile : binaryResultFiles) {
            // parse all clusters
            ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(binaryResultFile)));
            BinaryClusterIterable iterable = new BinaryClusterIterable(inputStream);

            for (ICluster cluster : iterable) {
                nClusters++;
            }

            nFiles++;

            // clean up
            binaryResultFile.delete();
        }

        Assert.assertEquals(40, nFiles);
        Assert.assertEquals(7397, nClusters);
        tmpDir.toFile().delete();
    }

    @Test
    public void testIdentifiedOnlyBinner() throws Exception {
            String[] filenames = {mgfTestFile.getAbsolutePath(), clusteringTestFile.getAbsolutePath()};

            // create temporary directory
            Path tmpDir = Files.createTempDirectory("spectra_cluster_cli_test");

            BinningSpectrumConverter spectrumConverter = new BinningSpectrumConverter(tmpDir.toFile(), 1, false);
            spectrumConverter.setLoadingMode(ClusteringSettings.LOADING_MODE.ONLY_IDENTIFIED);

            spectrumConverter.processPeaklistFiles(filenames);

            // read all created files
            int nClusters = 0;
            int nFiles = 0;
            File[] binaryResultFiles = tmpDir.toFile().listFiles();

            for (File binaryResultFile : binaryResultFiles) {
                // parse all clusters
                ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(binaryResultFile)));
                BinaryClusterIterable iterable = new BinaryClusterIterable(inputStream);

                for (ICluster cluster : iterable) {
                    nClusters++;
                }

                nFiles++;

                // clean up
                binaryResultFile.delete();
            }

            Assert.assertEquals(40, nFiles);
            Assert.assertEquals(7187, nClusters);
            tmpDir.toFile().delete();
    }

    @Test
    public void testUnidentifiedOnlyBinner() throws Exception {
        String[] filenames = {mgfTestFile.getAbsolutePath(), clusteringTestFile.getAbsolutePath()};

        // create temporary directory
        Path tmpDir = Files.createTempDirectory("spectra_cluster_cli_test");

        BinningSpectrumConverter spectrumConverter = new BinningSpectrumConverter(tmpDir.toFile(), 1, false);
        spectrumConverter.setLoadingMode(ClusteringSettings.LOADING_MODE.ONLY_UNIDENTIFIED);

        spectrumConverter.processPeaklistFiles(filenames);

        // read all created files
        int nClusters = 0;
        int nFiles = 0;
        File[] binaryResultFiles = tmpDir.toFile().listFiles();

        for (File binaryResultFile : binaryResultFiles) {
            // parse all clusters
            ObjectInputStream inputStream = new ObjectInputStream(new BufferedInputStream(new FileInputStream(binaryResultFile)));
            BinaryClusterIterable iterable = new BinaryClusterIterable(inputStream);

            for (ICluster cluster : iterable) {
                nClusters++;
            }

            nFiles++;

            // clean up
            binaryResultFile.delete();
        }

        Assert.assertEquals(40, nFiles);
        Assert.assertEquals(7187, nClusters);
        tmpDir.toFile().delete();
    }
}
