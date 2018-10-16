package uk.ac.ebi.pride.spectracluster.binning;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.util.BinaryFileScanner;

import java.io.File;
import java.util.List;
import java.util.Objects;

public class BinaryFileBinnerTest {
    File testFile;
    List<BinaryClusterFileReference> references;
    File outputDir;

    @Before
    public void setUp() throws Exception {
        testFile = getResource("convertedSpectra_0406_0415.cls");
        File[] testFiles = {testFile};
        references = BinaryFileScanner.scanBinaryFiles(testFiles);

        File tmpFile = File.createTempFile("rebinner_test", "");
        if (!tmpFile.delete()) {
            throw new Exception("Failed to delete temporary file");
        }
        if (!tmpFile.mkdir()) {
            throw new Exception("Failed to create temporary directory");
        }
        outputDir = tmpFile;
    }

    @Test
    public void testRebinningSingleFile() throws Exception {
        Assert.assertEquals(1, references.size());
        List<BinaryClusterFileReference> rebinned = BinaryFileRebinner.rebinBinaryFiles(references, outputDir, 0.01);
        Assert.assertEquals(1, rebinned.size());
    }

    @Test
    public void testRebinning() throws Exception {
        File file1 = getResource("clustered_binary/convertedSpectra_0350_0375.cls");
        File file2 = getResource("clustered_binary/convertedSpectra_0375_0385.cls");
        File[] testFiles = {file1, file2};
        List<BinaryClusterFileReference> references = BinaryFileScanner.scanBinaryFiles(testFiles);
        Assert.assertEquals(2, references.size());

        List<BinaryClusterFileReference> rebinned = BinaryFileRebinner.rebinBinaryFiles(references, outputDir, 2.0);
        Assert.assertEquals(2, rebinned.size());
        // next file always has to start 1 precursor tolerance lower
        Assert.assertTrue(references.get(0).getMaxMz() - 2.0 >= rebinned.get(0).getMaxMz());
        Assert.assertTrue(references.get(1).getMaxMz() - 2.0 >= rebinned.get(1).getMaxMz());
        Assert.assertTrue(references.get(1).getMaxMz() - 2.0 <= rebinned.get(2).getMinMz());

        // make sure the number of spectra stayed the same
        int specIn = references.stream().mapToInt(BinaryClusterFileReference::getnSpectra).sum();
        int specOut = rebinned.stream().mapToInt(BinaryClusterFileReference::getnSpectra).sum();
        Assert.assertEquals(specIn, specOut);
    }

    @After
    public void tearDown() {
        outputDir.delete();
    }

    private File getResource(String name) throws Exception {
        return new File(Objects.requireNonNull(BinaryFileBinnerTest.class.getClassLoader().getResource(name)).toURI());
    }
}
