package uk.ac.ebi.pride.spectracluster.spectra_list;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.clustering.BinaryClusterFileReference;
import uk.ac.ebi.pride.spectracluster.util.BinaryFileScanner;

import java.io.File;
import java.util.List;

/**
 * Created by jg on 14.01.18.
 */
public class BinaryFileScannerTest {
    private File testFile;

    @Before
    public void setUp() throws Exception {
        testFile = new File(BinaryFileScannerTest.class.getClassLoader().getResource("merged_0200.cls").toURI());
    }

    @Test
    public void scanAll() throws Exception {
        List<BinaryClusterFileReference> references = BinaryFileScanner.scanBinaryFiles(testFile);

        Assert.assertEquals(1, references.size());

        BinaryClusterFileReference fileReference = references.get(0);
        Assert.assertEquals(11615, fileReference.getnSpectra());
    }

    @Test
    public void scanIdentified() throws Exception {
        List<BinaryClusterFileReference> references = BinaryFileScanner.scanBinaryFiles(null,
                BinaryFileScanner.ONLY_IDENTIFIED_LOADING_PREDICATE, testFile);

        Assert.assertEquals(1, references.size());

        BinaryClusterFileReference fileReference = references.get(0);
        Assert.assertEquals(0, fileReference.getnSpectra());
    }

    @Test
    public void scanUnidentified() throws Exception {
        List<BinaryClusterFileReference> references = BinaryFileScanner.scanBinaryFiles(null,
                BinaryFileScanner.ONLY_UNIDENTIFIED_LOADING_PREDICATE, testFile);

        Assert.assertEquals(1, references.size());

        BinaryClusterFileReference fileReference = references.get(0);
        Assert.assertEquals(11615, fileReference.getnSpectra());
    }
}
