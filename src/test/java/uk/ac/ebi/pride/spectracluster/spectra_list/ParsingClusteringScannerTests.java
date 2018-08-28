package uk.ac.ebi.pride.spectracluster.spectra_list;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.List;

/**
 * Created by jg on 13.10.17.
 */
public class ParsingClusteringScannerTests {
    private File clusteringTestFile;

    @Before
    public void setUp() throws Exception {
        clusteringTestFile = new File(ParsingClusteringScannerTests.class.getClassLoader().getResource("test.clustering").toURI());
    }

    @Test
    public void testGetSpectrumReferences() throws Exception {
        ParsingClusteringScanner scanner = new ParsingClusteringScanner();
        String[] filenames = {clusteringTestFile.getAbsolutePath()};
        List<SpectrumReference> references = scanner.getSpectrumReferences(filenames);

        Assert.assertEquals(197, references.size());
        for (SpectrumReference reference : references) {
            Assert.assertEquals(0, reference.getFileId());
            Assert.assertEquals(SpectrumReference.IS_CLUSTER, reference.getSpectrumIndex());
        }

        Assert.assertEquals(197, scanner.getClusteringFileIndices().get(0).getIndex().size());
    }
}
