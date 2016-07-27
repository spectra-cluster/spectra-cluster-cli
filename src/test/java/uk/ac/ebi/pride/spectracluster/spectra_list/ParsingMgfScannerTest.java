package uk.ac.ebi.pride.spectracluster.spectra_list;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.io.File;
import java.util.List;

/**
 * Created by jg on 27.07.16.
 */
public class ParsingMgfScannerTest {
    private File testFile;

    @Before
    public void setUp() throws Exception {
        testFile = new File(ParsingMgfScannerTest.class.getClassLoader().getResource("header_test.mgf").toURI());
    }

    @Test
    public void testFileWithHeader() throws Exception {
        ParsingMgfScanner scanner = new ParsingMgfScanner();

        String[] filenames = {testFile.getAbsolutePath()};
        List<SpectrumReference> spectrumReferences = scanner.getSpectrumReferences(filenames);

        Assert.assertEquals(3, spectrumReferences.size());

        List<IndexElement> fileIndex = scanner.getFileIndices().get(0);
        Assert.assertEquals(3, fileIndex.size());

        Assert.assertEquals(43, fileIndex.get(0).getStart());
    }
}
