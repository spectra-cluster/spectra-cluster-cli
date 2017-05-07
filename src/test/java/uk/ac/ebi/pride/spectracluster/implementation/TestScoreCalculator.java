package uk.ac.ebi.pride.spectracluster.implementation;

import junit.framework.Assert;
import junit.framework.TestCase;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveIonContaminantsPeaksFunction;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.RemoveWindowPeaksFunction;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by jg on 03.05.17.
 */
public class TestScoreCalculator extends TestCase {
    private File testFile;
    private File newTestFile;

    @Before
    public void setUp() throws Exception {
        testFile = new File(TestScoreCalculator.class.getClassLoader().getResource("test.clustering").toURI());
        newTestFile = new File(TestScoreCalculator.class.getClassLoader().getResource("test_new.clustering").toURI());
    }

    @Test
    public void testGetMgfFiles() throws Exception {
        ScoreCalculator scoreCalculator = new ScoreCalculator();

        Set<String> mgfFilenames = scoreCalculator.extractMgfFilenames(testFile);

        Assert.assertEquals(2, mgfFilenames.size());
        Assert.assertTrue(mgfFilenames.contains("spectra_400.0_4.0.mgf"));
        Assert.assertTrue(mgfFilenames.contains("F001257.mgf"));
    }

    @Test
    public void testFindMgfFiles() throws Exception {
        ScoreCalculator scoreCalculator = new ScoreCalculator();

        // get the path by using the parent of the test file
        File[] mgfFileDirs = {new File(testFile.getParent())};

        Map<String, Path> filePaths = scoreCalculator.findMgfFiles(testFile, mgfFileDirs);

        Assert.assertEquals(2, filePaths.size());
        Assert.assertTrue(filePaths.containsKey("spectra_400.0_4.0.mgf"));
        Assert.assertTrue(filePaths.containsKey("F001257.mgf"));
    }

    @Test
    public void testIndexMgfFiles() throws Exception {
        ScoreCalculator scoreCalculator = new ScoreCalculator();

        // get the path by using the parent of the test file
        File[] mgfFileDirs = {new File(testFile.getParent())};

        Map<String, Path> filePaths = scoreCalculator.findMgfFiles(testFile, mgfFileDirs);

        Map<String, List<IndexElement>> indexes = scoreCalculator.indexMgfFiles(filePaths);

        Assert.assertEquals(2, indexes.size());
        Assert.assertTrue(indexes.containsKey("spectra_400.0_4.0.mgf"));
        Assert.assertTrue(indexes.containsKey("F001257.mgf"));
        Assert.assertEquals(210, indexes.get("spectra_400.0_4.0.mgf").size());
        Assert.assertEquals(10, indexes.get("F001257.mgf").size());
    }

    @Test
    public void testAddSimilarityScores() throws Exception {
        ScoreCalculator scoreCalculator = new ScoreCalculator();

        // get the path by using the parent of the test file
        File[] mgfFileDirs = {new File(testFile.getParent())};

        Map<String, Path> filePaths = scoreCalculator.findMgfFiles(testFile, mgfFileDirs);
        Map<String, List<IndexElement>> indexes = scoreCalculator.indexMgfFiles(filePaths);

        ClusteringSettings.addIntitalSpectrumFilter(new RemoveWindowPeaksFunction(200.0F, Float.MAX_VALUE));

        Path tmpFile = Files.createTempFile("debug_", ".clustering");

        scoreCalculator.addSimilarityScores(testFile.toPath(), tmpFile, filePaths, indexes);

        // test the first SPEC line
        Stream<String> lines = Files.lines(tmpFile);
        Optional<String> isSpecLine = lines.filter(s -> s.startsWith("SPEC\t")).findFirst();

        if (isSpecLine.isPresent()) {
            String line = isSpecLine.get();

            String[] fields = line.split("\t");
            Assert.assertEquals(9, fields.length);
            Assert.assertEquals("111.58301981855692", fields[8]);
        }

        // get all score fields
        lines = Files.lines(tmpFile);
        List<String> scoreStrings = lines.filter(s -> s.startsWith("SPEC\t")).map(s -> s.split("\t")[8]).collect(Collectors.toList());

        for (String scoreString : scoreStrings) {
            Assert.assertFalse("Invalid score", scoreString.equals("NaN"));

            try {
                Double score = Double.parseDouble(scoreString);
                Assert.assertNotNull(score);
            }
            catch (Exception e) {
                Assert.fail("Invalid score encountered: " + scoreString);
            }
        }
    }

    /**
     * Tests the scores on the clustering output created
     * with the adapted consensus spectrum builder.
     * @throws Exception
     */
    @Test
    public void testAddSimilarityScoresNew() throws Exception {
        ScoreCalculator scoreCalculator = new ScoreCalculator();

        // get the path by using the parent of the test file
        File[] mgfFileDirs = {new File(newTestFile.getParent())};

        Map<String, Path> filePaths = scoreCalculator.findMgfFiles(newTestFile, mgfFileDirs);
        Map<String, List<IndexElement>> indexes = scoreCalculator.indexMgfFiles(filePaths);

        ClusteringSettings.addIntitalSpectrumFilter(new RemoveWindowPeaksFunction(200.0F, Float.MAX_VALUE));

        Path tmpFile = Files.createTempFile("debug_", ".clustering");

        scoreCalculator.addSimilarityScores(newTestFile.toPath(), tmpFile, filePaths, indexes);

        // test the first SPEC line
        Stream<String> lines = Files.lines(tmpFile);
        Optional<String> isSpecLine = lines.filter(s -> s.startsWith("SPEC\t")).findFirst();

        if (isSpecLine.isPresent()) {
            String line = isSpecLine.get();

            String[] fields = line.split("\t");
            Assert.assertEquals(9, fields.length);
            Assert.assertEquals("111.58301981855692", fields[8]);
        }

        // get all score fields
        lines = Files.lines(tmpFile);
        List<String> scoreStrings = lines.filter(s -> s.startsWith("SPEC\t")).map(s -> s.split("\t")[8]).collect(Collectors.toList());

        for (String scoreString : scoreStrings) {
            Assert.assertFalse("Invalid score", scoreString.equals("NaN"));

            try {
                Double score = Double.parseDouble(scoreString);
                Assert.assertNotNull(score);
            }
            catch (Exception e) {
                Assert.fail("Invalid score encountered: " + scoreString);
            }
        }
    }

    @Test
    public void testParseFloatString() throws Exception {
        String testString = "1,2,3.456,4.564\n";
        ScoreCalculator scoreCalculator = new ScoreCalculator();

        List<Float> floats = scoreCalculator.parseFloatString(testString);

        Assert.assertEquals(4, floats.size());
        Assert.assertEquals(1.0F, floats.get(0));
        Assert.assertEquals(3.456F, floats.get(2));
        Assert.assertEquals(4.564F, floats.get(3));
    }
}
