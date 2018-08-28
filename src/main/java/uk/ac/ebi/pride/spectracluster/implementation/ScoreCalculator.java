package uk.ac.ebi.pride.spectracluster.implementation;

import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.similarity.ISimilarityChecker;
import uk.ac.ebi.pride.spectracluster.spectra_list.ParsingMgfScanner;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.ebi.pride.spectracluster.util.SpectrumConverter;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Stream;

/**
 * This class is used to calculate the similarity score of every spectrum
 * with the cluster's consensus spectrum. Since the consensus spectrum
 * changes during clustering, this can only be done afterwards.
 *
 * Created by jg on 03.05.17.
 */
public class ScoreCalculator {
    public void processClusteringResult(File clusteringFile, File[] mgfPaths) throws Exception {
        // make sure the input file can be overwritten
        if (!clusteringFile.canWrite()) {
            throw new Exception("Cannot write to " + clusteringFile.getName());
        }

        // make sure all MGF files exist and get their paths
        Map<String, Path> mgfPathMap = findMgfFiles(clusteringFile, mgfPaths);

        // index all MGF files
        Map<String, List<IndexElement>> mgfFileIndices = indexMgfFiles(mgfPathMap);

        // add the scores
        Path tmpFile = Files.createTempFile(null, ".clustering");

        addSimilarityScores(clusteringFile.toPath(), tmpFile, mgfPathMap, mgfFileIndices);

        // replace the original file
        Files.move(tmpFile, clusteringFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
    }

    /**
     * Add the similarity scores to the .clustering file.
     *
     * @param clusteringIn .clustering file to create the similarity scores for
     * @param clusteringOut Path to the newly created .clustering file.
     * @param pathMap A Map with the MGF filename as key and the complete Path as value
     * @param fileIndices A Map with the MGF filename as key and the List of Indices as value
     * @throws Exception
     */
    protected void addSimilarityScores(Path clusteringIn, Path clusteringOut, Map<String, Path> pathMap, Map<String,
            List<IndexElement>> fileIndices) throws Exception {
        // process the .clustering file line by line
        try (BufferedReader reader = Files.newBufferedReader(clusteringIn);
             BufferedWriter writer = Files.newBufferedWriter(clusteringOut)) {
            // initialize consensus spec variables
            float precursor = 0F;
            List<Float> mz = null;
            List<Float> intens = null;
            ISpectrum consensus = null;

            ISimilarityChecker similarityChecker = Defaults.getDefaultSimilarityChecker();

            // process the file line by line
            String line = null;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith("av_precursor_mz=")) {
                    precursor = Float.parseFloat(line.substring(16));
                }
                if (line.startsWith("consensus_mz=")) {
                    mz = parseFloatString(line.substring(13));
                }
                if (line.startsWith("consensus_intens=")) {
                    intens = parseFloatString(line.substring(17));
                    consensus = createConsensusSpectrum(precursor, mz, intens);
                }

                if (line.startsWith("SPEC\t")) {
                    ISpectrum spectrum = loadSpectrum(line, pathMap, fileIndices);

                    double score = similarityChecker.assessSimilarity(consensus, spectrum);

                    // replace the respective line
                    line = line.substring(0, ordinalIndexOf(line, "\t", 8) + 1) +
                            String.valueOf(score) + "\n";
                }

                writer.write(line + "\n");
            }
        }
    }

    /**
     * Finds the nth occurrence of a substring
     * @param str The string in which to search for.
     * @param substr The substring to search for
     * @param n The nth occurrence.
     * @return The index.
     */
    public static int ordinalIndexOf(String str, String substr, int n) {
        int pos = str.indexOf(substr);
        while (--n > 0 && pos != -1)
            pos = str.indexOf(substr, pos + 1);
        return pos;
    }

    protected ISpectrum loadSpectrum(String line, Map<String, Path> pathMap, Map<String, List<IndexElement>> fileIndices) throws Exception {
        // extract the filename and spec id from the line
        int idIndex = line.indexOf("#id=");

        if (idIndex < 0) {
            throw new Exception("Invalid spectrum line encountered: " + line);
        }

        String filePath = line.substring(11, idIndex);
        // idIndex + length of "#id=" + length of "index="
        String specIdString = line.substring(idIndex + 10, line.indexOf("#title="));

        int specId = Integer.parseInt(specIdString);

        Path file = Paths.get(filePath);
        String filename = file.getFileName().toString();

        if (!pathMap.containsKey(filename) || !fileIndices.containsKey(filename)) {
            throw new IllegalStateException("Missing " + filename + " in internal maps");
        }

        JMzReader peakListReader = new MgfFile(pathMap.get(filename).toFile(),
                fileIndices.get(filename), true);
        uk.ac.ebi.pride.tools.jmzreader.model.Spectrum spectrum = peakListReader.getSpectrumByIndex(specId);
        ISpectrum convertedSpectrum = SpectrumConverter.convertJmzReaderSpectrum(spectrum, "", filename);

        ISpectrum processedSpectrum = ClusteringSettings.getInitialSpectrumFilter().apply(convertedSpectrum);
        // normalize the spectrum
        processedSpectrum = new Spectrum(processedSpectrum,
                ClusteringSettings.getIntensityNormalizer().normalizePeaks(processedSpectrum.getPeaks()));

        // process the spectrum by creating a consensus spectrum
        IConsensusSpectrumBuilder consensusSpectrumBuilder = Defaults.getDefaultConsensusSpectrumBuilder();
        consensusSpectrumBuilder.addSpectra(processedSpectrum);

        ISpectrum consensusSpectrum = consensusSpectrumBuilder.getConsensusSpectrum();

        // only retain the peaks for comparison
        ISpectrum filteredSpectrum = new Spectrum(consensusSpectrum,
                ClusteringSettings.getComparisonFilterFunction().apply(consensusSpectrum.getPeaks()));

        return filteredSpectrum;
    }

    /**
     * Create a cluster representing the consensus spectrum based on the
     * precursor m/z value and the m/z and intensity values.
     *
     * @param precursor The precursor m/z.
     * @param mz List of m/z values.
     * @param intens List of intensity values.
     * @return
     */
    private ISpectrum createConsensusSpectrum(float precursor, List<Float> mz, List<Float> intens) {
        List<IPeak> peakList = new ArrayList<>(mz.size());

        for (int i = 0; i < mz.size(); i++) {
            peakList.add(new Peak(mz.get(i), intens.get(i)));
        }

        ISpectrum spec = new Spectrum("", 0, precursor, Defaults.getDefaultQualityScorer(),
                peakList);

        ISpectrum filteredSpectrum = new Spectrum(spec, ClusteringSettings.getComparisonFilterFunction().apply(spec.getPeaks()));

        return filteredSpectrum;
    }

    /**
     * Convert a string of "," delimited floats into a List of Floats.
     *
     * @param string String to parse.
     * @return List of Floats
     */
    protected List<Float> parseFloatString(String string) {
        StringTokenizer tokenizer = new StringTokenizer(string, ",");

        List<Float> floats = new ArrayList<>();
        while(tokenizer.hasMoreTokens()) {
            floats.add(Float.parseFloat(tokenizer.nextToken()));
        }

        return floats;
    }

    /**
     * Indexes the passed MGF files and returns a Map containing the
     * MGF filename as key and the list of indexes as values.
     * @param mgfPathMap Map with the MGF filenames as key and the Path as value
     * @return Map with the MGF filename as key and the List of Indexes as value
     * @throws Exception In case parsing fails.
     */
    protected Map<String, List<IndexElement>> indexMgfFiles(Map<String, Path> mgfPathMap) throws Exception {
        Map<String, List<IndexElement>> fileIndices = new HashMap<>();

        for (String mgfFilename : mgfPathMap.keySet()) {

            ParsingMgfScanner scanner = new ParsingMgfScanner();
            String[] filename = { mgfPathMap.get(mgfFilename).toString() };
            scanner.getSpectrumReferences(filename);

            if (scanner.getFileIndices().size() != 1) {
                throw new IllegalStateException("ParsingMgfScanner contains multiple indices.");
            }

            fileIndices.put(mgfFilename, scanner.getFileIndices().get(0));
        }

        return fileIndices;
    }

    /**
     * Tests whether all MGF files specified in the .clustering result
     * file are present in the supplied path. If a file is missing, an
     * Exception is thrown. Returns a Map with the MGF filename as key and
     * the matching Path object as value.
     *
     * @param clusteringFile Path to the clustering file to process
     * @param mgfPaths An array of directories where to find the MGF files.
     * @return A Map with the MGF filenames as keys and the matching Path objects
     *         as values.
     * @throws Exception Thrown if an MGF file cannot be found.
     */
    protected Map<String, Path> findMgfFiles(File clusteringFile, File[] mgfPaths) throws Exception {
        Set<String> mgfFilenames = extractMgfFilenames(clusteringFile);
        Map<String, Path> pathPerFilename = new HashMap<>(mgfFilenames.size());

        for (String mgfFilename : mgfFilenames) {
            boolean fileFound = false;

            for (File mgfPath : mgfPaths) {
                Path file = mgfPath.toPath().resolve(mgfFilename);

                if (Files.isRegularFile(file)) {
                    fileFound = true;
                    pathPerFilename.put(mgfFilename, file);
                    break;
                }
            }

            if (!fileFound) {
                throw new Exception("Cannot find " + mgfFilename + " in the supplied directories.");
            }
        }

        return pathPerFilename;
    }

    /**
     * Extracts the filenames of all MGF files referenced in the
     * clustering file.
     * @param clusteringFile The clustering file to process.
     */
    protected Set<String> extractMgfFilenames(File clusteringFile) throws Exception {
        Path path = clusteringFile.toPath();
        Set<String> mgfFilenames = new HashSet<>();

        try (Stream<String> lines = Files.lines(path).filter(s -> s.contains("SPEC\t"))) {
            lines.forEach((String line) -> {
                String[] fields = line.split("\t");

                if (fields.length < 5) {
                    throw new IllegalStateException("Unexpected SPEC line encountered: " + line);
                }

                String idField = fields[1];

                String filename = idField.substring(idField.indexOf("#file=") + 6,
                        idField.indexOf("#id="));
                Path filePath = Paths.get(filename);
                mgfFilenames.add(filePath.getFileName().toString());
            });
        }

        return mgfFilenames;
    }
}
