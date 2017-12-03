package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.tools.braf.BufferedRandomAccessFile;
import uk.ac.ebi.pride.tools.jmzreader.model.IndexElement;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.IndexElementImpl;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by jg on 28.05.15.
 */
public class ParsingMgfScanner implements IPeaklistScanner {
    private List<List<IndexElement>> fileIndices;
    private boolean ignoreEmptySpectra;

    public ParsingMgfScanner(boolean ignoreEmptySpectra) {
        this.ignoreEmptySpectra = ignoreEmptySpectra;
    }

    public ParsingMgfScanner() {
        this(true);
    }

    @Override
    public Map<Integer, List<SpectrumReference>> getSpectraPerMajorPeaks(String[] filenames, int nMajorPeaks) throws Exception {
        throw new UnsupportedOperationException();
    }

    @Override
    public List<SpectrumReference> getSpectrumReferences(String[] filenames) throws Exception {
        List<SpectrumReference> spectrumReferences = new ArrayList<SpectrumReference>();
        List<List<IndexElement>> indexElements = new ArrayList<List<IndexElement>>();

        for (int i = 0; i < filenames.length; i++) {
            List<IndexElement> fileIndex = new ArrayList<IndexElement>();
            List<SpectrumReference> fileSpectrumReferences = parseMgfFile(filenames[i], i, fileIndex);

            spectrumReferences.addAll(fileSpectrumReferences);
            indexElements.add(fileIndex);
        }

        fileIndices = indexElements;

        return spectrumReferences;
    }

    /**
     * Prases a MGF file to extract the SpectrumReferences. Also creates a file index and stores
     * it in the passed list.
     * @param filename
     * @param fileId
     * @param fileIndex Array to hold the file index. This array should be empty.
     * @return
     * @throws Exception
     */
    private List<SpectrumReference> parseMgfFile(String filename, int fileId, List<IndexElement> fileIndex) throws Exception {
        // open the file
        BufferedRandomAccessFile randomAccessFile = new BufferedRandomAccessFile(filename, "r", 1024 * 100);

        List<SpectrumReference> spectrumReferences = new ArrayList<SpectrumReference>();

        // process the file line by line
        String line;
        long currentStart = 0;
        // the end position of the last line = the starting position of the current line
        long lastLineEnd = 0;
        int spectrumIndex = 1; // 1-based index
        float precursorMz = 0;
        boolean inHeader = true;
        boolean hasPeaks = false;


        while ((line = randomAccessFile.readLine()) != null) {
            if (Thread.currentThread().isInterrupted()) {
                randomAccessFile.close();
                throw new InterruptedException();
            }

            // ignore all header fields
            if (line.startsWith("BEGIN IONS")) {
                currentStart = lastLineEnd;
            }

            // save the end position of a spectrum as the current last position
            if (line.startsWith("END IONS")) {
                // save the index element - the index has to be complete
                IndexElement indexElement = new IndexElementImpl(currentStart, (int) (randomAccessFile.getFilePointer() - currentStart));
                fileIndex.add(indexElement);

                // save the spectrum reference only if defined
                if (hasPeaks || !this.ignoreEmptySpectra) {
                    SpectrumReference spectrumReference = new SpectrumReference(fileId, spectrumIndex, precursorMz);
                    spectrumReferences.add(spectrumReference);
                }

                // move to the next spectrum
                spectrumIndex++;

                precursorMz = 0; // to detect any problems
                hasPeaks = false;
            }
            else if (line.startsWith("PEPMASS=")) {
                int index = line.indexOf("=");
                String value = line.substring(index + 1);
                String[] fields = value.split("\\s+");

                precursorMz = Float.parseFloat(fields[0]);
            }
            else if (line.length() > 0 && Character.isDigit(line.charAt(0))) {
                hasPeaks = true;
            }

            lastLineEnd = randomAccessFile.getFilePointer();
        }

        randomAccessFile.close();

        return spectrumReferences;
    }

    public List<List<IndexElement>> getFileIndices() {
        return fileIndices;
    }
}
