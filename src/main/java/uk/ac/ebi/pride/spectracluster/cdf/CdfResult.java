package uk.ac.ebi.pride.spectracluster.cdf;

import java.util.ArrayList;
import java.util.List;

/**
 * This class stores the number of spectra from different peptides
 * for multiple score thresholds and is used to create CDFs
 * from
 *
 * Created by jg on 04.05.15.
 */
public class CdfResult {
    private long totalComparisons = 0;
    /**
     * The "score" windows to use to save the different thresholds
     * for. F.e. for the dot product an increment of 0.5 would mean
     * that values are stored for "0.5 > 1 > 1.5 >..."
     */
    private final double scoreIncrements;
    private List<Long> lowerPeptidesPerScoreIncrement = new ArrayList<Long>();

    public CdfResult(double scoreIncrements) {
        this.scoreIncrements = scoreIncrements;
    }

    /**
     * Saves the similarity score observed for a random (= incorrect)
     * match.
     * @param similarity
     */
    public void saveRandomMatchResult(double similarity) {
        totalComparisons++;

        int bin = getBinForScore(similarity);

        // ignore scores that cannot be stored correctly
        if (bin < 0) {
            return;
        }

        if (lowerPeptidesPerScoreIncrement.size() < bin + 1) {
            int oldSize = lowerPeptidesPerScoreIncrement.size();
            for (int i = oldSize; i < bin + 1; i++) {
                lowerPeptidesPerScoreIncrement.add(0L);
            }
        }

        lowerPeptidesPerScoreIncrement.set(bin, lowerPeptidesPerScoreIncrement.get(bin) + 1L);
    }

    private int getBinForScore(double score) {
        if (score < 0)
            score = 0;

        double doubleBin = score / scoreIncrements;

        int bin = (int) Math.floor(doubleBin);

        return bin;
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder("max_score\tlower_diff_matches\tcum_lower_diff_matches\trel_cum_lower_matches\ttotal_matches\n");

        long totalMatches = 0;

        for (int i = 0; i < lowerPeptidesPerScoreIncrement.size(); i++) {
            totalMatches += lowerPeptidesPerScoreIncrement.get(i);

            stringBuilder
                    .append(String.format("%.2f", scoreIncrements * (i + 1)))
                    .append("\t")
                    .append(lowerPeptidesPerScoreIncrement.get(i))
                    .append("\t")
                    .append(totalMatches)
                    .append("\t")
                    .append((double) totalMatches / totalComparisons)
                    .append("\t")
                    .append(totalComparisons)
                    .append("\n");
        }

        return stringBuilder.toString();
    }

    public void addCdfResult(CdfResult other) throws Exception {
        if (this == other)
            throw new Exception("Cannot join same object");

        if (this.scoreIncrements != other.scoreIncrements)
            throw new Exception("Cannot add cdf result with different score increment (this = " + this.scoreIncrements + ", other = " + other.scoreIncrements + ")");

        int maxSize = Math.max(other.lowerPeptidesPerScoreIncrement.size(), this.lowerPeptidesPerScoreIncrement.size());

        for (int i = 0; i < maxSize; i++) {
            if (i >= this.lowerPeptidesPerScoreIncrement.size()) {
                this.lowerPeptidesPerScoreIncrement.add(other.lowerPeptidesPerScoreIncrement.get(i));
                continue;
            }

            if (i >= other.lowerPeptidesPerScoreIncrement.size()) {
                continue;
            }

            long thisCount = this.lowerPeptidesPerScoreIncrement.get(i);
            long otherCount = other.lowerPeptidesPerScoreIncrement.get(i);

            // override this number
            this.lowerPeptidesPerScoreIncrement.set(i, thisCount + otherCount);
        }

        this.totalComparisons += other.totalComparisons;
    }

    public long getTotalComparisons() {
        return totalComparisons;
    }
}
