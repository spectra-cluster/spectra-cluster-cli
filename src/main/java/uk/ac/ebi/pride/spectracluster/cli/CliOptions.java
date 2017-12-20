package uk.ac.ebi.pride.spectracluster.cli;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

/**
 * Created with IntelliJ IDEA.
 * User: jg
 * Date: 9/15/13
 * Time: 11:35 AM
 * To change this template use File | Settings | File Templates.
 */
public class CliOptions {
    public enum OPTIONS {
        OUTPUT_PATH("output_path"),
        PRECURSOR_TOLERANCE("precursor_tolerance"),
        PRECURSOR_TOELRANCE_UNIT("precursor_tolerance_unit"),
        FRAGMENT_TOLERANCE("fragment_tolerance"),
        MAJOR_PEAK_JOBS("major_peak_jobs"),
        START_THRESHOLD("threshold_start"),
        END_THRESHOLD("threshold_end"),
        ROUNDS("rounds"),
        BINARY_TMP_DIR("binary_directory"),
        KEEP_BINARY_FILE("keep_binary_files"),
        REUSE_BINARY_FILES("reuse_binary_files"),
        REMOVE_REPORTER_PEAKS("remove_reporters"),
        FAST_MODE("fast_mode"),
        ADD_SCORES("add_scores"),
        FILTER("filter"),
        HELP("help"),
        VERBOSE("verbose"),

        // Advanced options
        ADVANCED_MIN_NUMBER_COMPARISONS("x_min_comparisons"),
        ADVANCED_MIN_ADAPTIVE_COMPARISONS("x_min_adapt_comparisons"),
        ADVANCED_NUMBER_PREFILTERED_PEAKS("x_n_prefiltered_peaks"),
        ADVANCED_LEARN_CDF("x_learn_cdf"),
        ADVANCED_LOAD_CDF_FILE("x_load_cdf"),
        ADVANCED_DISABLE_MGF_COMMENTS("x_disable_mgf_comments"),
        ADVANCED_MIN_CONSENSUS_PEAKS_TO_KEEP("x_min_consensus_peaks_to_keep"),
        ADVANCED_MERGE_BINARY_FILES("x_merge_binary_files"),
        ADVANCED_CONVERT_CGF("x_convert_cgf"),
        ADVANCED_FILTER_PEAKS_PER_MZ("x_filter_peaks_mz");

        private String value;

        OPTIONS(String value) {
            this.value = value;
        }

        public String getValue() {
            return value;
        }

        @Override
        public String toString() {
            return value;
        }
    }

    private static final Options options = new Options();

    static {
        Option fragmentTolerance = OptionBuilder
                .hasArg()
                .withDescription("fragment ion tolerance in m/z to use for fragment peak matching")
                .create(OPTIONS.FRAGMENT_TOLERANCE.getValue());
        options.addOption(fragmentTolerance);

        Option precursorTolerance = OptionBuilder
                .hasArg()
                .withDescription("precursor tolerance (clustering window size) in m/z (default) or ppm used during matching." +
                        " The unit can be changed by setting the -" + OPTIONS.PRECURSOR_TOELRANCE_UNIT.getValue() + " parameter")
                .create(OPTIONS.PRECURSOR_TOLERANCE.getValue());
        options.addOption(precursorTolerance);

        Option precursorToleranceUnit = OptionBuilder
                .hasArg()
                .withDescription("sets the precursor tolerance unit. Allowed values are \"mz\" and \"ppm\". " +
                                "Default is \"mz\".")
                .create(OPTIONS.PRECURSOR_TOELRANCE_UNIT.getValue());
        options.addOption(precursorToleranceUnit);

        Option outputPath = OptionBuilder
                .hasArg()
                .withDescription("path to the outputfile. Outputfile must not exist.")
                .create(OPTIONS.OUTPUT_PATH.getValue());
        options.addOption(outputPath);

        Option startThreshold = OptionBuilder
                .hasArg()
                .withDescription("(highest) starting threshold")
                .create(OPTIONS.START_THRESHOLD.getValue());
        options.addOption(startThreshold);

        Option endThreshold = OptionBuilder
                .hasArg()
                .withDescription("(lowest) final clustering threshold")
                .create(OPTIONS.END_THRESHOLD.getValue());
        options.addOption(endThreshold);

        Option rounds = OptionBuilder
                .hasArg()
                .withDescription("number of clustering rounds to use.")
                .create(OPTIONS.ROUNDS.getValue());
        options.addOption(rounds);

        Option majorPeakJobs = OptionBuilder
                .hasArg()
                .withDescription("number of threads to use for major peak clustering.")
                .create(OPTIONS.MAJOR_PEAK_JOBS.getValue());
        options.addOption(majorPeakJobs);

        Option binaryDirectory = OptionBuilder
                .hasArg()
                .withDescription("path to the directory to (temporarily) store the binary files. By default a temporary directory is being created")
                .create(OPTIONS.BINARY_TMP_DIR.getValue());
        options.addOption(binaryDirectory);

        Option keepBinary = OptionBuilder
                .withDescription("if this options is set, the binary files are not deleted after clustering.")
                .create(OPTIONS.KEEP_BINARY_FILE.getValue());
        options.addOption(keepBinary);

        Option reuseBinaryFiles = OptionBuilder
                .withDescription("if this option is set, the binary files found in the binary file directory will be used for clustering.")
                .create(OPTIONS.REUSE_BINARY_FILES.getValue());
        options.addOption(reuseBinaryFiles);

        Option fastMode = OptionBuilder
                .withDescription("if this option is set the 'fast mode' is enabled. In this mode, the radical peak filtering used for the comparison function is already applied during spectrum conversion. Thereby, the clustering and consensus spectrum quality is slightly decreased but speed increases 2-3 fold.")
                .create(OPTIONS.FAST_MODE.getValue());
        options.addOption(fastMode);

        Option filter = OptionBuilder
                .hasArg()
                .withDescription("adds a filter to be applied to the input spectrum. Available values are ['immonium_ions', 'mz_150', 'mz_200']")
                .create(OPTIONS.FILTER.getValue());
        options.addOption(filter);

        Option verbose = OptionBuilder
                .withDescription("if set additional status information is printed.")
                .create(OPTIONS.VERBOSE.getValue());
        options.addOption(verbose);

        Option removeReporters = OptionBuilder
                .hasArg()
                .withArgName("QUANTITATION TYPE")
                .withDescription("remove reporter ion peaks in quantitation experiments. Possible QUANTIATION TYPES are 'ITRAQ', 'TMT' and 'ALL' ('TMT' and 'ITRAQ' peaks are removed.")
                .create(OPTIONS.REMOVE_REPORTER_PEAKS.getValue());
        options.addOption(removeReporters);

        Option addScores = OptionBuilder
                .withDescription("if set, the similarity scores of each spectrum to the cluster's consensus spectrum is" +
                        "added to the output file.")
                .create(OPTIONS.ADD_SCORES.getValue());
        options.addOption(addScores);

        Option help = new Option(
                OPTIONS.HELP.toString(),
                "print this message.");
        options.addOption(help);

        /**
         * ADVANCED OPTIONS
         */
        Option xMinComparisons = OptionBuilder
                .withDescription("(Experimental option) Sets the minimum number of comparisons used to calculate the probability that incorrect spectra are clustered.")
                .hasArg()
                .create(OPTIONS.ADVANCED_MIN_NUMBER_COMPARISONS.getValue());
        options.addOption(xMinComparisons);

        Option xMinAdaptComparisons = OptionBuilder
                .withDescription("(Experimental option) Uses an adaptive function for the minimum comparisons together with a minimum number always to return.")
                .hasArg()
                .create(OPTIONS.ADVANCED_MIN_ADAPTIVE_COMPARISONS.getValue());
        options.addOption(xMinAdaptComparisons);

        Option xLearnCdf = OptionBuilder
                .hasArg()
                .withArgName("output filename")
                .withDescription("(Experimental option) Learn the used cumulative distribution function directly from the processed data. This is only recommended for high-resolution data. The result will be written to the defined file.")
                .create(OPTIONS.ADVANCED_LEARN_CDF.getValue());
        options.addOption(xLearnCdf);

        Option xLoadCdf = OptionBuilder
                .hasArg()
                .withArgName("CDF filename")
                .withDescription("(Experimental option) Loads the cumulative distribution function to use from the specified file. These files can be created using the " + OPTIONS.ADVANCED_LEARN_CDF.getValue() + " parameter")
                .create(OPTIONS.ADVANCED_LOAD_CDF_FILE.getValue());
        options.addOption(xLoadCdf);

        Option xNumberPrefilteredPeaks = OptionBuilder
                .hasArg()
                .withArgName("number peaks")
                .withDescription("(Experimental option) Set the number of highest peaks that are kept per spectrum during loading.")
                .create(OPTIONS.ADVANCED_NUMBER_PREFILTERED_PEAKS.getValue());
        options.addOption(xNumberPrefilteredPeaks);

        Option xDisableMgfComments = OptionBuilder
                .withDescription("(Advanced option) If set, MGF comment strings are NOT supported. This will increase performance but only works for MGF files that do not contain any comments")
                .create(OPTIONS.ADVANCED_DISABLE_MGF_COMMENTS.getValue());
        options.addOption(xDisableMgfComments);

        Option xMinConsensusPeaksToKeep = OptionBuilder
                .withDescription("(Advanced option) Sets the minimum number of peaks to keep in a consensus spectrum. If" +
                        " the consensus spectrum contains fewer than these peaks, the noise filtering which retains N" +
                        " peaks per M m/z is not used but all peaks are kept")
                .hasArg()
                .withArgName("number peaks")
                .create(OPTIONS.ADVANCED_MIN_CONSENSUS_PEAKS_TO_KEEP.getValue());
        options.addOption(xMinConsensusPeaksToKeep);

        Option xMergeBinaryfiles = OptionBuilder
                .withDescription("(Advanced option) If this option is set, the input files should be binary" +
                                " result files from the clustering process. In this mode, the spectra-cluster-cli " +
                                "only merges the binary files and creates the final output file.")
                .create(OPTIONS.ADVANCED_MERGE_BINARY_FILES.getValue());

        Option xConvertCgf = OptionBuilder
                .withDescription("(Advanced option) It this option is set, the input file in CGF format will be " +
                        "converted to the .clustering file. No other processing will be performed.")
                .create(OPTIONS.ADVANCED_CONVERT_CGF.getValue());
        options.addOption(xConvertCgf);

        Option xFilterPeaksMz = OptionBuilder
                .withDescription("(Advanced option) If this option is set, the top 10 peaks per 100 m/z are retained instead of the top N peaks overall.")
                .create(OPTIONS.ADVANCED_FILTER_PEAKS_PER_MZ.getValue());
        options.addOption(xFilterPeaksMz);

        options.addOption(xMergeBinaryfiles);
    }

    public static Options getOptions() {
        return options;
    }
}
