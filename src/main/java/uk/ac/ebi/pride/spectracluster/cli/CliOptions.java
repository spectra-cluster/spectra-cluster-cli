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
        FRAGMENT_TOLERANCE("fragment_tolerance"),
        MAJOR_PEAK_JOBS("major_peak_jobs"),
        START_THRESHOLD("threshold_start"),
        END_THRESHOLD("threshold_end"),
        ROUNDS("rounds"),
        BINARY_TMP_DIR("binary_directory"),
        KEEP_BINARY_FILE("keep_binary_files"),
        REUSE_BINARY_FILES("reuse_binary_files"),
        CLUSTER_BINARY_FILE("cluster_binary_file"),
        MERGE_BINARY_RESULTS("merge_binary_results"),
        CONVERT_CGF("convert_cgf"),
        FAST_MODE("fast_mode"),
        HELP("help"),

        // Advanced options
        ADVANCED_MIN_NUMBER_COMPARISONS("x_min_comparisons"),
        ADVANCED_NUMBER_PREFILTERED_PEAKS("x_n_prefiltered_peaks"),
        ADVANCED_LEARN_CDF("x_learn_cdf"),
        ADVANCED_LOAD_CDF_FILE("x_load_cdf");

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
                .withDescription("precursor tolerance (clustering window size) in m/z used during matching.")
                .create(OPTIONS.PRECURSOR_TOLERANCE.getValue());
        options.addOption(precursorTolerance);

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

        Option clusterBinaryFile = OptionBuilder
                .withDescription("if this option is set, only the passed binary file will be clustered and the result written to the file specified in '-output_path' in the binary format")
                .hasArg()
                .create(OPTIONS.CLUSTER_BINARY_FILE.getValue());
        options.addOption(clusterBinaryFile);

        Option mergeBinaryResuls = OptionBuilder
                .withDescription("if this option is set, the passed binary results files are merged into a single .cgf file and written to '-output_path'")
                .create(OPTIONS.MERGE_BINARY_RESULTS.getValue());
        options.addOption(mergeBinaryResuls);

        Option convertCgf = OptionBuilder
                .withDescription("if this option is set the passed CGF file is converted into a .clustering file")
                .create(OPTIONS.CONVERT_CGF.getValue());
        options.addOption(convertCgf);

        Option fastMode = OptionBuilder
                .withDescription("if this option is set the 'fast mode' is enabled. In this mode, the radical peak filtering used for the comparison function is already applied during spectrum conversion. Thereby, the clustering and consensus spectrum quality is slightly decreased but speed increases 2-3 fold.")
                .create(OPTIONS.FAST_MODE.getValue());
        options.addOption(fastMode);

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
    }

    public static Options getOptions() {
        return options;
    }
}
