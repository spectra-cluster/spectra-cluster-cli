package uk.ac.ebi.pride.spectracluster.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;

/**
 * Created with IntelliJ IDEA.
 * User: jg
 * Date: 9/15/13
 * Time: 11:19 AM
 * To change this template use File | Settings | File Templates.
 */
public class SpectraClusterCliMain {
    public static void main(String[] args) {
        CommandLineParser parser = new GnuParser();

        try {
            CommandLine commandLine = parser.parse(CliOptions.getOptions(), args);

            // 1.) Pre-process spectra and store in list per highest peak
            /**
             * pre-scan all files and create list of SpectrumReferences for all spectra
             * > already store each spectrum N times for each highest peak
             *   > NOTE: this already requires the precursor filter to be applied
             * sort each list according to precursor m/z
             * write spectra to files in this specific order (use ObjectOutputStream)
             */
            String[] filenames = commandLine.getArgs();

            // 2.) Cluster each list in one thread
            /**
             * For each binary file one thread is started
             * Cluster spectra using decreasing similarity thresholds
             * ?write to CGF file or return complete list > may cause memory issues
             */

            // 3.) merge duplicated clusters
            /**
             * pre-scan all generated clusters
             *  > store in one big list, sort according to m/z
             * run incremental clustering directly on this list / big file
             */


            // HELP
            if (commandLine.hasOption(CliOptions.OPTIONS.HELP.getValue())) {
                printUsage();
                return;
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Error: " + e.getMessage());

            System.exit(1);
        }
    }

   private static void printUsage() {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Spectra Cluster - Clusterer",
                "Clusters the spectra found in an MGF file and writes the results in a text-based file.\n",
                CliOptions.getOptions(), "\n\n", true);
    }


}
