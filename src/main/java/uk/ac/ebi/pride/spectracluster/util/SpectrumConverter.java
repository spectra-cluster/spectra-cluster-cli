package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.cluster.GreedySpectralCluster;
import uk.ac.ebi.pride.spectracluster.cluster.ICluster;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ClusteringFileSpectrumReference;
import uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ISpectrumReference;
import uk.ac.ebi.pride.spectracluster.consensus.BinnedGreedyConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.ConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.GreedyConsensusSpectrum;
import uk.ac.ebi.pride.spectracluster.consensus.IConsensusSpectrumBuilder;
import uk.ac.ebi.pride.spectracluster.spectra_list.SpectrumReference;
import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.spectracluster.util.function.spectrum.BinSpectrumMaxFunction;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.CvParam;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.ParamGroup;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.UserParam;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Created by jg on 13.05.15.
 */
public final class SpectrumConverter {
    private SpectrumConverter() {

    }

    /**
     * Convert a cluster read from a .clustering file using the ClusteringFileReader class
     * to a spectra-cluster API cluster object. Currently, this function only creates
     * GreedyClusters.
     *
     * @param readerCluster The cluster read using the ClusteringFileReader class.
     * @return An ICluster object.
     */
    public static ICluster convertClusteringFileReaderCluster(uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster readerCluster) throws Exception {
        // indicates whether the .clustering file contains peak lists
        boolean hasPeakList = false;

        // create the list of spectra
        List<ISpectrum> clusteredSpectra = new ArrayList<>(readerCluster.getSpecCount());

        for (ISpectrumReference specRef : readerCluster.getSpectrumReferences()) {
            List<IPeak> peakList = convertClusteringFileReaderPeaklist(specRef);
            ISpectrum spectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(
                    specRef.getSpectrumId(), specRef.getCharge(), specRef.getPrecursorMz(),
                    Defaults.getDefaultQualityScorer(), peakList);

            if (!hasPeakList && specRef.hasPeaks()) {
                hasPeakList = true;
            }

            clusteredSpectra.add(spectrum);
        }

        // get the consensus spectrum  - force the creation of GreedyConsensusSpectrumBuilder
        IConsensusSpectrumBuilder consensusSpectrum = getClusterFileReaderConsensusSpectrum(readerCluster, false);

        GreedySpectralCluster greedyCluster = new GreedySpectralCluster(
                readerCluster.getId(), clusteredSpectra, (BinnedGreedyConsensusSpectrum) consensusSpectrum,
                Collections.emptyList());

        return greedyCluster;
    }

    /**
     * Create the ConsensusSpectrumBuilder object representing the ClusteringFileReader's cluster's consensus spectrum.
     * @param readerCluster An ICluster object from the ClusteringFileReader.
     * @param hasPeaklist Indicates whether the cluster stores peak lists.
     * @return
     * @throws Exception If the .clustering file does not contain any consensus peak counts
     */
    private static IConsensusSpectrumBuilder getClusterFileReaderConsensusSpectrum(uk.ac.ebi.pride.spectracluster.clusteringfilereader.objects.ICluster readerCluster, boolean hasPeaklist) throws Exception {
        IConsensusSpectrumBuilder consensusSpectrumBuilder = null;

        // create the peaklist
        List<IPeak> peaklist = new ArrayList<>(readerCluster.getConsensusMzValues().size());

        for (int i = 0; i < readerCluster.getConsensusMzValues().size(); i++) {
            int count = 1;
            if (readerCluster.getConsensusCountValues().size() > 0) {
                count = readerCluster.getConsensusCountValues().get(i);
            } else {
                throw new Exception("Missing consensus spectrum peak frequency in .clustering file. Available since spectra-cluster API version 1.0.11");
            }

            IPeak peak = new Peak(readerCluster.getConsensusMzValues().get(i),
                    readerCluster.getConsensusIntensValues().get(i),
                    count);

            peaklist.add(peak);
        }

        // get the sum of the charges
        int sumCharge = 0;

        for (ISpectrumReference specRef : readerCluster.getSpectrumReferences()) {
            sumCharge += specRef.getCharge();
        }

        if (hasPeaklist) {
            consensusSpectrumBuilder = new ConsensusSpectrum(readerCluster.getId(),
                    readerCluster.getSpecCount(),
                    readerCluster.getAvPrecursorMz() * readerCluster.getSpecCount(),
                    readerCluster.getAvPrecursorIntens() * readerCluster.getSpecCount(),
                    sumCharge,
                    peaklist,
                    Defaults.getFragmentIonTolerance());
        } else {
            consensusSpectrumBuilder = new BinnedGreedyConsensusSpectrum(
                    Defaults.getFragmentIonTolerance(),
                    readerCluster.getId(),
                    readerCluster.getSpecCount(),
                    readerCluster.getAvPrecursorMz() * readerCluster.getSpecCount(),
                    readerCluster.getAvPrecursorIntens() * readerCluster.getSpecCount(),
                    sumCharge,
                    peaklist,
                    new BinSpectrumMaxFunction(Defaults.getFragmentIonTolerance()));
        }

        return consensusSpectrumBuilder;
    }

    /**
     * Convert a list of ClusteringFileReader peaks into a list of IPeakS
     * @param specRef The SpectrumReference which peaks should be converted
     * @return
     */
    private static List<IPeak> convertClusteringFileReaderPeaklist(ISpectrumReference specRef) {
        if (!specRef.hasPeaks()) {
            return Collections.emptyList();
        }

        List<IPeak> peaklist = new ArrayList<>(specRef.getPeaks().size());

        for (ClusteringFileSpectrumReference.Peak readerPeak : specRef.getPeaks()) {
            IPeak peak = new Peak(readerPeak.getMz(), readerPeak.getIntensity());
            peaklist.add(peak);
        }

        return peaklist;
    }

    public static ISpectrum convertJmzReaderSpectrum(Spectrum jmzReaderSpectrum, String spectrumId, String peakListFilename) {
        // create the peak list first
        List<IPeak> peaks = new ArrayList<IPeak>();

        for (double mz : jmzReaderSpectrum.getPeakList().keySet()) {
            Peak peak = new Peak((float) mz, (float) jmzReaderSpectrum.getPeakList().get(mz).doubleValue(), 1);
            peaks.add(peak);
        }

        // encode missing charges using "0"
        int charge = 0;
        if (jmzReaderSpectrum.getPrecursorCharge() != null) {
            charge = jmzReaderSpectrum.getPrecursorCharge();
        }

        // create the spectrum
        ISpectrum convertedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(spectrumId,
                charge, (float) jmzReaderSpectrum.getPrecursorMZ().doubleValue(),
                Defaults.getDefaultQualityScorer(), peaks);

        // set the original title if available
        String spectrumTitle = null;

        for (CvParam param : jmzReaderSpectrum.getAdditional().getCvParams()) {
            if ("MS:1000796".equals(param.getAccession())) {
                spectrumTitle = param.getValue();
                convertedSpectrum.setProperty(KnownProperties.SPECTRUM_TITLE,
                        String.format("#file=%s#id=index=%s#title=%s",
                                peakListFilename, jmzReaderSpectrum.getId(), spectrumTitle));
                break;
            }
        }

        // if the title wasn't found, format the default version
        if (spectrumTitle == null) {
            convertedSpectrum.setProperty(KnownProperties.SPECTRUM_TITLE, String.format("" +
                    "#file=%s#id=index=%s#title=Unknown", peakListFilename, jmzReaderSpectrum.getId()));
        }

        // if a spectrum title was found, try to extract the id
        boolean sequenceFound = false;
        if (spectrumTitle != null) {
            int index = spectrumTitle.indexOf("sequence=");
            if (index >= 0) {
                int end = spectrumTitle.indexOf(",", index);
                if (end < 0)
                    end = spectrumTitle.length();

                String sequence = spectrumTitle.substring(index + "sequence=".length(), end);
                convertedSpectrum.setProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY, sequence);
                sequenceFound = true;
            }
            else {
                index = spectrumTitle.indexOf("splib_sequence=");

                if (index >= 0) {
                    int end = spectrumTitle.indexOf(",", index);
                    if (end < 0)
                        end = spectrumTitle.length();

                    String sequence = spectrumTitle.substring(index + "splib_sequence=".length(), end);
                    convertedSpectrum.setProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY, sequence);
                    sequenceFound = true;
                }
            }
        }

        // if there is no peptide sequence encoded in the title, check for "SEQ=" tags
        if (!sequenceFound) {
            ParamGroup paramGroup = jmzReaderSpectrum.getAdditional();

            for (UserParam userParam : paramGroup.getUserParams()) {
                if ("Sequence".equals(userParam.getName())) {
                    convertedSpectrum.setProperty(KnownProperties.IDENTIFIED_PEPTIDE_KEY, userParam.getValue());
                    break;
                }
            }
        }

        // add the retention time if known
        if (Ms2Query.class.isInstance(jmzReaderSpectrum)) {
            Ms2Query query = (Ms2Query) jmzReaderSpectrum;
            if (query.getRetentionTime() != null) {
                convertedSpectrum.setProperty(KnownProperties.RETENTION_TIME, query.getRetentionTime());
            }

            // add the new PRIDE fields
            Map<Integer, String> userTags = query.getUserTags();

            if (userTags != null && userTags.size() > 0) {
                if (userTags.containsKey(3)) {
                    String value = userTags.get(3);
                    if (value.contains("USER03=")) {
                        value = value.substring(7);
                    }
                    convertedSpectrum.setProperty(KnownProperties.MODIFICATION_KEY, value);
                }
                if (userTags.containsKey(5)) {
                    convertedSpectrum.setProperty(KnownProperties.PSM_DECOY_STATUS, userTags.get(5));
                }
                if (userTags.containsKey(6)) {
                    convertedSpectrum.setProperty(KnownProperties.PSM_FDR_SCORES, userTags.get(6));
                }
            }
        }

        return convertedSpectrum;
    }
}
