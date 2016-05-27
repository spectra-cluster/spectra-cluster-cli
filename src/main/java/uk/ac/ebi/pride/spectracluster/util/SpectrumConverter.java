package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.tools.jmzreader.model.Param;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.CvParam;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.ParamGroup;
import uk.ac.ebi.pride.tools.jmzreader.model.impl.UserParam;

import java.util.ArrayList;
import java.util.List;
import java.util.UUID;

/**
 * Created by jg on 13.05.15.
 */
public final class SpectrumConverter {
    private SpectrumConverter() {

    }

    public static ISpectrum convertJmzReaderSpectrum(Spectrum jmzReaderSpectrum, String spectrumId, String peakListFilename) {
        // create the peak list first
        List<IPeak> peaks = new ArrayList<IPeak>();

        for (double mz : jmzReaderSpectrum.getPeakList().keySet()) {
            Peak peak = new Peak((float) mz, (float) jmzReaderSpectrum.getPeakList().get(mz).doubleValue(), 1);
            peaks.add(peak);
        }

        String.format("#file=%s#id=%s#title=%s");

        // create the spectrum
        ISpectrum convertedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(spectrumId,
                jmzReaderSpectrum.getPrecursorCharge(), (float) jmzReaderSpectrum.getPrecursorMZ().doubleValue(),
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

        return convertedSpectrum;
    }
}
