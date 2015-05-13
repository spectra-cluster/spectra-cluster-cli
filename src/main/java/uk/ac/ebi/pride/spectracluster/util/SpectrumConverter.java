package uk.ac.ebi.pride.spectracluster.util;

import uk.ac.ebi.pride.spectracluster.spectrum.IPeak;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.Peak;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 13.05.15.
 */
public final class SpectrumConverter {
    private SpectrumConverter() {

    }

    public static ISpectrum convertJmzReaderSpectrum(Spectrum jmzReaderSpectrum) {
        // create the peak list first
        List<IPeak> peaks = new ArrayList<IPeak>();

        for (double mz : jmzReaderSpectrum.getPeakList().keySet()) {
            Peak peak = new Peak((float) mz, (float) jmzReaderSpectrum.getPeakList().get(mz).doubleValue(), 1);
            peaks.add(peak);
        }

        // create the spectrum
        ISpectrum convertedSpectrum = new uk.ac.ebi.pride.spectracluster.spectrum.Spectrum(jmzReaderSpectrum.getId(), jmzReaderSpectrum.getPrecursorCharge(), (float) jmzReaderSpectrum.getPrecursorMZ().doubleValue(), Defaults.getDefaultQualityScorer(), peaks);

        return convertedSpectrum;
    }
}
