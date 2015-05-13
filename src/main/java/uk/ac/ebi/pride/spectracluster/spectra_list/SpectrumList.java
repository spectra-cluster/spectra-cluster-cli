package uk.ac.ebi.pride.spectracluster.spectra_list;

import uk.ac.ebi.pride.spectracluster.spectrum.Spectrum;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import java.util.List;

/**
 * Created by jg on 13.05.15.
 */
public class SpectrumList implements Externalizable{
    private Spectrum[] spectra;

    public SpectrumList(List<Spectrum> references) {
        spectra = references.toArray(new Spectrum[references.size()]);
    }

    @Override
    public void writeExternal(ObjectOutput out) throws IOException {
        out.writeInt(spectra.length);

        for (Spectrum spectrum : spectra) {
            out.writeObject(spectrum);
        }
    }

    @Override
    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
        int nSpectra = in.readInt();
        spectra = new Spectrum[nSpectra];

        for (int i = 0; i < nSpectra; i++) {
            Spectrum spectrum = (Spectrum) in.readObject();
            spectra[i] = spectrum;
        }
    }

    public Spectrum[] getSpectra() {
        return spectra;
    }
}
