package uk.ac.ebi.pride.spectracluster.util;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import uk.ac.ebi.pride.spectracluster.spectrum.ISpectrum;
import uk.ac.ebi.pride.spectracluster.spectrum.KnownProperties;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.io.File;

/**
 * Created by jg on 21.06.17.
 */
public class TestPrideMgfFields {
    private File testfile;

    @Before
    public void setUp() throws Exception {
        testfile = new File(TestPrideMgfFields.class.getClassLoader().getResource("pride_export_2017.mgf").toURI());
    }

    @Test
    public void testPrideFields() throws Exception {
        MgfFile reader = new MgfFile(testfile);
        reader.setDisableCommentSupport(true);

        Ms2Query testQuery = reader.getMs2Query(1);

        ISpectrum convertedSpectrum = SpectrumConverter.convertJmzReaderSpectrum(testQuery,
                "no_id", "testfile.mgf");

        Assert.assertNotNull(convertedSpectrum.getProperty(KnownProperties.MODIFICATION_KEY));
        Assert.assertNotNull(convertedSpectrum.getProperty(KnownProperties.PSM_FDR_SCORES));
        Assert.assertNotNull(convertedSpectrum.getProperty(KnownProperties.PSM_DECOY_STATUS));
    }
}
