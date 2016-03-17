package uk.ac.ebi.pride.spectracluster.util;

/**
 * Created by jg on 17.03.16.
 */
public class MissingParameterException extends Exception {
    public MissingParameterException() {
        super();
    }

    public MissingParameterException(String message) {
        super(message);
    }
}
