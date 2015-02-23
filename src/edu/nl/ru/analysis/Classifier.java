package edu.nl.ru.analysis;

import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.miscellaneous.ParameterChecker;

/**
 * Created by Pieter on 23-2-2015.
 */
public class Classifier {

    boolean verbose;
    String type, welchWindowType, aveType;
    boolean detrend;
    double badChannelThreshold, badTrThreshold;
    int[] timeIdx, freqIdx;
    int dimension;
    Matrix W, b, spatialFilter, spMx;

    public Classifier(boolean detrend, double badChannelThreshold, double badTrThreshold, String welchWindowType,
                      String aveType, boolean verbose, int[] timeIdx, int[] freqIdx, int dimension, Matrix spatialFilter,
                      Matrix spMx) {
        this.verbose = verbose;
        this.type = "ersp";
        this.detrend = detrend;
        this.badChannelThreshold = badChannelThreshold;
        this.badTrThreshold = badTrThreshold;
        this.welchWindowType = welchWindowType;
        ParameterChecker.checkString(aveType, new String[]{"amp", "power", "db"});
        this.aveType = aveType;
        this.timeIdx = timeIdx;
        this.freqIdx = freqIdx;
        this.dimension = dimension;
        this.spatialFilter = spatialFilter;
        this.spMx = spMx;
    }

    public ClassifierResult apply(Matrix data) {
        // 0) Bad channel removal

        // 1) Detrend

        // 2) Again, bad channel removal

        // 3) Time range selection

        // 4) Spatial filter

        // 5) Bad trials

        // 6) Convert to spectral

        // 7) Select subset of frequency range (we care about)

        // 8) Apply linear classifier

        // 9) Correct classifier output for bad channels

        // 10) Get probability of the positive class

        return null;
    }
}
