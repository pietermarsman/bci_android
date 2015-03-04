package edu.nl.ru.analysis;

import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.miscellaneous.ParameterChecker;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;

/**
 * Created by Pieter on 23-2-2015.
 */
public class Classifier {

    private Matrix W, b, spatialFilter, filter, spectrumMx, spectrumKey, windowFn;
    private String[] spectrumDescription;
    private String type, welchWindowType, welchAveType;
    private Boolean verbose, detrend;
    private Double badChannelThreshold, badTrialThreshold;
    private int[] timeIdx, freqIdx, binsp, isBad;
    private Integer dimension, fs, windowLength;
    private Object dvStats, outsz;

    public Classifier(boolean detrend, double badChannelThreshold, double badTrialThreshold, String welchWindowType,
                      String welchAveType, boolean verbose, int[] timeIdx, int[] freqIdx, int dimension, Matrix spatialFilter,
                      Matrix spMx, int windowLength) {
        this.verbose = verbose;
        this.type = "ersp";
        this.detrend = detrend;
        this.badChannelThreshold = badChannelThreshold;
        this.badTrialThreshold = badTrialThreshold;
        this.welchWindowType = welchWindowType;
        ParameterChecker.checkString(welchAveType, new String[]{"amp", "power", "db"});
        this.welchAveType = welchAveType;
        this.timeIdx = timeIdx;
        this.freqIdx = freqIdx;
        this.dimension = dimension;
        this.spatialFilter = spatialFilter;
        this.spectrumMx = spMx;

        // Additional parameters
        this.windowLength = windowLength;
    }

    public ClassifierResult apply(Matrix data) {
        // 0) Bad channel removal
        if (isBad != null) {
            int[] columns = Matrix.range(0, data.getRowDimension(), 1);
            data = new Matrix(data.getSubMatrix(isBad, columns));
        }

        // 1) Detrend
        if (detrend != null && detrend)
            data = data.detrend(2, "linear");

        // 2) Again, bad channel removal
        if (badChannelThreshold != null) {

            // FIXME get the bad channels. Matlab is using t product

            int[] badChannels = new int[]{};
            Matrix car = data.mean(0);
            for (int channel : badChannels) {
                data.setRow(channel, car.getColumn(0));
            }
        }

        // 3) Time range selection
        if (timeIdx != null) {
            int[] rows = Matrix.range(0, data.getRowDimension(), 1);
            data = new Matrix(data.getSubMatrix(rows, timeIdx));
        }

        // 4) Spatial filter
        if (spatialFilter != null) {

            // FIXME do spatial filter. Matlab is using t product.
        }

        // 5) Check for bad trials
        if (badTrialThreshold != null) {

            double x2 = 0.0;
            // FIXME select bad channels. Matlab is using t product.

            boolean isBadTr = x2 > badTrialThreshold;
            // TODO verbose output
        }

        // 6) Convert to spectral
        if (data.getColumnDimension() > windowFn.getRowDimension()) {
            // FIXME other parameters used for welch
//            data = data.welch(windowLength, "density", "linear");
        }

        // 7) Select subset of frequency range (we care about)
        if (freqIdx != null) {
            int[] allRows = Matrix.range(0, data.getRowDimension(), 1);
            data = new Matrix(data.getSubMatrix(allRows, freqIdx));
        }

        // 8) Apply linear classifier
        // FIXME is 1 the right dimension?
        Matrix fraw = linearClassifier(data, 1);
        // FIXME mean of sum over whole matrix is used to get the multi class f??

        // 9) Correct classifier output for bad channels
        // FIXME using a boolean to set the whole matrix to zero?

        // 10) Get probability of the positive class
        fraw.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            public double visit(int row, int column, double value) {
                return 1. / (1. + Math.exp(-value));
            }
        });

        return null;
    }

    private Matrix linearClassifier(Matrix data, int dim) {
        int size = data.getDimension(dim);
        double[][] result = new double[size][this.b.getRowDimension()];
        for (int i = 0; i < size; i++) {
            double[] rowOrColumn;
            if (dim == 0)
                rowOrColumn = data.getRow(i);
            else
                rowOrColumn = data.getColumn(i);
            result[i] = this.W.multiply(new Matrix(rowOrColumn)).add(this.b).getRow(0);
        }
        return new Matrix(result);
    }
}
