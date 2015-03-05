package edu.nl.ru.analysis;

import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.ParameterChecker;
import edu.nl.ru.miscellaneous.Windows;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealVector;

import java.util.LinkedList;
import java.util.List;

/**
 * Created by Pieter on 23-2-2015.
 */
public class Classifier {

    // Spatial (channel) filter
    double[] spatialFilter;
    // Welch param
    double[] startMs;
    private Matrix W, filter, spectrumMx, spectrumKey;
    private RealVector b;
    private String[] spectrumDescription;
    private String type;
    private WelchOutputType welchAveType;
    private Boolean verbose, detrend;
    private Double badChannelThreshold, badTrialThreshold;
    private int[] timeIdx, freqIdx, binsp, isBad;
    private Integer dimension, fs, windowLength;
    private Object dvStats, outsz;
    // Data acquisition properties
    private double samplingFrequency;
    // Window
    private Windows.WindowType windowType;
    private double[] windowFn;

    public Classifier(Matrix W, RealVector b, boolean detrend, double badChannelThreshold, double badTrialThreshold,
                      Windows.WindowType windowType, WelchOutputType welchAveType, boolean verbose, int[] timeIdx,
                      int[] freqIdx, int dimension, double[] spatialFilter, Matrix spMx, int windowLength, double
                              samplingFrequency, double[] startMs) {

        this.verbose = verbose;
        this.type = "ersp";
        this.detrend = detrend;
        ParameterChecker.checkString(welchAveType.toString(), new String[]{"AMPLITUDE", "power", "db"});
        this.welchAveType = welchAveType;
        this.dimension = dimension;
        this.spectrumMx = spMx;

        // Classifier
        this.W = W;
        this.b = b;

        // Thresholds
        this.badChannelThreshold = badChannelThreshold;
        this.badTrialThreshold = badTrialThreshold;

        // Data acquisition properties
        this.samplingFrequency = samplingFrequency;

        // Additional parameters
        this.windowLength = windowLength;

        // Spatial filter
        this.spatialFilter = spatialFilter;

        // Welch param
        this.startMs = startMs;

        // Create window
        this.timeIdx = timeIdx;
        this.freqIdx = freqIdx;
        this.windowType = windowType;
        windowFn = Windows.getWindow(windowLength, windowType);
    }

    public ClassifierResult apply(Matrix data) {
        // 0) Bad channel removal
        if (isBad != null) {
            int[] columns = Matrix.range(0, data.getRowDimension(), 1);
            data = new Matrix(data.getSubMatrix(isBad, columns));
        }

        // 1) Detrend
        if (detrend != null && detrend)
            data = data.detrend(1, "linear");

        // 2) Again, bad channel removal
        List<Integer> badChannels = null;
        if (badChannelThreshold != null) {
            // TODO use more efficient multiplication of only rows
            Matrix norm = new Matrix(data.multiply(data.transpose()).scalarMultiply(1. / data.getColumnDimension()));
            badChannels = new LinkedList<Integer>();
            for (int r = 0; r < data.getRowDimension(); r++)
                badChannels.add(r);

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
            for (int r = 0; r < data.getRowDimension(); r++)
                data.setRowVector(r, data.getRowVector(r).mapMultiply(spatialFilter[r]));
        }

        // 5) Bad trial selection
        // Not needed because only considering a single epoch

        // 6) Convert to spectral
        if (data.getColumnDimension() > windowFn.length) {
            double[] startMs = new double[]{0};
            int[] start = computeSampleStarts(samplingFrequency, startMs);
            data = data.welch(dimension, windowFn, start, windowLength, true);
        }

        // 7) Select subset of frequency range (we care about)
        if (freqIdx != null) {
            int[] allRows = Matrix.range(0, data.getRowDimension(), 1);
            data = new Matrix(data.getSubMatrix(allRows, freqIdx));
        }

        // 8) Apply linear classifier
        // FIXME is 0 the right dimension?
        Matrix fraw = linearClassifier(data, 0);
        // FIXME mean of sum over whole matrix is used to get the multi class f??

        // 9) Correct classifier output for bad channels
        Matrix f = new Matrix(fraw.copy());
        if (badChannels != null) {
            for (int channel : badChannels)
                f.setRowMatrix(channel, Matrix.zeros(1, f.getColumnDimension()));
        }

        // 10) Get probability of the positive class
        Matrix p = new Matrix(f.copy());
        p.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            public double visit(int row, int column, double value) {
                return 1. / (1. + Math.exp(-value));
            }
        });

        return new ClassifierResult(f, fraw, p, data);
    }

    private Matrix linearClassifier(Matrix data, int dim) {
        int size = data.getDimension(dim);
        double[][] result = new double[size][this.b.getDimension()];
        for (int i = 0; i < size; i++) {
            double[] rowOrColumn;
            if (dim == 0)
                rowOrColumn = data.getRow(i);
            else
                rowOrColumn = data.getColumn(i);
            result[i] = this.W.operate(new ArrayRealVector(rowOrColumn)).add(this.b).toArray();
        }
        return new Matrix(result);
    }

    public int computeSampleWidth(double samplingFrequency, double widthMs) {
        return (int) Math.floor(widthMs * (samplingFrequency / 1000.));
    }

    public int[] computeSampleStarts(double samplingFrequency, double[] startMs) {
        int[] sampleStarts = new int[startMs.length];
        for (int i = 0; i < startMs.length; i++)
            sampleStarts[i] = (int) Math.floor(startMs[i] * (samplingFrequency / 1000.));
        return sampleStarts;
    }
}
