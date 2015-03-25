package edu.nl.ru.analysis;

import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.ArrayFunctions;
import edu.nl.ru.miscellaneous.ParameterChecker;
import edu.nl.ru.miscellaneous.Windows;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealVector;
import org.apache.log4j.Logger;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by Pieter Marsman on 23-2-2015.
 * Classifies a piece of the data using a linear classifier on the welch spectrum.
 */
public class Classifier {

    private static final Logger log = Logger.getLogger(Classifier.class);

    private final double[] windowFn;
    private final Double[] startMs;
    private final Double[] spatialFilter;
    private final List<Matrix> Ws;
    private final Matrix spectrumMx;
    private final RealVector b;
    private final String[] spectrumDescription;
    private final String type;
    private final WelchOutputType welchAveType;
    private final Boolean detrend;
    private final Double badChannelThreshold;
    private final Integer[] timeIdx;
    private final Integer[] freqIdx;
    private final Integer[] isBad;
    private final Integer[] outSize;
    private final Integer dimension;
    private final Integer windowLength;
    private final Double samplingFrequency;
    private final Windows.WindowType windowType;

    public Classifier(List<Matrix> W, RealVector b, Boolean detrend, Double badChannelThreshold, Windows.WindowType
            windowType, WelchOutputType welchAveType, Integer[] timeIdx, Integer[] freqIdx, Integer dimension,
                      Double[] spatialFilter, Matrix spMx, Integer windowLength, Double samplingFrequency, Double[]
            startMs, String[] spectrumDescription, Integer[] isBad) {
        log.setLevel(null);

        ParameterChecker.checkString(welchAveType.toString(), new String[]{"AMPLITUDE", "power", "db"});

        this.type = "ERsP";
        this.detrend = detrend;
        this.welchAveType = welchAveType;
        this.dimension = dimension;
        this.spectrumMx = spMx;
        this.spectrumDescription = spectrumDescription;
        this.isBad = isBad;

        // Classifier
        this.Ws = W;
        this.b = b;

        // Thresholds
        this.badChannelThreshold = badChannelThreshold;

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
        windowFn = Windows.getWindow(windowLength, windowType, true);

        outSize = null;

        log.debug("Finished initializing");
    }

    public ClassifierResult apply(Matrix data) {
        log.debug("Applying classifier on data with shape " + data.shapeString());
        log.debug("With windows size " + windowFn.length);

        if (isBad != null) {
            log.debug("\tBad channel removal");
            int[] columns = Matrix.range(0, data.getColumnDimension(), 1);
            int[] rows = new int[isBad.length];
            int index = 0;
            for (int i = 0; i < isBad.length; i++)
                if (isBad[i] == 0) {
                    rows[index] = i;
                    index++;
                }
            rows = Arrays.copyOf(rows, index);
            data = new Matrix(data.getSubMatrix(rows, columns));
            log.debug("Data shape after bad channel removal: " + data.shapeString());
        }

        if (detrend != null && detrend) {
            data = data.detrend(1, "linear");
            log.debug("Detrending. New data shape: " + data.shapeString());
        }

        List<Integer> badChannels = null;
        if (badChannelThreshold != null) {
            log.debug("\tSecond bad channel removal");
            // TODO use more efficient multiplication of only rows
            Matrix norm = new Matrix(data.multiply(data.transpose()).scalarMultiply(1. / data.getColumnDimension()));
            badChannels = new LinkedList<Integer>();
            for (int r = 0; r < data.getRowDimension(); r++)
                if (norm.getEntry(r, 0) > badChannelThreshold) {
                    log.debug("Channel " + r + " norm is higher than " + badChannelThreshold);
                    badChannels.add(r);
                }

            Matrix car = data.mean(0);
            for (int channel : badChannels) {
                data.setRow(channel, car.getColumn(0));
            }
            log.debug("Second bad channel removal done. New data shape: " + data.shapeString());
        }

        if (timeIdx != null) {
            int[] rows = Matrix.range(0, data.getRowDimension(), 1);
            log.debug("\tTime range selection. Selecting rows: " + Arrays.toString(rows));
            data = new Matrix(data.getSubMatrix(rows, ArrayFunctions.toPrimitiveArray(timeIdx)));
            log.debug("Time range selection done. New data shape: " + data.shapeString());
        }

        if (spatialFilter != null) {
            log.debug("\tSpatial filtering with: " + Arrays.toString(spatialFilter));
            for (int r = 0; r < data.getRowDimension(); r++)
                data.setRowVector(r, data.getRowVector(r).mapMultiply(spatialFilter[r]));
            log.debug("Spatial filtering done. New data shape: " + data.shapeString());
        }

        if (data.getColumnDimension() >= windowFn.length) {
            log.debug("Spectral filtering with welch method");
            log.debug("Input size: " + data.shapeString());
            log.debug("Window size: " + windowFn.length);
            double[] startMs = new double[]{0};
            int[] start = computeSampleStarts(samplingFrequency, startMs);
            data = data.welch(1, windowFn, start, windowLength, true);
            log.debug("Data shape: " + data.shapeString());
        }

        if (freqIdx != null) {
            log.debug("Selecting subset of frequencies: " + Arrays.toString(freqIdx));
            int[] allRows = Matrix.range(0, data.getRowDimension(), 1);
            data = new Matrix(data.getSubMatrix(allRows, ArrayFunctions.toPrimitiveArray(freqIdx)));
            log.debug("Data shape: " + data.shapeString());
        }

        // FIXME is 0 the right dimension?
        Matrix fraw = linearClassifier(data, 0);
        log.debug("Classifying with linear classifier");
        // FIXME mean of sum over whole matrix is used to get the multi class f??
        log.debug("fraw shape: " + fraw.shapeString());

        Matrix f = new Matrix(fraw.copy());
        if (badChannels != null) {
            log.debug("Removing bad channels from classification: " + Arrays.toString(badChannels.toArray()));
            for (int channel : badChannels)
                f.setRowMatrix(channel, Matrix.zeros(1, f.getColumnDimension()));
            log.debug("f shape: " + f.shapeString());
        }

        Matrix p = new Matrix(f.copy());
        p.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            public double visit(int row, int column, double value) {
                return 1. / (1. + Math.exp(-value));
            }
        });
        log.debug("Computed class probability: " + p);
        log.debug("p shape: " + p.shapeString());
        log.debug("Done classifying");
        return new ClassifierResult(f, fraw, p, data);
    }

    private Matrix linearClassifier(Matrix data, int dim) {
        log.debug("Applying linear classifier (" + Ws.get(0).shapeString() + "->" + getOutputSize() + ") on data with" +
                " shape " + data.shapeString() + " on dimension " + dim);
        double[] results = new double[Ws.size()];
        for (int i = 0; i < Ws.size(); i++)
            results[i] = this.Ws.get(i).multiplyElements(data).sum().getEntry(0, 0) + b.getEntry(i);
        return new Matrix(results);
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

    public Integer getSampleTrialLength(Integer sampleTrialLength) {
        if (outSize != null)
            return Math.max(sampleTrialLength, outSize[0]);
        else if (timeIdx != null)
            return Math.max(sampleTrialLength, timeIdx[1]);
        else if (windowFn != null)
            return Math.max(sampleTrialLength, windowFn.length);
        throw new RuntimeException("Either outSize, timeIdx or windowFn should be defined");
    }

    public int getOutputSize() {
        return Ws.size();
    }

    public String toString() {
        return "Classifier with parameters:" + "\nWindow Fn length:  \t" + windowFn.length + "\nstartMs            " +
                "\t" + Arrays.toString(startMs) + "\nSpatial filter     \t" + Arrays.toString(spatialFilter) + "\nWs shape            " +
                "\t" + (Ws != null ? Ws.get(0).shapeString() : "null") + "\nSpectrum mx shape  \t" + (spectrumMx !=
                null ? spectrumMx.shapeString() : "null") + "\nb                  \t" + b + "\nSpectrum desc      \t"
                + Arrays.toString(spectrumDescription) + "\nType               \t" + type + "\nWelch ave type     \t"
                + welchAveType + "\nDetrend            \t" + detrend + "\nBad channel thres  \t" +
                badChannelThreshold + "\nTime idx           \t" + Arrays.toString(timeIdx) + "\nFrequency idx      \t" + Arrays
                .toString(freqIdx) + "\nIs bad channel    \t" + Arrays.toString(isBad) + "\nDimension          \t" +
                dimension + "\nWindow length      \t" + windowLength + "\nSampling frequency \t" + samplingFrequency
                + "\nWindow type        \t" + windowType;
    }
}
