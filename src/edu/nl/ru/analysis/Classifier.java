package edu.nl.ru.analysis;

import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.ArrayFunctions;
import edu.nl.ru.miscellaneous.ParameterChecker;
import edu.nl.ru.miscellaneous.Windows;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealVector;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by Pieter Marsman on 23-2-2015.
 * Classifies a piece of the data using a linear classifier on the welch spectrum.
 */
public class Classifier {

    private static Logger log = Logger.getLogger(Classifier.class);

    private double[] windowFn;
    private Double[] startMs, spatialFilter;
    private List<Matrix> Ws;
    private Matrix filter, spectrumMx, spectrumKey;
    private RealVector b;
    private String[] spectrumDescription;
    private String type;
    private WelchOutputType welchAveType;
    private Boolean detrend;
    private Double badChannelThreshold, badTrialThreshold;
    private Integer[] timeIdx, freqIdx, binsp, isBad, outSize;
    private Integer dimension, fs, windowLength;
    private Object dvStats;
    private Double samplingFrequency;
    private Windows.WindowType windowType;

    public Classifier(List<Matrix> W, RealVector b, Boolean detrend, Double badChannelThreshold, Double badTrialThreshold,
                      Windows.WindowType windowType, WelchOutputType welchAveType, Integer[] timeIdx, Integer[]
            freqIdx, Integer dimension, Double[] spatialFilter, Matrix spMx, Integer windowLength, Double
            samplingFrequency, Double[] startMs, String[] spectrumDescription, Integer[] isBad) {
        log.setLevel(Level.DEBUG);

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
        windowFn = Windows.getWindow(windowLength, windowType, true);

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

        log.info(this);
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
        log.debug(f);
        log.debug(fraw);
        log.debug(p);
        return new ClassifierResult(f, fraw, p, data);
    }

    private Matrix linearClassifier(Matrix data, int dim) {
        log.debug("Applying linear classifier (" + Ws.get(0).shapeString() + "->" + getOutputSize() + ") on data with shape " + data.shapeString() + " on dimension " + dim);
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
        StringBuilder sb = new StringBuilder();
        sb.append("Classifier with parameters:");
        sb.append("\nWindow Fn length:  \t").append(windowFn.length);
        sb.append("\nstartMs            \t").append(Arrays.toString(startMs));
        sb.append("\nSpatial filter     \t").append(spatialFilter);
        sb.append("\nWs shape            \t").append(Ws != null ? Ws.get(0).shapeString() : "null");
        sb.append("\nFilter shape       \t").append(filter != null ? filter.shapeString() : "null");
        sb.append("\nSpectrum mx shape  \t").append(spectrumMx != null ? spectrumMx.shapeString() : "null");
        sb.append("\nSpectrum key       \t").append(spectrumKey != null ? spectrumKey.shapeString() : "null");
        sb.append("\nb                  \t").append(b);
        sb.append("\nSpectrum desc      \t").append(Arrays.toString(spectrumDescription));
        sb.append("\nType               \t").append(type);
        sb.append("\nWelch ave type     \t").append(welchAveType);
        sb.append("\nDetrend            \t").append(detrend);
        sb.append("\nBad channel thres  \t").append(badChannelThreshold);
        sb.append("\nTime idx           \t").append(timeIdx);
        sb.append("\nFrequency idx      \t").append(Arrays.toString(freqIdx));
        sb.append("\nBinsp              \t").append(binsp);
        sb.append("\nIs bad channeld    \t").append(Arrays.toString(isBad));
        sb.append("\nOut size           \t").append(outSize);
        sb.append("\nDimension          \t").append(dimension);
        sb.append("\nFs                 \t").append(fs);
        sb.append("\nWindow length      \t").append(windowLength);
        sb.append("\ndvStats            \t").append(dvStats);
        sb.append("\nSampling frequency \t").append(samplingFrequency);
        sb.append("\nWindow type        \t").append(windowType);
        return sb.toString();
    }
}
