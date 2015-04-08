package bmird.radboud.fieldtripclientsservice.threads;

import android.content.Context;
import android.util.Log;
import bmird.radboud.fieldtripclientsservice.base.AndroidHandle;
import bmird.radboud.fieldtripclientsservice.base.Argument;
import bmird.radboud.fieldtripclientsservice.base.ThreadBase;
import bmird.radboud.fieldtripclientsservice.threads.analysis.Classifier;
import bmird.radboud.fieldtripclientsservice.threads.analysis.ClassifierResult;
import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.ArrayFunctions;
import edu.nl.ru.miscellaneous.Tuple;
import edu.nl.ru.miscellaneous.Windows;
import nl.fcdonders.fieldtrip.bufferclient.BufferClientClock;
import nl.fcdonders.fieldtrip.bufferclient.BufferEvent;
import nl.fcdonders.fieldtrip.bufferclient.Header;
import nl.fcdonders.fieldtrip.bufferclient.SamplesEventsCount;
import org.apache.commons.math3.linear.RealVector;

import java.io.*;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by Pieter on 23-2-2015.
 * Continuous classifying of data from the buffer and sending events back
 */
public class ContinuousClassifier extends ThreadBase {

    private static final String TAG = ContinuousClassifier.class.toString();

    private String bufferHost;
    private String endValue;
    private String predictionEventType;
    private String baselineEventType;
    private String endType;
    private String baselineEnd;
    private String baselineStart;
    private Integer nBaselineStep;
    private int bufferPort;
    private Double overlap;
    private Double predictionFilter;
    private Integer sampleTrialMs;
    private Integer sampleStepMs;
    private Integer timeoutMs;
    private boolean normalizeLatitude;
    private List<Classifier> classifiers;
    private BufferClientClock C;
    private Integer sampleTrialLength;
    private Integer sampleStep;
    private Float fs;
    private Header header;

    private static List<Classifier> createClassifiers(AndroidHandle android) {
        List<Matrix> Ws = loadWFromFile(3, 56, android.getContext());
        RealVector b = Matrix.zeros(5, 1).getColumnVector(0);
        Integer[] freqIdx = ArrayFunctions.toObjectArray(Matrix.range(0, 56, 1));
        String[] spectrumDescription = new String[]{"alphaL", "alphaR", "baddness", "badChL", "badChR"};
        Integer[] isBad = new Integer[]{0, 0, 0};
        Classifier classifier = new Classifier(Ws, b, true, null, Windows.WindowType.HANNING, WelchOutputType.AMPLITUDE, null, freqIdx, 1, null, null, 128, 100., new Double[]{0.}, spectrumDescription, isBad);
        List<Classifier> classifiers = new LinkedList<Classifier>();
        classifiers.add(classifier);
        return classifiers;
    }

    private static List<Matrix> loadWFromFile(int rows, int columns, Context context) {

        InputStream is = null;
        try {
            is = context.getAssets().open("w.csv");
        } catch (IOException e) {
            Log.e(TAG, Log.getStackTraceString(e));
        }

        List<Matrix> matrices = new LinkedList<Matrix>();
        BufferedReader br = null;
        String line;
        String cvsSplitBy = ",";

        try {
            br = new BufferedReader(new InputStreamReader(is));
            while ((line = br.readLine()) != null) {
                // use comma as separator
                String[] items = line.split(cvsSplitBy);
                Matrix m = new Matrix(ArrayFunctions.fromString(items)).reshape(rows, columns);
                matrices.add(m);
            }
        } catch (FileNotFoundException e) {
            Log.e(TAG, Log.getStackTraceString(e));
        } catch (IOException e) {
            Log.e(TAG, Log.getStackTraceString(e));
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    Log.e(TAG, Log.getStackTraceString(e));
                }
            }
        }
        Log.i(ContinuousClassifier.class.toString(), matrices.get(0).toString());
        return matrices;
    }

    private void setNullFields() {
        // Set trial length
        if (header != null) {
            fs = header.fSample;
        } else {
            Log.e(getClass().toString(), "First connect to the buffer");
        }
        if (sampleTrialLength == null) {
            sampleTrialLength = 0;
            if (sampleTrialMs != null) {
                Float ret = sampleTrialMs / 1000 * fs;
                sampleTrialLength = ret.intValue();
            }
        }

        // Set windows
        for (Classifier c : classifiers) {
            sampleTrialLength = c.getSampleTrialLength(sampleTrialLength);
        }

        // Set wait time
        if (sampleStepMs != null) {
            sampleStep = Double.valueOf(Math.round(sampleStepMs / 1000.0 * fs)).intValue();
        } else {
            sampleStep = Long.valueOf(Math.round(sampleTrialLength * overlap)).intValue();
        }
    }

    private void connect() {
        while (header == null) {
            try {
                Log.i(getClass().toString(), "Connecting to " + bufferHost + ":" + bufferPort);
                C.connect(bufferHost, bufferPort);
                //C.setAutoReconnect(true);
                if (C.isConnected()) {
                    header = C.getHeader();
                }
            } catch (IOException e) {
                header = null;
            }
            if (header == null) {
                Log.w(getClass().toString(), "Invalid Header... waiting");
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    Log.e(TAG, Log.getStackTraceString(e));
                }
            }
        }
    }

    @Override
    public Argument[] getArguments() {
        final Argument[] arguments = new Argument[17];
        arguments[0] = new Argument("Buffer address", "localhost");
        arguments[1] = new Argument("Buffer port", 1972, true);
        arguments[2] = new Argument("Header", null);
        arguments[3] = new Argument("End type", "stimulus.test");
        arguments[4] = new Argument("End value", "end");
        arguments[5] = new Argument("Prediction event type", "classifier.prediction");
        arguments[6] = new Argument("Baseline event type", "stimulus.startbaseline");
        arguments[7] = new Argument("Baseline end", null);
        arguments[8] = new Argument("Baseline start", "start");
        arguments[9] = new Argument("N Baseline step", 5000, true);
        arguments[10] = new Argument("Overlap", .5, true);
        arguments[11] = new Argument("Timeout ms", 1000, true);
        arguments[12] = new Argument("Sample step ms", 500, true);
        arguments[13] = new Argument("Prediction filter", 1.0, true);
        arguments[14] = new Argument("Sample trial length", 25, true);
        arguments[15] = new Argument("Sample trial ms", null);
        arguments[16] = new Argument("Normalize latitude", true);
        return arguments;
    }

    @Override
    public String getName() {
        return "ContinuousClassifier";
    }

    @Override
    public void mainloop() {

        for (Argument a : arguments) {
            Log.i(TAG, "ARGUMENT: " + a.getDescription() + ", " + a.getDefault());
        }

        this.bufferHost = arguments[0].getString();
        this.bufferPort = arguments[1].getInteger();
        this.header = null; //arguments[2];
        this.endType = arguments[3].getString();
        this.endValue = arguments[4].getString();
        this.predictionEventType = arguments[5].getString();
        this.baselineEventType = arguments[6].getString();
        this.baselineEnd = arguments[7].getString();
        this.baselineStart = arguments[8].getString();
        this.nBaselineStep = arguments[9].getInteger();
        this.overlap = arguments[10].getDouble();
        this.timeoutMs = arguments[11].getInteger();
        this.sampleStepMs = arguments[12].getInteger();
        this.predictionFilter = arguments[13].getDouble();
        this.sampleTrialLength = arguments[14].getInteger();
        this.sampleTrialMs = arguments[15].getInteger();
        this.normalizeLatitude = arguments[16].getBoolean();

        this.classifiers = createClassifiers(android);
        this.C = new BufferClientClock();

        this.connect();
        this.setNullFields();

        Log.i(TAG, "Done preparing for mainloop");
        Log.i(TAG, this.toString());

        // Now do the echo-server
        int nEvents = header.nEvents, nSamples = header.nSamples;

        Matrix baseLineVal = Matrix.zeros(classifiers.get(0).getOutputSize() - 1, 1);
        Matrix baseLineVar = Matrix.ones(classifiers.get(0).getOutputSize() - 1, 1);
        boolean baselinephase = false;
        int nBaseline = 0;
        Matrix dvBaseline = null;
        Matrix dv2Baseline = null;
        Matrix dv = null;
        boolean endExpected = false;
        long t0 = 0;

        run = true;

        while (!endExpected && run) {
            // Getting data from buffer
            SamplesEventsCount status = null; // Block until there are new events
            try {
                Log.d(getClass().toString(), "Waiting for " + (nSamples + sampleTrialLength + 1) + " samples");
                status = C.waitForSamples(nSamples + sampleTrialLength + 1, this.timeoutMs);
            } catch (IOException e) {
                Log.e(TAG, Log.getStackTraceString(e));
            }
            if (status.nSamples < header.nSamples) {
                Log.i(getClass().toString(), "Buffer restart detected");
                nSamples = status.nSamples;
                dv = null;
                continue;
            }

            // Logging stuff when nothing is happening
            if (System.currentTimeMillis() - t0 > 5000) {
                Log.i(getClass().toString(), String.format("%5.3f seconds, %d samples, %d events", System.currentTimeMillis() / 1000., status.nSamples, status.nEvents));
                t0 = System.currentTimeMillis();
            }

            // Process any new data
            int onSamples = nSamples;
            int[] start = Matrix.range(onSamples, status.nSamples - sampleTrialLength - 1, sampleStep);
            if (start.length > 0) nSamples = start[start.length - 1] + sampleStep;
//
            for (int from : start) {
                // Get the data
                int to = from + sampleTrialLength - 1;
                Matrix data = null;
                try {
                    data = new Matrix(new Matrix(C.getDoubleData(from, to)).transpose());
                } catch (IOException e) {
                    Log.e(TAG, Log.getStackTraceString(e));
                }
                Log.d(getClass().toString(), String.format("Got data @ %d->%d samples", from, to));

                // Apply classification
                Matrix f = new Matrix(classifiers.get(0).getOutputSize(), 1);
                Matrix fraw = new Matrix(classifiers.get(0).getOutputSize(), 1);
                ClassifierResult result;
                for (Classifier c : classifiers) {
                    result = c.apply(data);
                    Log.d(TAG, result.toString());
                    f = new Matrix(f.add(result.f));
                    fraw = new Matrix(fraw.add(result.fraw));
                }

                // Postprocessing of alpha lat score
                double[] dvColumn = f.getColumn(0);
                f = new Matrix(dvColumn.length - 1, 1);
                f.setColumn(0, Arrays.copyOfRange(dvColumn, 1, dvColumn.length));
                if (normalizeLatitude) f.setEntry(0, 0, (dvColumn[0] - dvColumn[1]) / (dvColumn[0] + dvColumn[1]));
                else f.setEntry(0, 0, dvColumn[0] - dvColumn[1]);

                // Smooth the classifiers
                if (dv == null || predictionFilter == null)
                    dv = f;
                else {
                    if (predictionFilter > 0.) {
                        dv = new Matrix(dv.scalarMultiply(predictionFilter).add(f.scalarMultiply(1. - predictionFilter)));
                    }
                }
                Log.d(TAG, "Result from classifiers: " + Arrays.toString(dv.getColumn(0)));

                // Update baseline
                if (baselinephase) {
                    nBaseline++;
                    dvBaseline = new Matrix(dvBaseline.add(dv));
                    dv2Baseline = new Matrix(dv2Baseline.add(dv.multiplyElements(dv)));
                    if (nBaselineStep != null && nBaseline > nBaselineStep) {
                        Log.i(TAG, "Baseline timeout\n");
                        baselinephase = false;
                        Tuple<Matrix, Matrix> ret = baselineValues(dvBaseline, dv2Baseline, nBaseline);
                        baseLineVal = ret.x;
                        baseLineVar = ret.y;
                    }
                }

                dv = new Matrix(dv.subtract(baseLineVal)).divideElements(baseLineVar);
                Log.d(TAG, "dv: " + dv.toString());

                // Send prediction event
                float[] floatArray = new float[dv.getColumn(0).length];
                Log.d(TAG, "SEND event value: " + Arrays.toString(floatArray));
                for (int i = 0; i < dv.getRowDimension(); i++)
                    try {
                        BufferEvent event = new BufferEvent(predictionEventType, floatArray[i], from);
                        C.putEvent(event);
                    } catch (IOException e) {
                        Log.e(TAG, Log.getStackTraceString(e));
                    }
            }

            // Deal with events
            if (status.nEvents > nEvents) {
                BufferEvent[] events = null;
                try {
                    events = C.getEvents(nEvents, status.nEvents - 1);
                } catch (IOException e) {
                    Log.e(TAG, Log.getStackTraceString(e));
                }

                for (BufferEvent event : events) {
                    String type = event.getType().toString();
                    String value = event.getValue().toString();
                    Log.i(TAG, "GET EVENT (" + event.sample + "): " + type + ", value: " + value);
                    if (type.equals(endType) && value.equals(endValue)) {
                        Log.i(TAG, "End expected");
                        endExpected = true;
                    } else if (type.equals(baselineEventType) && value.equals(baselineEnd)) {
                        Log.i(TAG, "Baseline end event received");
                        baselinephase = false;
                        Tuple<Matrix, Matrix> ret = baselineValues(dvBaseline, dv2Baseline, nBaseline);
                        baseLineVal = ret.x;
                        baseLineVar = ret.y;
                    } else if (type.equals(baselineEventType) && value.equals(baselineStart)) {
                        Log.i(TAG, "Baseline start event received");
                        baselinephase = true;
                        nBaseline = 0;
                        dvBaseline = Matrix.zeros(classifiers.get(0).getOutputSize() - 1, 1);
                        dv2Baseline = Matrix.ones(classifiers.get(0).getOutputSize() - 1, 1);
                    }

                }
                nEvents = status.nEvents;
            }
        }
        try {
            C.disconnect();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Tuple<Matrix, Matrix> baselineValues(Matrix dvBaseline, Matrix dv2Baseline, int nBaseline) {
        double scale = 1. / nBaseline;
        Matrix baseLineVal = new Matrix(dvBaseline.scalarMultiply(scale));
        Matrix baseLineVar = new Matrix(new Matrix(dv2Baseline.subtract(dvBaseline.multiplyElements(dvBaseline).scalarMultiply(scale))).abs().scalarMultiply(scale)).sqrt();
        Log.i(TAG, "Baseline val: " + Arrays.toString(baseLineVal.getColumn(0)));
        Log.i(TAG, "Baseline var: " + Arrays.toString(baseLineVar.getColumn(0)));
        return new Tuple<Matrix, Matrix>(baseLineVal, baseLineVar);
    }

    @Override
    public void validateArguments(Argument[] arguments) {

    }

    public String toString() {
        return "\nContinuousClassifier with parameters:" + "\nBuffer host:  \t" + bufferHost + "\nBuffer port:  \t" +
                bufferPort + "\nHeader:\n     \t" + header + "\nEnd type:     \t" + endType + "\nEnd value:    \t" +
                endValue + "\npredictionEventType:\t" + predictionEventType +
                "\nSampletrialms:\t" + sampleTrialMs + "\nsampleTrialLength:\t" + sampleTrialLength + "\nOverlap:\t" +
                overlap + "\nsampleStepMs:      \t" + sampleStepMs + "\npredictionFilter:\t" + predictionFilter +
                "\nTimeoutMs:\t" + timeoutMs +
                "\nBaselineEnd \t" + baselineEnd +
                "\nBaselineStart\t" + baselineStart +
                "\nBaselineStep\t" + nBaselineStep +
                "\nNormalizeLat\t" + normalizeLatitude +
                "\nFs           \t" + fs;
    }
}
