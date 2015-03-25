package edu.nl.ru.analysis;

import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.Windows;
import nl.fcdonders.fieldtrip.bufferclient.BufferClientClock;
import nl.fcdonders.fieldtrip.bufferclient.BufferEvent;
import nl.fcdonders.fieldtrip.bufferclient.Header;
import nl.fcdonders.fieldtrip.bufferclient.SamplesEventsCount;
import org.apache.commons.math3.linear.RealVector;
import org.apache.log4j.Logger;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by Pieter on 23-2-2015.
 * Continuous classifying of data from the buffer and sending events back
 */
public class ContinuousClassifier implements Runnable {

    private static final Logger log = Logger.getLogger(ContinuousClassifier.class);
    private final String bufferHost;
    private final String endValue;
    private final String predictionEventType;
    private final String endType;
    private final int bufferPort;
    private final Double overlap;
    private final Double predictionFilter;
    private Integer sampleTrialLength;
    private final Integer sampleTrialMs;
    private final Integer sampleStepMs;
    private Integer sampleStep;
    private final Integer timeoutMs;
    private Float fs;
    private Header header;
    private final List<Classifier> classifiers;
    private final BufferClientClock C;

    public ContinuousClassifier(String bufferHost, int bufferPort, Header header, String endType, String endValue,
                                String predictionEventType, Integer sampleTrialLength, Integer sampleTrialMs, Double
            overlap, Integer timeoutMs, List<Classifier> classifiers, Double predictionFilter) {
        log.setLevel(null);
        this.bufferHost = bufferHost;
        this.bufferPort = bufferPort;
        this.header = header;
        this.endType = endType;
        this.endValue = endValue;
        this.predictionEventType = predictionEventType;
        this.sampleTrialLength = sampleTrialLength;
        this.sampleTrialMs = sampleTrialMs;
        this.overlap = overlap;
        this.timeoutMs = timeoutMs;
        this.classifiers = classifiers;
        this.predictionFilter = predictionFilter;

        // Compute parameters
        //        this.sampleStepMs = (int) Math.round(this.sampleTrialLength * this.overlap);
        this.sampleStepMs = null;

        C = new BufferClientClock();
        connect();
        setNullFields();
        log.info(this);
    }

    public static void main(String[] args) {
        Integer[] timeIdx = new Integer[]{0, 1, 2};
        Integer[] freqIdx = new Integer[]{0, 1};
        Matrix W = Matrix.zeros(2, 2);
        List<Matrix> Ws = new LinkedList<Matrix>();
        Ws.add(W);
        RealVector b = Matrix.zeros(2, 1).getColumnVector(0);
        Double[] startMs = new Double[]{0.};
        String[] spectrumDescription = new String[]{"alphaL", "alphaR", "baddness", "badChL", "badChR"};
        Integer[] isBad = new Integer[]{0, 0, 0};
        Classifier classifier = new Classifier(Ws, b, true, null, Windows.WindowType.HANNING, WelchOutputType
                .AMPLITUDE, timeIdx, freqIdx, 1, null, null, 2, 100., startMs, spectrumDescription, isBad);
        List<Classifier> classifiers = new LinkedList<Classifier>();
        classifiers.add(classifier);
        ContinuousClassifier c = new ContinuousClassifier("localhost", 1973, null, "stimulus.test", "end",
                "classifiers.prediction", 25, null, .5, 1000, classifiers, 1.);
        Thread t = new Thread(c);
        t.start();
    }

    private void setNullFields() {
        // Set trial length
        if (sampleTrialLength == null) {
            sampleTrialLength = 0;
            if (sampleTrialMs != null) {
                if (header != null) {
                    fs = header.fSample;
                } else {
                    log.error("First connect to the buffer");
                }
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
            sampleStep = Math.round(sampleStepMs / 1000 * fs);
        } else {
            sampleStep = new Long(Math.round(sampleTrialLength * overlap)).intValue();
        }
    }

    private void connect() {
        while (header == null) {
            try {
                log.info("Connecting to " + bufferHost + ":" + bufferPort);
                C.connect(bufferHost, bufferPort);
                //C.setAutoReconnect(true);
                if (C.isConnected()) {
                    header = C.getHeader();
                }
            } catch (IOException e) {
                header = null;
            }
            if (header == null) {
                log.warn("Invalid Header... waiting");
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    @Override
    public void run() {

        log.info("Running the data processor");
        log.info("#channels....: " + header.nChans);
        log.info("#samples.....: " + header.nSamples);
        log.info("#events......: " + header.nEvents);
        log.info("Sampling Freq: " + header.fSample);
        log.info("data type....: " + header.dataType);

        for (int n = 0; n < header.nChans; n++) {
            if (header.labels[n] != null) {
                log.debug("Ch. " + n + ": " + header.labels[n]);
            }
        }

        // Now do the echo-server
        int nEvents = header.nEvents, nSamples = header.nSamples;

        Matrix dv = null;
        boolean endExpected = false;
        long t0 = 0;

        log.info(this);

        while (!endExpected) {
            // Getting data from buffer
            SamplesEventsCount status = null; // Block until there are new events
            try {
                log.debug("Waiting for " + (nSamples + sampleTrialLength + 1) + " samples");
                status = C.waitForSamples(nSamples + sampleTrialLength + 1, this.timeoutMs);
            } catch (IOException e) {
                e.printStackTrace();
            }
            if (status.nSamples < header.nSamples) {
                log.info("Buffer restart detected");
                nSamples = status.nSamples;
                dv = null;
                continue;
            }

            // Logging stuff when nothing is happening
            if (System.currentTimeMillis() - t0 > 5000) {
                log.info(String.format("%5.3f seconds, %d samples, %d events", System.currentTimeMillis() / 1000.,
                        status.nSamples, status.nEvents));
                t0 = System.currentTimeMillis();
            }

            // Process any new data
            int onSamples = nSamples;
            int[] start = Matrix.range(onSamples, status.nSamples - sampleTrialLength - 1, sampleStep);
            if (start.length > 0)
                nSamples = start[start.length - 1] + sampleStep;

            for (int from : start) {
                // Get the data
                int to = from + sampleTrialLength - 1;
                Matrix data = null;
                try {
                    data = new Matrix(new Matrix(C.getDoubleData(from, to)).transpose());
                } catch (IOException e) {
                    e.printStackTrace();
                }
                log.debug(String.format("Got data @ %d->%d samples", from, to));

                // Apply classification
                Matrix f = new Matrix(classifiers.get(0).getOutputSize(), 1);
                Matrix fraw = new Matrix(classifiers.get(0).getOutputSize(), 1);
                ClassifierResult result;
                for (Classifier c : classifiers) {
                    result = c.apply(data);
                    log.debug(result);
                    f = new Matrix(f.add(result.f));
                    fraw = new Matrix(fraw.add(result.fraw));
                }

                // Smooth the classifiers
                if (dv == null || predictionFilter == null)
                    dv = f;
                else {
                    if (predictionFilter > 0.) {
                        dv = new Matrix(dv.scalarMultiply(predictionFilter).add(f.scalarMultiply(1. -
                                predictionFilter)));
                    }
                }

                // Send prediction event
                BufferEvent event = new BufferEvent(predictionEventType + "_JAVA", dv.getColumn(0), from);
                log.debug("SEND EVENT: " + event.getType() + ", value: " + event.getValue());
                try {
                    C.putEvent(event);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // Deal with events
            if (status.nEvents > nEvents) {
                BufferEvent[] events = null;
                try {
                    events = C.getEvents(nEvents, status.nEvents - 1);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                boolean any = false;
                for (BufferEvent event : events) {
                    any = any || (event.getType().getArray().equals(endType) && event.getValue().getArray().equals
                            (endValue));
                    log.info("GET EVENT (" + event.sample + "): " + event.getType() + ", value: " + event.getValue());
                }
                if (any) {
                    log.info("Got exit event. Stopping");
                    endExpected = true;
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

    public String toString() {
        return "\nContinuousClassifier with parameters:" + "\nBuffer host:  \t" + bufferHost + "\nBuffer port:  \t" +
                bufferPort + "\nHeader:\n     \t" + header + "\nEnd type:     \t" + endType + "\nEnd value:    \t" +
                endValue + "\nLogger level: \t" + log.getLevel() + "\npredEventType:\t" + predictionEventType +
                "\nTrial ms:     \t" + "null" + "\nTrial samples:\t" + sampleTrialLength + "\nOverlap:      \t" +
                overlap + "\nStep ms:      \t" + sampleStepMs + "\nPred filter:    \t" + predictionFilter + "\nTimeout " +
                "ms:   \t" + timeoutMs;
    }
}
