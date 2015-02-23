package edu.nl.ru.analysis;

/**
 * Created by Pieter on 23-2-2015.
 */
public class DataProcessor implements Runnable {

    String buffHost, endValue;
    NFEvent endType, predEventType, baselineEventType;
    int buffPort, baselineEnd, timeout_ms;
    double overlap;
    double[] lengthTrMs, lengthTrSample, stepMs, predFilt;
    boolean normLat, verbose;
    Object[] hdr;


    public DataProcessor() {
        buffHost = "localhost";
        buffPort = 1972;
        hdr = new Object[]{};
        endType = NFEvent.STIMULUS_TEST;
        endValue = "end";
        verbose = false;
        predEventType = NFEvent.CLASSIFIER_PREDICTION;
        baselineEventType = NFEvent.STIMULUS_BASELINE;
        baselineEnd = 1000;
        lengthTrMs = new double[]{};
        lengthTrSample = new double[]{};
        overlap = .5;
        stepMs = new double[]{};
        predFilt = new double[]{};
        timeout_ms = 1000;
        normLat = true;
    }

    @Override
    public void run() {

    }
}
