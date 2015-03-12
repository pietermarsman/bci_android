package edu.nl.ru.test;

import edu.nl.ru.analysis.Classifier;
import edu.nl.ru.analysis.ContinuousClassifier;
import edu.nl.ru.analysis.NFEvent;
import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.Windows;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.RealVector;
import org.apache.log4j.Level;

/**
 * Created by Pieter on 12-3-2015.
 */
public class ContinousClassifierTest extends TestCase {

    public void setUp() throws Exception {

    }

    public void testRun() throws Exception {
        double badChannelThreshold = -1.;
        double badTrialThreshold = -1.;
        int[] timeIdx = new int[]{0, 1, 2};
        int[] freqIdx = new int[]{0, 1};
        Matrix W = Matrix.zeros(2, 2);
        RealVector b = Matrix.zeros(2, 1).getColumnVector(0);
        double[] startMs = new double[]{0};
        Classifier classfier = new Classifier(W, b, true, badChannelThreshold, badTrialThreshold, Windows.WindowType
                .HANNING, WelchOutputType.AMPLITUDE, Level.ALL, timeIdx, freqIdx, 1, null, null, 2, 100, startMs);
        ContinuousClassifier c = new ContinuousClassifier("localhost", 1973, null, NFEvent.STIMULUS_TEST, "end",
                "classifier.prediction", 25, .5, 1000, classfier, 1.);
        Thread t = new Thread(c);
        t.start();
        Thread.sleep(10 * 1000);
    }
}
