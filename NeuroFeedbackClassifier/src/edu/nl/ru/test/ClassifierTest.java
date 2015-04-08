package edu.nl.ru.test;

import edu.nl.ru.analysis.Classifier;
import edu.nl.ru.analysis.ClassifierResult;
import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.Windows;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.RealVector;

import java.util.LinkedList;
import java.util.List;

/**
 * Created by Pieter on 5-3-2015.
 * Testing the classifier
 */
public class ClassifierTest extends TestCase {

    private Matrix e;

    protected void setUp() throws Exception {
        //        double[][] dataA = {{1.0, 2.0}, {5.0, 4.0}};
        //        double[][] dataB = {{0.5, 0.4, 0.2}, {0.3, 0.2, 0.2}, {.2, .2, .7}};
        //        double[][] dataC = {{1.2, -1232.0}, {-67.5, .232}};
        //        double[][] dataD = {{1.2, 23., 12., 12.}, {12., 43., 432., 23.}, {3., 23.2, -12., -3.2}, {1., 1., 2
        // ., 2.}};
        double[][] dataE = {{13.21, 32., 432., .324, 43., .1}, {234., 56., 56.6, 765., 876., .1}, {345., 34., 2123.,
                76., 34., .2}};
        e = new Matrix(dataE);
    }

    public void testApply() throws Exception {
        double badChannelThreshold = -1.;
        Integer[] timeIdx = new Integer[]{0, 1, 2};
        Integer[] freqIdx = new Integer[]{0, 1};
        Matrix W = Matrix.zeros(5, 2);
        List<Matrix> Ws = new LinkedList<Matrix>();
        Ws.add(W);
        RealVector b = Matrix.zeros(5, 1).getColumnVector(0);
        Double[] startMs = new Double[]{0.};
        String[] spectrumDescription = new String[]{"alphaL", "alphaR", "baddness", "badChL", "badChR"};
        Integer[] isBad = new Integer[]{0, 0, 0};
        Classifier classifier = new Classifier(Ws, b, true, badChannelThreshold, Windows.WindowType.HANNING,
                WelchOutputType.AMPLITUDE, timeIdx, freqIdx, 1, null, null, 2, 100., startMs, spectrumDescription,
                isBad);
        ClassifierResult ret = classifier.apply(e);
    }
}
