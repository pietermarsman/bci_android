package edu.nl.ru.test;

import edu.nl.ru.analysis.Classifier;
import edu.nl.ru.analysis.ClassifierResult;
import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.Windows;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.RealVector;

/**
 * Created by Pieter on 5-3-2015.
 */
public class ClassifierTest extends TestCase {

    Matrix a, b, c, d, e;

    protected void setUp() throws Exception {
        double[][] dataA = {{1.0, 2.0}, {5.0, 4.0}};
        double[][] dataB = {{0.5, 0.4, 0.2}, {0.3, 0.2, 0.2}, {.2, .2, .7}};
        double[][] dataC = {{1.2, -1232.0}, {-67.5, .232}};
        double[][] dataD = {{1.2, 23., 12., 12.}, {12., 43., 432., 23.}, {3., 23.2, -12., -3.2}, {1., 1., 2., 2.}};
        double[][] dataE = {{13.21, 32., 432., .324, 43., .1}, {234., 56., 56.6, 765., 876., .1}, {345., 34., 2123.,
                76., 34., .2}};
        a = new Matrix(dataA);
        b = new Matrix(dataB);
        c = new Matrix(dataC);
        d = new Matrix(dataD);
        e = new Matrix(dataE);
    }

    public void testApply() throws Exception {
        double badChannelThreshold = -1.;
        double badTrialThreshold = -1.;
        int[] timeIdx = new int[]{0, 1, 2};
        int[] freqIdx = new int[]{0, 1};
        Matrix W = Matrix.zeros(5, 2);
        RealVector b = Matrix.zeros(5, 1).getColumnVector(0);
        double[] startMs = new double[]{0};
        Classifier classfier = new Classifier(W, b, true, badChannelThreshold, badTrialThreshold, Windows.WindowType
                .HANNING, WelchOutputType.AMPLITUDE, true, timeIdx, freqIdx, 1, null, null, 2, 100, startMs);
        ClassifierResult ret = classfier.apply(e);
    }
}
