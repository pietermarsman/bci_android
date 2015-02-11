package edu.nl.ru.miscellaneous;

import edu.nl.ru.linalg.Matrix;
import org.apache.commons.math3.analysis.function.Gaussian;

/**
 * Created by Pieter on 11-2-2015.
 */
public class Windows {

    public static Matrix gaussianWindow(int size, double sigma) {
        double[] gaussian = new double[size];
        Gaussian distribution = new Gaussian(((double) size - 1.) / 2.0, sigma);
        for (int i = 0; i < gaussian.length; i++) {
            gaussian[i] = distribution.value(i);
        }
        Matrix ret = new Matrix(new Matrix(gaussian).scalarMultiply(1.0 / gaussian[(size - 1) / 2]));
        return ret;
    }
}
