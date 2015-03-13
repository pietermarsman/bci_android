package edu.nl.ru.analysis;

import edu.nl.ru.linalg.Matrix;

/**
 * Created by Pieter on 23-2-2015.
 */
public class ClassifierResult {

    public final Matrix f, fraw, p, X;

    public ClassifierResult(Matrix f, Matrix fraw, Matrix p, Matrix X) {
        this.f = f;
        this.fraw = fraw;
        this.p = p;
        this.X = X;
    }
}
