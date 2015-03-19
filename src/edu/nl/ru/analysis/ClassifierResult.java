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

    public ClassifierResult(ClassifierResult classifierResult) {
        this.f = classifierResult.f;
        this.fraw = classifierResult.fraw;
        this.p = classifierResult.p;
        this.X = classifierResult.X;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Classifier result: ");
        sb.append("f").append(f.shapeString());
        sb.append(", fraw").append(fraw.shapeString());
        sb.append(", p").append(p.shapeString());
        sb.append(", X").append(X.shapeString());
        return sb.toString();
    }
}
