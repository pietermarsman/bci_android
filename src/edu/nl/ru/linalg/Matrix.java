package edu.nl.ru.linalg;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 * Created by Pieter on 27-1-2015.
 */
public class Matrix extends Array2DRowRealMatrix {
    // Methods for
    // [x] n-d arrays,
    // [x] matrix-matrix and                            > multiply
    // [x] matrix-vector products,                      > operate (M*v) and preMultiply (v*M)
    // [x] mean computation,                            > this.mean()
    // [x] sum-along a dimension,                       > this.sum()
    // [ ] FFT,
    // [ ] inverse FFT,
    // [ ] temporal filtering,
    // [x] convolution,
    // [ ] detrending,
    // [ ] covariance computation,
    // [ ] data selection (subsetting),
    // [ ] outlier-detection,
    // [ ] welch's method for spectrum estimation,
    // [ ] fft-based spectral filtering,
    // [ ] eigen-decompositions (SVD/EIG).


    public Matrix() {
        super();
    }

    public Matrix(double[] v) {
        super(v);
    }

    public Matrix(double[][] d) {
        super(d);
    }

    public Matrix(double[][] d, boolean copyArray) {
        super(d, copyArray);
    }

    public Matrix(int rowDimension, int columnDimension) {
        super(rowDimension, columnDimension);
    }

    public Matrix(RealMatrix m) {
        super(m.getData());
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append(this.getClass().getSimpleName()).append("\n");
        for (int i = 0; i < getRowDimension(); i++) {
            for (int j = 0; j < getColumnDimension(); j++)
                sb.append(getData()[i][j]).append(" ");
            if (i < getRowDimension() - 1)
                sb.append("\n");
        }
        return sb.toString();
    }

    public int getDimension(int axis) {
        Matrix.checkAxis(axis);
        return axis == 0 ? this.getRowDimension() : this.getColumnDimension();
    }

    public Matrix round(int decimals) {
        if (decimals < 0)
            throw new IllegalArgumentException("The number of decimals should be bigger than 0");
        double[][] data = this.getData();
        double factor = Math.pow(10, decimals);
        for (int i = 0; i < data.length; i++)
            for (int j = 0; j < data[i].length; j++)
                data[i][j] = Math.round(data[i][j] * factor) / factor;
        return new Matrix(data);
    }

    public Matrix repeat(int times, int axis) {
        if (times < 1)
            throw new IllegalArgumentException("Times should be bigger than 0 but it is " + times);
        Matrix.checkAxis(axis);
        int rows = this.getRowDimension();
        int columns = this.getColumnDimension();
        if (axis == 0)
            rows = rows * times;
        else
            columns = columns * times;
        Matrix repeated;
        if (axis == 0) {
            double[][] newData = new double[rows][columns];
            for (int r = 0; r < this.getRowDimension(); r++)
                for (int t = 0; t < times; t++)
                    newData[t + times * r] = this.getRow(r);
            repeated = new Matrix(newData);
        } else {
            double[][] newData = new double[columns][rows];
            for (int c = 0; c < this.getColumnDimension(); c++)
                for (int t = 0; t < times; t++)
                    newData[t + times * c] = this.getColumn(c);
            repeated = new Matrix(new Matrix(newData).transpose());
        }
        return repeated;
    }

    public Matrix mean() {
        return mean(-1);
    }

    public Matrix mean(int axis) {
        double scalar;
        if (axis == 0)
            scalar = this.getRowDimension();
        else if (axis == 1)
            scalar = this.getColumnDimension();
        else if (axis == -1)
            scalar = this.getRowDimension() * this.getColumnDimension();
        else
            throw new IllegalArgumentException("Wrong axis selected. Should be either -1, 0 or 1 but is " + axis);
        scalar = 1.0 / scalar;
        return new Matrix(sum(axis).scalarMultiply(scalar));
    }

    public Matrix sum() {
        return sum(-1);
    }

    public Matrix sum(int axis) {
        double[][] data = this.getData();
        if (axis < 0) {
            double[] mean = {0.0};
            for (double[] row : data) {
                for (double elem : row) {
                    mean[0] += elem;
                }
            }

            return new Matrix(mean);
        } else if (axis == 0) {
            double[] mean = new double[this.getColumnDimension()];
            for (int i = 0; i < this.getRowDimension(); i++) {
                for (int j = 0; j < this.getColumnDimension(); j++) {
                    mean[j] += data[i][j];
                }
            }
            return new Matrix(mean);
        } else if (axis == 1) {
            double[] mean = new double[this.getRowDimension()];
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[i].length; j++) {
                    mean[i] += data[i][j];
                }
            }
            return new Matrix(mean);
        } else {
            throw new IllegalArgumentException("Wrong axis selected. Should be either -1, 0 or 1 but is " + axis);
        }
    }

    public Matrix covariance() {
        Covariance cov = new Covariance(this, true);
        return new Matrix(cov.getCovarianceMatrix());

    }

    public Matrix detrend(int axis, String type) {
        Matrix.checkAxis(axis);
        if (!(type.equalsIgnoreCase("linear") || type.equalsIgnoreCase("constant")))
            throw new IllegalArgumentException("Type should be `linear' or `constant' but is `" + type + "'");
        if (type.equalsIgnoreCase("constant")) {
            Matrix mean = this.mean(axis);
            int otherAxis = axis == 0 ? 1 : 0;
            mean = mean.repeat(this.getDimension(axis), 1);
            if (axis == 0)
                mean = new Matrix(mean.transpose());
            return new Matrix(this.subtract(mean));
        } else {
            double[][] ret = this.getData();
            if (axis == 0) {
                for (int c = 0; c < this.getColumnDimension(); c++) {
                    SimpleRegression regression = new SimpleRegression();
                    double[] column = this.getColumn(c);
                    for (int i = 0; i < column.length; i++)
                        regression.addData(i, column[i]);
                    double slope = regression.getSlope();
                    double intercept = regression.getIntercept();
                    for (int r = 0; r < this.getRowDimension(); r++)
                        ret[r][c] = ret[r][c] - intercept - r * slope;
                }
            } else {
                for (int r = 0; r < this.getRowDimension(); r++) {
                    SimpleRegression regression = new SimpleRegression();
                    double[] row = this.getRow(r);
                    for (int i = 0; i < row.length; i++)
                        regression.addData(i, row[i]);
                    double slope = regression.getSlope();
                    double intercept = regression.getIntercept();
                    for (int c = 0; c < this.getColumnDimension(); c++)
                        ret[r][c] = ret[r][c] - intercept - c * slope;
                }
            }
            return new Matrix(ret);
        }
    }

    public static void checkAxis(int axis) throws IllegalArgumentException {
        if (axis < 0 || axis > 1)
            throw new IllegalArgumentException("Axis should be 0 or 1 but is " + axis);
    }
}
