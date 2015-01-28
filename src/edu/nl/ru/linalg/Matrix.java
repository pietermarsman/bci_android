package edu.nl.ru.linalg;

import edu.nl.ru.miscellaneous.Triple;
import edu.nl.ru.miscellaneous.Tuple;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import java.util.Arrays;

import static edu.nl.ru.miscellaneous.DoubleArrayFunctions.reverseDoubleArrayInPlace;
import static edu.nl.ru.miscellaneous.DoubleArrayFunctions.getSortIdx;
import static org.apache.commons.math3.stat.StatUtils.max;

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
    // [x] FFT,
    // [x] inverse FFT,
    // [x] temporal filtering,
    // [ ] convolution,
    // [x] detrending,
    // [x] covariance computation,
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

    public static void checkString(String given, String[] options) {
        if (!(Arrays.asList(options).contains(given))) {
            StringBuilder sb = new StringBuilder();
            sb.append("Type should be in [");
            for (String option : options)
                sb.append(option).append(", ");
            sb.append("] but is ").append(given);
            throw new IllegalArgumentException(sb.toString());
        }
    }

    public static void checkRepeats(int repeats) {
        if (repeats < 1)
            throw new IllegalArgumentException("Times should be bigger than 0 but it is " + repeats);
    }

    public static void checkDecimals(int decimals) {
        if (decimals < 0)
            throw new IllegalArgumentException("The number of decimals should be bigger than 0");
    }

    public static void checkAxis(int axis) throws IllegalArgumentException {
        if (axis < 0 || axis > 1)
            throw new IllegalArgumentException("Axis should be 0 or 1 but is " + axis);
    }

    public static Matrix zeros(int dim0, int dim1) {
        double[][] zeros = new double[dim0][dim1];
        return new Matrix(zeros);
    }

    public static Matrix ones(int dim0, int dim1) {
        double[][] zeros = new double[dim0][dim1];
        return new Matrix(new Array2DRowRealMatrix(zeros).scalarAdd(1.0));
    }

    public static Matrix eye(int dim) {
        double[] ones = new double[dim];
        for (int i = 0; i < dim; i++)
            ones[i] = 1.0;
        return new Matrix(new DiagonalMatrix(ones));
    }

    public static Matrix car(int size) {
        return new Matrix(Matrix.eye(size).scalarAdd(-1.0 / ((double) size)));
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
        Matrix.checkDecimals(decimals);
        double[][] data = this.getData();
        double factor = Math.pow(10, decimals);
        for (int i = 0; i < data.length; i++)
            for (int j = 0; j < data[i].length; j++)
                data[i][j] = Math.round(data[i][j] * factor) / factor;
        return new Matrix(data);
    }

    public Matrix flipUD() {
        double[][] newMatrix = this.getData();
        double[][] oldMatrix = this.getData();
        for (int r = 0; r < oldMatrix.length; r++)
            for (int c = 0; c < oldMatrix[r].length; c++)
                newMatrix[newMatrix.length - r - 1][c] = oldMatrix[r][c];
        return new Matrix(newMatrix);
    }

    public Matrix flipLR() {
        double[][] newMatrix = this.getData();
        double[][] oldMatrix = this.getData();
        for (int r = 0; r < oldMatrix.length; r++)
            for (int c = 0; c < oldMatrix[r].length; c++)
                newMatrix[r][newMatrix[r].length - c - 1] = oldMatrix[r][c];
        return new Matrix(newMatrix);
    }

    public Matrix repeat(int repeats, int axis) {
        Matrix.checkAxis(axis);
        Matrix.checkRepeats(repeats);
        int rows = this.getRowDimension();
        int columns = this.getColumnDimension();
        if (axis == 0)
            rows = rows * repeats;
        else
            columns = columns * repeats;
        Matrix repeated;
        if (axis == 0) {
            double[][] newData = new double[rows][columns];
            for (int r = 0; r < this.getRowDimension(); r++)
                for (int t = 0; t < repeats; t++)
                    newData[t + repeats * r] = this.getRow(r);
            repeated = new Matrix(newData);
        } else {
            double[][] newData = new double[columns][rows];
            for (int c = 0; c < this.getColumnDimension(); c++)
                for (int t = 0; t < repeats; t++)
                    newData[t + repeats * c] = this.getColumn(c);
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
        if (axis == -1) {
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
        Covariance cov = new Covariance(this.transpose(), true);
        return new Matrix(cov.getCovarianceMatrix());

    }

    public Matrix detrend(int axis, String type) {
        Matrix.checkAxis(axis);
        Matrix.checkString(type, new String[]{"constant", "linear"});
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

    public Matrix fft(int axis) {
        return fft(axis, TransformType.FORWARD);
    }

    public Matrix ifft(int axis) {
        return fft(axis, TransformType.INVERSE);
    }

    public Matrix fft(int axis, TransformType direction) {
        Matrix.checkAxis(axis);
        // FIXME which normalization to use?
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        double[][] ft = new double[this.getRowDimension()][this.getColumnDimension()];
        if (axis == 0) {
            for (int c = 0; c < this.getColumnDimension(); c++) {
                Complex[] complexResult = fft.transform(this.getColumn(c), direction);
                for (int i = 0; i < complexResult.length; i++)
                    ft[i][c] = Math.pow(complexResult[i].abs(), 2.0);
            }
        } else {
            for (int r = 0; r < this.getRowDimension(); r++) {
                Complex[] complexResult = fft.transform(this.getRow(r), direction);
                double[] power = new double[complexResult.length];
                for (int i = 0; i < complexResult.length; i++)
                    ft[r][i] = Math.pow(complexResult[i].abs(), 2.0);
            }
        }
        return new Matrix(ft);
    }

    public Tuple<Matrix, RealVector> eig() {
        return eig("descending");
    }

    public Tuple<Matrix, RealVector> eig(String order) {
        // FIXME is returning negative of python implementation
        // FIXME slightly different from python implementation
        Matrix.checkString(order, new String[]{"descending", "ascending"});
        EigenDecomposition eig = new EigenDecomposition(this);
        RealMatrix oldVecs = new Array2DRowRealMatrix(eig.getV().getData());
        double[] oldVals = eig.getRealEigenvalues();
        RealMatrix newVecs = new Array2DRowRealMatrix(eig.getV().getData());
        double[] newVals = eig.getRealEigenvalues();
        Integer[] sortIdx = getSortIdx(oldVals);
        for (int i = 0; i < sortIdx.length; i++) {
            Integer id = sortIdx[i];
            newVecs.setColumn(i, oldVecs.getColumn(id));
            newVals[i] = oldVals[id];
        }
        if (order.equalsIgnoreCase("descending")) {
            return new Tuple<Matrix, RealVector>(new Matrix(newVecs), new ArrayRealVector(newVals));
        } else {
            Matrix eigenVecs = new Matrix(newVecs).flipLR();
            reverseDoubleArrayInPlace(newVals);
            RealVector eigenValues = new ArrayRealVector(newVals);
            return new Tuple<Matrix, RealVector>(eigenVecs, eigenValues);
        }
    }

    public Triple<Matrix, Matrix, Matrix> svd() {
        SingularValueDecomposition svd = new SingularValueDecomposition(this);
        return new Triple<Matrix, Matrix, Matrix>(new Matrix(svd.getU()), new Matrix(svd.getS()), new Matrix(svd
                .getVT()));
    }

    public Matrix spatialFilter(String type) {
        return spatialFilter(type, 1e-15);
    }

    public Matrix spatialFilter(String type, double whitenThres) {
        Matrix.checkString(type, new String[]{"car", "whiten"});
        if (type.equalsIgnoreCase("car")) {
            return new Matrix(this.preMultiply(Matrix.car(this.getRowDimension())));
        } else {
            // Get eigen decomposition of the covariance matrix
            Tuple<Matrix, RealVector> eigenDecomposition = this.covariance().eig("ascending");
            Matrix eigenVecs = eigenDecomposition.x;
            double[] eigenVals = eigenDecomposition.y.toArray();
            // Use the decomposition to create the multiplication matrix for the whiten spatial filter
            double[] diagVals = new double[eigenVals.length];
            double max = max(diagVals);
            for (int i = 0; i < eigenVals.length; i++)
                diagVals[i] = eigenVals[i] > max * whitenThres ? Math.pow(eigenVals[i], -.5) : 0.0;
            RealMatrix diag = new DiagonalMatrix(diagVals);
            RealMatrix transform = eigenVecs.multiply(diag).multiply(eigenVecs.transpose());
            return new Matrix(this.preMultiply(transform));
        }
    }
}
