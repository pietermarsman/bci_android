package edu.nl.ru.linalg;

import edu.nl.ru.miscellaneous.ExtraMath;
import edu.nl.ru.miscellaneous.ParameterChecker;
import edu.nl.ru.miscellaneous.Triple;
import edu.nl.ru.miscellaneous.Tuple;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.AbstractUnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.MathArrays;
import org.apache.log4j.Logger;

import java.util.ArrayList;

import static edu.nl.ru.miscellaneous.DoubleArrayFunctions.getSortIdx;
import static edu.nl.ru.miscellaneous.DoubleArrayFunctions.reverseDoubleArrayInPlace;
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
    // [x] FFT,                                         > this.fft()
    // [x] inverse FFT,                                 > this.ifft()
    // [x] temporal filtering,                          > this.spatialFilter()
    // [x] convolution,                                 > this.convolve()
    // [x] detrending,                                  > this.detrend()
    // [x] covariance computation,                      > this.covariance()
    // [x] data selection (subsetting)                  > getSubMatrix()
    // [x] outlier-detection,                           > removeOutliers()
    // [x] welch's method for spectrum estimation,      > welch()
    // [ ] fft-based spectral filtering,
    // [x] eigen-decompositions (SVD/EIG)               > this.eig() and this.svd()

    static Logger log = Logger.getLogger(Matrix.class);

    public Matrix() {
        super();
        log.setLevel(null);
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

    public static Matrix zeros(int dim0, int dim1) {
        ParameterChecker.checkNonNegative(dim0);
        ParameterChecker.checkNonNegative(dim1);

        double[][] zeros = new double[dim0][dim1];
        return new Matrix(zeros);
    }

    public static Matrix ones(int dim0, int dim1) {
        ParameterChecker.checkNonNegative(dim0);
        ParameterChecker.checkNonNegative(dim1);

        double[][] zeros = new double[dim0][dim1];
        return new Matrix(new Array2DRowRealMatrix(zeros).scalarAdd(1.0));
    }

    public static Matrix eye(int dim) {
        ParameterChecker.checkNonNegative(dim);

        double[] ones = new double[dim];
        for (int i = 0; i < dim; i++)
            ones[i] = 1.0;
        return new Matrix(new DiagonalMatrix(ones));
    }

    public static Matrix car(int size) {
        ParameterChecker.checkNonNegative(size);

        return new Matrix(Matrix.eye(size).scalarAdd(-1.0 / ((double) size)));
    }

    public static int[] range(int start, int end, int step) {
        ParameterChecker.checkNonZero(step);
        int size = (end - start) / step;
        int[] arr = new int[size];
        int index = 0;
        for (int i = start; i < end; i += step) {
            arr[index] = i;
            index++;
        }
        return arr;
    }

    public static Matrix range(double start, double end, double step) {
        ParameterChecker.checkNonZero(step);

        int size = (int) ((end - start) / step);
        double[] arr = new double[size];
        int index = 0;
        for (double i = start; i < end; i += step)
            arr[index] = i;
        return new Matrix(arr);
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

    public String shapeString() {
        return "(" + this.getRowDimension() + ", " + this.getColumnDimension() + ")";
    }

    public int getDimension(int axis) {
        ParameterChecker.checkAxis(axis);

        return axis == 0 ? this.getRowDimension() : this.getColumnDimension();
    }

    public Matrix round(int decimals) {
        ParameterChecker.checkNonNegative(decimals);

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
        ParameterChecker.checkAxis(axis);
        ParameterChecker.checkRepeats(repeats);

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

    public Matrix abs() {
        Matrix m = new Matrix(this.copy());
        m.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return Math.abs(value);
            }
        });
        return m;
    }

    public Matrix sqrt() {
        Matrix m = new Matrix(this.copy());
        m.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return Math.sqrt(value);
            }
        });
        return m;
    }

    public Matrix mean() {
        return mean(-1);
    }

    public Matrix mean(int axis) {
        ParameterChecker.checkAxis(axis, true);

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

    public Matrix median(int axis) {
        Median med = new Median();
        return this.evaluateUnivariateStatistic(axis, med);
    }

    public Matrix multiplyElements(final Matrix b) {
        ParameterChecker.checkEquals(this.getRowDimension(), b.getRowDimension());
        ParameterChecker.checkEquals(this.getColumnDimension(), b.getColumnDimension());

        RealMatrix c = this.copy();
        c.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            public double visit(int row, int column, double value) {
                return value * b.getEntry(row, column);
            }
        });
        return new Matrix(c);
    }

    public Matrix evaluateUnivariateStatistic(int axis, AbstractUnivariateStatistic stat) {
        ParameterChecker.checkAxis(axis, true);

        double[] data;
        if (axis == -1) {
            return this.flatten().evaluateUnivariateStatistic(0, stat);
        } else if (axis == 0) {
            data = new double[this.getColumnDimension()];
            for (int c = 0; c < this.getColumnDimension(); c++)
                data[c] = stat.evaluate(this.getColumn(c));
        } else {
            data = new double[this.getRowDimension()];
            for (int r = 0; r < this.getRowDimension(); r++)
                data[r] = stat.evaluate(this.getRow(r));
        }
        return new Matrix(data);
    }

    public Matrix sum() {
        return sum(-1);
    }

    public Matrix sum(int axis) {
        ParameterChecker.checkAxis(axis, true);

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

    public Matrix variance(int axis) {
        Variance var = new Variance(false);
        return this.evaluateUnivariateStatistic(axis, var);
    }

    public Matrix std(int axis) {
        StandardDeviation std = new StandardDeviation();
        return this.evaluateUnivariateStatistic(axis, std);
    }

    public Matrix flatten() {
        double[] data = new double[this.getRowDimension() * this.getColumnDimension()];
        for (int r = 0; r < this.getRowDimension(); r++)
            for (int c = 0; c < this.getColumnDimension(); c++)
                data[r * this.getColumnDimension() + c] = this.getEntry(r, c);
        return new Matrix(data);
    }

    public Matrix detrend(int axis, String type) {
        ParameterChecker.checkAxis(axis);
        ParameterChecker.checkString(type, new String[]{"constant", "linear"});

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
        ParameterChecker.checkAxis(axis);

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
                for (int i = 0; i < complexResult.length; i++)
                    ft[r][i] = Math.pow(complexResult[i].abs(), 2.0);
            }
        }
        return new Matrix(ft);
    }

    public Complex[][] fftComplex(int axis, TransformType direction) {
        ParameterChecker.checkAxis(axis);

        // FIXME which normalization to use?
        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[][] ft = new Complex[this.getRowDimension()][this.getColumnDimension()];
        if (axis == 0) {
            for (int c = 0; c < this.getColumnDimension(); c++) {
                // FIXME use data lengths with powers of 2
                Complex[] complexResult = fft.transform(this.getColumn(c), direction);
                for (int i = 0; i < complexResult.length; i++)
                    ft[i][c] = complexResult[i];
            }
        } else {
            for (int r = 0; r < this.getRowDimension(); r++) {
                double row[] = new double[(int) Math.pow(Math.ceil(ExtraMath.log(this.getRowDimension(), 2)), 2)];
                for (int c = 0; c < this.getColumnDimension(); c++)
                    row[c] = this.getEntry(r, c);
                Complex[] fftResult = fft.transform(row, direction);
                for (int c = 0; c < this.getColumnDimension(); c++)
                    ft[r][c] = fftResult[c];
            }
        }
        return ft;
    }

    public Tuple<Matrix, RealVector> eig() {
        return eig("descending");
    }

    public Tuple<Matrix, RealVector> eig(String order) {
        ParameterChecker.checkString(order, new String[]{"descending", "ascending"});

        // FIXME is returning negative of python implementation
        // FIXME slightly different from python implementation
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
        ParameterChecker.checkString(type, new String[]{"car", "whiten"});
        ParameterChecker.checkNonNegative(whitenThres);

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

    public Matrix convolve(double[] function, int axis) {
        ParameterChecker.checkNonZero(function.length);
        ParameterChecker.checkAxis(axis);

        if (axis == 0) {
            int newLength = this.getColumnDimension() + function.length - 1;
            double[][] data = new double[this.getRowDimension()][newLength];
            for (int r = 0; r < this.getRowDimension(); r++) {
                data[r] = MathArrays.convolve(this.getRow(r), function);
            }
            return new Matrix(data);
        } else {
            return new Matrix(new Matrix(this.transpose()).convolve(function, 0).transpose());
        }
    }

    public Matrix removeOutliers(int axis, double lowerThreshold, double upperThreshold, int maxIter, String feat) {
        ParameterChecker.checkAxis(axis);
        ParameterChecker.checkNonZero(maxIter);
        ParameterChecker.checkNonNegative(maxIter);
        ParameterChecker.checkString(feat, new String[]{"var", "mu"});

        Matrix m = this;
        if (maxIter > 1)
            m = m.removeOutliers(axis, lowerThreshold, upperThreshold, maxIter - 1, feat);
        if (m == null)
            return m;

        Matrix feature;
        if (feat.equalsIgnoreCase("var")) {
            feature = m.variance(axis).abs().sqrt();
        } else {
            feature = m.mean(axis);
        }

        double median = feature.median(-1).getData()[0][0];
        double std = feature.std(-1).getData()[0][0];
        double high = median + upperThreshold * std;
        double low = median + lowerThreshold * std;

        int inlierCount = 0;
        for (double val : feature.getColumn(0)) {
            inlierCount += (low < val) && (val < high) ? 1 : 0;
        }

        if (inlierCount <= 0)
            return null;
        else {
            int index = 0;
            Matrix ret;
            if (axis == 0) {
                ret = new Matrix(m.getRowDimension(), inlierCount);
                for (int c = 0; c < m.getColumnDimension(); c++) {
                    double featureValue = feature.getColumn(0)[c];
                    if (low < featureValue && featureValue < high) {
                        ret.setColumn(index, m.getColumn(c));
                        index++;
                    }
                }
            } else {
                ret = new Matrix(inlierCount, m.getColumnDimension());
                for (int r = 0; r < m.getRowDimension(); r++) {
                    double featureValue = feature.getColumn(0)[r];
                    if (low < featureValue && featureValue < high) {
                        ret.setRow(index, m.getRow(r));
                        index++;
                    }
                }
            }

            return ret;
        }
    }

    public Matrix welch(final int dim, final double[] taper, int[] start, int width, boolean detrendp) {
        log.debug("Applying welch algorithm: dim=" + dim + ", width=" + width);
        log.debug("Checking parameters");
        ParameterChecker.checkAxis(dim);
        ParameterChecker.checkNonZero(taper.length);
        ParameterChecker.checkNonZero(start.length);
        ParameterChecker.checkNonZero(width);
        ParameterChecker.checkNonNegative(width);

        // TODO add outtype, default is 'amp'
        WelchOutputType outType = WelchOutputType.AMPLITUDE;

        // Abbreviations
        int otherDim = dim == 0 ? 1 : 0;
        int sizeDim = this.getDimension(dim);
        int sizeOtherDim = this.getDimension(otherDim);

        // Create W
        int wWidth, wHeight;
        int reducedDim = (int) Math.round(((Math.ceil(((double) width - 1) / 2) + 1)));
        if (dim == 0) {
            wHeight = reducedDim;
            wWidth = this.getColumnDimension();
        } else {
            wHeight = this.getRowDimension();
            wWidth = reducedDim;
        }
        Matrix W = Matrix.zeros(wHeight, wWidth);

        // Create indexes
        ArrayList<int[]> idx = new ArrayList<int[]>();
        idx.add(Matrix.range(0, this.getRowDimension(), 1));
        idx.add(Matrix.range(0, this.getColumnDimension(), 1));
        ArrayList<int[]> wIdx = new ArrayList<int[]>();
        wIdx.add(Matrix.range(0, wHeight, 1));
        wIdx.add(Matrix.range(0, wWidth, 1));

        // Sum over the windows
        for (int wi : start) {
            // Window the dimension
            idx.set(dim, Matrix.range(wi, wi + width, 1));

            // Get submatrix
            Matrix wX = new Matrix(this.getSubMatrix(idx.get(0), idx.get(1)));

            // TODO add centerp (subtracting mean from the sample p

            // Detrend submatrix
            if (detrendp)
                wX = wX.detrend(dim, "linear");

            // Window
            wX.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int row, int column, double value) {
                    if (dim == 0)
                        return value * taper[row];
                    else
                        return value * taper[column];
                }
            });

            // Fourier
            wX = new Matrix(wX.fft(dim).scalarMultiply(2.0));

            // Positive frequency only
            wX = new Matrix(wX.getSubMatrix(wIdx.get(0), wIdx.get(1)));

            switch (outType) {
                case AMPLITUDE:
                    W = new Matrix(W.add(wX.sqrt()));
                    break;
                default:
                    throw new IllegalArgumentException("Only amp is supported");
            }
        }
        W = new Matrix(W.scalarMultiply(1. / (start.length * new Matrix(taper).sum().getEntry(0, 0))));
        return W;
    }
}
