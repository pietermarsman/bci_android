package edu.nl.ru.test;

import edu.nl.ru.analysis.Classifier;
import edu.nl.ru.analysis.ContinuousClassifier;
import edu.nl.ru.linalg.Matrix;
import edu.nl.ru.linalg.WelchOutputType;
import edu.nl.ru.miscellaneous.ArrayFunctions;
import edu.nl.ru.miscellaneous.Windows;
import junit.framework.TestCase;
import org.apache.commons.math3.linear.RealVector;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by Pieter on 12-3-2015.
 * Test of the continuous classifier
 */
public class ContinuousClassifierTest extends TestCase {

    /*
    Start fileplayback
\cmd
     */

    private static List<Matrix> loadWFromFile(String file, int rows, int columns) {
        List<Matrix> matrices = new LinkedList<Matrix>();
        BufferedReader br = null;
        String line;
        String cvsSplitBy = ",";

        try {
            br = new BufferedReader(new FileReader(file));
            while ((line = br.readLine()) != null) {
                // use comma as separator
                String[] items = line.split(cvsSplitBy);
                Matrix m = new Matrix(ArrayFunctions.fromString(items)).reshape(rows, columns);
                matrices.add(m);
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        return matrices;
    }

    public void testRun() throws Exception {
        List<Matrix> Ws = ContinuousClassifierTest.loadWFromFile("resources/W.csv", 3, 56);
        RealVector b = Matrix.zeros(5, 1).getColumnVector(0);
        //        Integer[] timeIdx = new Integer[]{0, 1, 2};
        Integer[] freqIdx = ArrayFunctions.toObjectArray(Matrix.range(0, 56, 1));
        //        Double[] startMs = new Double[]{0.};
        List<Classifier> classifiers = new LinkedList<Classifier>();
        String[] spectrumDescription = new String[]{"alphaL", "alphaR", "baddness", "badChL", "badChR"};
        Integer[] isBad = new Integer[]{0, 0, 0};
        Classifier classifier = new Classifier(Ws, b, true, null, Windows.WindowType.HANNING, WelchOutputType
                .AMPLITUDE, null, freqIdx, 1, null, null, 128, 100., new Double[]{0.}, spectrumDescription, isBad);
        classifiers.add(classifier);
        ContinuousClassifier c = new ContinuousClassifier("localhost", 1973, null, "stimulus.test", "end",
                "classifier.prediction", "stimulus.startbaseline", "end", "start", .5, 1000, classifiers, null, 128, null, true);
        Thread t = new Thread(c);
        t.start();
        Thread.sleep(1000 * 1000);
    }
}
