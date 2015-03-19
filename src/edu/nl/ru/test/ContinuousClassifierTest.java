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
 */
public class ContinuousClassifierTest extends TestCase {

    /*
    Start fileplayback
    "C:\Program Files\Java\jdk1.7.0_02\bin\java" -Didea.launcher.port=7532 "-Didea.launcher.bin.path=C:\Program Files
     (x86)\JetBrains\IntelliJ IDEA Community Edition 14.0.2\bin" -Dfile.encoding=UTF-8 -classpath "C:\Program
     Files\Java\jdk1.7.0_02\jre\lib\charsets.jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\deploy.jar;C:\Program
     Files\Java\jdk1.7.0_02\jre\lib\javaws.jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\jce.jar;C:\Program
     Files\Java\jdk1.7.0_02\jre\lib\jsse.jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\management-agent.jar;
     C:\Program Files\Java\jdk1.7.0_02\jre\lib\plugin.jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\resources.jar;
     C:\Program Files\Java\jdk1.7.0_02\jre\lib\rt.jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\ext\dnsns.jar;
     C:\Program Files\Java\jdk1.7.0_02\jre\lib\ext\localedata.jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\ext\sunec
     .jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\ext\sunjce_provider.jar;C:\Program Files\Java\jdk1.7
     .0_02\jre\lib\ext\sunmscapi.jar;C:\Program Files\Java\jdk1.7.0_02\jre\lib\ext\zipfs.jar;
     C:\Users\Pieter\Documents\Projects\buffer_bci\java\out\production\java;
     C:\Users\Pieter\Documents\Projects\buffer_bci\dataAcq\buffer\java\BufferClient.jar;C:\Program Files (x86)
     \JetBrains\IntelliJ IDEA Community Edition 14.0.2\lib\idea_rt.jar" com.intellij.rt.execution.application.AppMain
      filePlayback
     */

    private static List<Matrix> loadWFromFile(String file, int rows, int columns) {
        List<Matrix> matrices = new LinkedList<Matrix>();
        BufferedReader br = null;
        String line = "";
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
        Classifier classifier = new Classifier(Ws, b, true, null, null, Windows.WindowType.HANNING,
                WelchOutputType.AMPLITUDE, null, freqIdx, 1, null, null, 128, 100., new Double[]{0.},
                spectrumDescription, isBad);
        classifiers.add(classifier);
        ContinuousClassifier c = new ContinuousClassifier("localhost", 1973, null, "stimulus.test", "end",
                "classifier.prediction", null, null, .5, 1000, classifiers, null);
        Thread t = new Thread(c);
        t.start();
        Thread.sleep(10 * 1000);
    }
}
