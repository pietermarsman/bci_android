package edu.nl.ru.linalg;

/**
 * Created by Pieter on 28-1-2015.
 */
public class Miscelancous {

    public static void reverseDoubleArrayInPlace(double[] arr) {
        for (int i = 0; i < arr.length / 2; i++) {
            double temp = arr[i];
            arr[i] = arr[arr.length - i - 1];
            arr[arr.length - i - 1] = temp;
        }
    }

    private static double max(double[] arr) {
        double max = Double.NEGATIVE_INFINITY;
        for (double val : arr)
            if (val > max)
                max = val;
        return max;
    }
}
