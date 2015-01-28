package edu.nl.ru.miscellaneous;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by Pieter on 28-1-2015.
 */
public class DoubleArrayFunctions {

    public static void reverseDoubleArrayInPlace(double[] arr) {
        for (int i = 0; i < arr.length / 2; i++) {
            double temp = arr[i];
            arr[i] = arr[arr.length - i - 1];
            arr[arr.length - i - 1] = temp;
        }
    }

    public static double max(double[] arr, int from) {
        double max = Double.NEGATIVE_INFINITY;
        for (double val : arr)
            if (val > max)
                max = val;
        return max;
    }

    public static Integer[] getSortIdx(final double[] arr) {
        Integer[] idx = new Integer[arr.length];
        for (int i = 0; i < idx.length; i++)
            idx[i] = i;
        Comparator<Integer> comp = new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return arr[o1] > arr[o2] ? -1 : (arr[o1] == arr[o2] ? 0 : 1);
            }
        };
        Arrays.sort(idx, comp);
        return idx;
    }
}
