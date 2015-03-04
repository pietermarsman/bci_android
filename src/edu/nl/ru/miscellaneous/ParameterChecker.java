package edu.nl.ru.miscellaneous;

import java.util.Arrays;

/**
 * Created by Pieter on 23-2-2015.
 */
public class ParameterChecker {
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

    public static void checkNonNegative(int decimals) {
        if (decimals < 0)
            throw new IllegalArgumentException("The number of decimals should be bigger than 0");
    }

    public static void checkNonZero(double number) {
        if (number == 0.0)
            throw  new IllegalArgumentException("Number should not be zoro.");
    }

    public static void checkAxis(int axis) throws IllegalArgumentException {
        checkAxis(axis, false);
    }

    public static void checkAxis(int axis, boolean allowMinOne) throws IllegalArgumentException {
        if (axis < -1 || axis > 1)
            if (!(allowMinOne && axis == -1))
                throw new IllegalArgumentException("Axis should be 0 or 1 but is " + axis);
    }

    public static void checkLowerUpperThreshold(double lowerThreshold, double upperThreshold) throws
            IllegalArgumentException {
        if (lowerThreshold > upperThreshold)
            throw new IllegalArgumentException("Lower threshold (=" + lowerThreshold + ") should be lower than upper " +
                    "threshold (=" + upperThreshold + ")");
    }

    public static void checkEquals(int a, int b) throws  IllegalArgumentException{
        if (a != b)
            throw new IllegalArgumentException("Should be equal but are not: " + a + " and " + b);
    }
}