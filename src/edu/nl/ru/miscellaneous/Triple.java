package edu.nl.ru.miscellaneous;

/**
 * Created by Pieter on 28-1-2015.
 */
public class Triple<X, Y, Z> {
    public final X x;
    public final Y y;
    public final Z z;

    public Triple(X x, Y y, Z z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}