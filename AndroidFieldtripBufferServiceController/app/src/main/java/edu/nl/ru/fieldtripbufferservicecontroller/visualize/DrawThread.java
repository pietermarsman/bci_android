package edu.nl.ru.fieldtripbufferservicecontroller.visualize;

import android.content.Context;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.os.Handler;
import android.util.Log;
import android.view.SurfaceHolder;

/**
 * Created by pieter on 18-5-15.
 */
public class DrawThread extends Thread {

    private static final String TAG = DrawThread.class.getSimpleName();
    private final Paint paint = new Paint(Paint.ANTI_ALIAS_FLAG);
    private final SurfaceHolder sh;
    private final Context ctx;
    private int canvasWidth;
    private int canvasHeight;
    private boolean run = false;

    private float bubbleX;
    private float bubbleY;
    private int color;
    private float size;
    private float[] max;

    float baddnessFilter;

    private BufferThread bufferThread;

    public DrawThread(SurfaceHolder surfaceHolder, Context context, Handler handler, BufferThread bufferThread) {
        sh = surfaceHolder;
        handler = handler;
        ctx = context;
        max = new float[3];
        this.bufferThread = bufferThread;
        baddnessFilter = 0.0f;
    }

    public void initializeModel() {
        synchronized (sh) {
            // Start bubble in centre and create some random motion
            bubbleX = canvasWidth / 2;
            bubbleY = canvasHeight / 2;
        }
        updateModel();
    }

    public void updateModel() {
        synchronized (sh) {
            float[] values = bufferThread.getValues();
            color = computeColor(values);
            size = computeSize(values);
        }
    }

    private int computeColor(float values[]) {
        baddnessFilter = (float) (0.3 * values[1] + 0.7 * baddnessFilter);
        int mean = 128;
        int std = 64;
        int red = (int) (mean + baddnessFilter * std);
        red = Math.max(Math.min(red, 255), 0);
        int green = 255 - red;
        return Color.rgb(red, green, 0);
    }

    private int computeSize(float values[]) {
        float meanSize = 200.f;
        float stdSize = 200.f;
        float newSize = meanSize - values[0] * stdSize;
        newSize = Math.max(Math.min(newSize, Math.min(canvasHeight, canvasWidth)), 20.f);
        return (int) newSize;
    }

    public void run() {
        while (run) {
            updateModel();
            Canvas c = null;
            try {
                c = sh.lockCanvas(null);
                synchronized (sh) {
                    doDraw(c);
                }
            } finally {
                if (c != null) {
                    sh.unlockCanvasAndPost(c);
                }
            }
        }
    }

    public void setRunning(boolean b) {
        run = b;
    }

    public void setSurfaceSize(int width, int height) {
        synchronized (sh) {
            canvasWidth = width;
            canvasHeight = height;
            initializeModel();
        }
    }

    private void doDraw(Canvas canvas) {
        String rgb = "(" + Color.red(color) + ", " + Color.green(color) + ", " + Color.blue(color) + ")";
        Log.v(TAG, "Draw circle at (" + bubbleX + ", " + bubbleY + " with size " + size + " and color " + rgb);
        canvas.save();
        canvas.restore();
        canvas.drawColor(Color.BLACK);
        paint.setColor(color);
        canvas.drawCircle(bubbleX, bubbleY, size, paint);
    }
}