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

    private BufferThread bufferThread;

    public DrawThread(SurfaceHolder surfaceHolder, Context context, Handler handler, BufferThread bufferThread) {
        sh = surfaceHolder;
        handler = handler;
        ctx = context;
        max = new float[3];
        this.bufferThread = bufferThread;
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
        max[0] = Math.max(values[1], max[0]);
        max[1] = Math.max(values[2], max[1]);
        max[2] = Math.max(values[3], max[2]);
        float[] colors = new float[]{values[1] / max[0] * 360, values[2] / max[1], values[3] / max[2]};
        return Color.HSVToColor(colors);
    }

    private int computeSize(float values[]) {
        return (int) (values[0] * 800.);
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