package bmird.radboud.fieldtripclientsservice.threads;

import android.util.Log;
import bmird.radboud.fieldtripclientsservice.base.Argument;
import bmird.radboud.fieldtripclientsservice.base.ThreadBase;
import com.interaxon.libmuse.*;
import nl.fcdonders.fieldtrip.bufferclient.BufferClient;
import nl.fcdonders.fieldtrip.bufferclient.DataType;
import nl.fcdonders.fieldtrip.bufferclient.Header;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

/**
 * Created by Pieter Marsman on 6-4-2015.
 */
public class MuseConnection extends ThreadBase {

    public final String TAG = MuseConnection.class.toString();

    private int BUFFERSIZE = 65500;

    private BufferClient client = null;
    private Muse muse = null;
    private ConnectionListener connectionListener = null;
    private DataListener dataListener = null;
    private boolean dataTransmission = true;
    private int nSamples;

    private String address;
    private int port;
    private int channels;
    private int samplingFrequency;
    private int dataType;

    @Override
    public Argument[] getArguments() {
        Argument[] arguments = new Argument[5];
        arguments[0] = new Argument("Buffer Address", "localhost");
        arguments[1] = new Argument("Buffer port", 1972, true);
        arguments[2] = new Argument("Channels", 4, true);
        arguments[3] = new Argument("Sampling frequency", 220, true);
        arguments[4] = new Argument("Datatype", DataType.FLOAT64, true);
        return arguments;
    }

    @Override
    public String getName() {
        return "MuseConnection";
    }

    @Override
    public void mainloop() {
        // Create listeners and pass reference to activity to them
        address = arguments[0].getString();
        port = arguments[1].getInteger();
        channels = arguments[2].getInteger();
        samplingFrequency = arguments[3].getInteger();
        dataType = arguments[4].getInteger();
        android.updateStatus("Address: " + address + ":" + String.valueOf(port));
        Log.i(TAG, "Buffer server: " + address + " : " + port);

        client = new BufferClient();
        connectionListener = new ConnectionListener();
        dataListener = new DataListener();
        Log.i("Muse Headband", "libmuse version=" + LibMuseVersion.SDK_VERSION);

        getMuse();
        connectToMuse();

        run = true;
        connectToBuffer();
        uploadHeaderToBuffer(channels, samplingFrequency, dataType);

        long startMs = Calendar.getInstance().getTimeInMillis();
        long elapsedMs = 0;
        nSamples = 0;

        while (run) {
            long now = Calendar.getInstance().getTimeInMillis();
            if (startMs + elapsedMs > now + 5000) {
                android.updateStatus(nSamples + "(" + elapsedMs / 1000 + ")");
                Log.i(TAG, "Elapsed time: " + elapsedMs + ". nSamples: " + nSamples);
            }
        }
    }

    private void getMuse() {
        while (muse == null) {
            Log.i(TAG, "Refreshing paired muses list");
            MuseManager.refreshPairedMuses();
            List<Muse> pairedMuses = MuseManager.getPairedMuses();
            if (pairedMuses.size() > 0) muse = pairedMuses.get(0);
        }
    }

    private void connectToMuse() {
        ConnectionState state = muse.getConnectionState();
        if (state == ConnectionState.CONNECTED || state == ConnectionState.CONNECTING) {
            Log.w("Muse Headband", "doesn't make sense to connect second time to the same muse");
            return;
        }
        configureLibrary();
        /**
         * In most cases libmuse native library takes care about
         * exceptions and recovery mechanism, but native code still
         * may throw in some unexpected situations (like bad bluetooth
         * connection). Print all exceptions here.
         */
        try {
            muse.runAsynchronously();
        } catch (Exception e) {
            Log.e("Muse Headband", Log.getStackTraceString(e));
        }
    }

    private void disconectMuse() {
        if (muse != null) {
            /**
             * true flag will force libmuse to unregister all listeners,
             * BUT AFTER disconnecting and sending disconnection event.
             * If you don't want to receive disconnection event (for ex.
             * you call disconnect when application is closed), then
             * unregister listeners first and then call disconnect:
             * muse.unregisterAllListeners();
             * muse.disconnect(false);
             */
            muse.disconnect(true);
        }
    }

    private void configureLibrary() {
        muse.registerConnectionListener(connectionListener);
        muse.registerDataListener(dataListener, MuseDataPacketType.EEG);
        muse.setPreset(MusePreset.PRESET_14);
        muse.enableDataTransmission(dataTransmission);
    }

    private void connectToBuffer() {
        while (!client.isConnected()) {
            Log.i(TAG, "Connecting to " + address + ":" + port);
            android.updateStatus("Connecting to " + address + ":" + port);
            try {
                client.connect(address, port);
            } catch (IOException ex) {
            }
            if (!client.isConnected()) {
                android.updateStatus("Couldn't connect. Waiting");
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private void uploadHeaderToBuffer(int channels, int samplingFrequency, int dataType) {
        Header hdr = new Header(channels, samplingFrequency, dataType);
        try {
            client.putHeader(hdr);
        } catch (IOException e) {
            e.printStackTrace();
        }
        Log.i(TAG, "Uploaded header to the buffer");
        Log.i(TAG, hdr.toString());
    }

    @Override
    public void validateArguments(Argument[] arguments) {

    }

    @Override
    public void stop() {
        disconectMuse();
        if (client != null) try {
            client.disconnect();
        } catch (IOException e) {
            Log.e(TAG, Log.getStackTraceString(e));
        } finally {
            client = null;
        }

        super.stop();
    }

    class ConnectionListener extends MuseConnectionListener {

        @Override
        public void receiveMuseConnectionPacket(MuseConnectionPacket p) {
            final ConnectionState current = p.getCurrentConnectionState();
            final String status = p.getPreviousConnectionState().toString() +
                    " -> " + current;
            final String full = "Muse " + p.getSource().getMacAddress() + " " + status;
            Log.i("Muse Headband", full);
        }
    }

    /**
     * Data listener will be registered to listen for: Accelerometer,
     * Eeg and Relative Alpha bandpower packets. In all cases we will
     * update UI with new values.
     * We also will log message if Artifact packets contains "blink" flag.
     * DataListener methods will be called from execution thread. If you are
     * implementing "serious" processing algorithms inside those listeners,
     * consider to create another thread.
     */
    class DataListener extends MuseDataListener {

        public final String TAG = DataListener.class.toString();

        @Override
        public void receiveMuseDataPacket(MuseDataPacket p) {
            switch (p.getPacketType()) {
                case EEG:
                    updateEeg(p.getValues());
                    break;
                case ACCELEROMETER:
                    updateAccelerometer(p.getValues());
                    break;
                case ALPHA_RELATIVE:
                    updateAlphaRelative(p.getValues());
                    break;
                default:
                    break;
            }
        }

        @Override
        public void receiveMuseArtifactPacket(MuseArtifactPacket p) {
            if (p.getHeadbandOn() && p.getBlink()) {
                Log.i("Artifacts", "blink");
            }
        }

        private void updateAccelerometer(final ArrayList<Double> data) {
            Log.v(TAG, "Accelerometer: " + data.toString());
        }

        private void updateEeg(final ArrayList<Double> data) {
            Log.v(TAG, "EEG: " + data.toString());
            double[] dataArray = new double[data.size()];
            for (int i = 0; i < data.size(); i++)
                dataArray[i] = data.get(i);
            double[][] dataMatrix = new double[][]{dataArray};
            try {
                client.putData(dataMatrix);
            } catch (IOException e) {
                Log.e(TAG, Log.getStackTraceString(e));
            }
            nSamples += 1;
        }

        private void updateAlphaRelative(final ArrayList<Double> data) {
            Log.v(TAG, "Alpha: " + data.toString());
        }
    }

}
