package edu.nl.ru.fieldtripbufferservicecontroller;

import android.app.Activity;
import android.content.BroadcastReceiver;
import android.content.Context;
import android.content.Intent;
import android.content.IntentFilter;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.TableLayout;
import android.widget.TableRow;
import android.widget.TextView;
import edu.nl.ru.fieldtripbufferservicecontroller.visualize.BubbleSurfaceView;
import edu.nl.ru.fieldtripclientsservice.ThreadInfo;
import edu.nl.ru.fieldtripclientsservice.base.Argument;
import edu.nl.ru.monitor.ClientInfo;

import java.util.HashMap;


public class MainActivity extends Activity {

    private static final int LARGEPRIME = 961748927;
    public static String TAG = MainActivity.class.getSimpleName();

    public ServerController serverController;
    public ClientsController clientsController;

    private BroadcastReceiver mMessageReceiver;

    // Gui
    private TextView textView;
    private TableLayout table;
    private HashMap<Integer, Integer> threadToView;
    private BubbleSurfaceView surfaceView;


    @Override
    protected void onCreate(final Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        setContentView(R.layout.main_activity);
        textView = (TextView) findViewById(R.id.textView);
        table = (TableLayout) findViewById(R.id.table);
        threadToView = new HashMap<Integer, Integer>();
        surfaceView = (BubbleSurfaceView) findViewById(R.id.surfaceView);

        serverController = new ServerController(this);
        clientsController = new ClientsController(this);

        if (savedInstanceState == null) {
            IntentFilter intentFilter = new IntentFilter(C.FILTER_FROM_SERVER);
            mMessageReceiver = new BroadcastReceiver() {
                @Override
                public void onReceive(final Context context, Intent intent) {
                    updateServerController(intent);
                }
            };
            this.registerReceiver(mMessageReceiver, intentFilter);
        }
    }

    @Override
    public void onStop() {
        super.onStop();
        this.unregisterReceiver(mMessageReceiver);
        stopClients(null);
        stopServer(null);
    }

    private void updateServerController(Intent intent) {
        if (intent.getBooleanExtra(C.IS_BUFFER_INFO, false)) {
            updateServerInfo(intent);
        }
        if (intent.getBooleanExtra(C.IS_CLIENT_INFO, false)) {
            updateClientInfoFromServer(intent);
        }
        if (intent.getBooleanExtra(C.IS_THREAD_INFO, false)) {
            updateThreadsInfo(intent);
        }
    }

    private void updateServerInfo(Intent intent) {
        serverController.buffer = intent.getParcelableExtra(C.BUFFER_INFO);
        if (!serverController.initalUpdateCalled) {
            serverController.initialUpdate();
        }
        serverController.updateBufferInfo();
        textView.setText(serverController.toString());
    }

    private void updateClientInfoFromServer(Intent intent) {
        int numOfClients = intent.getIntExtra(C.CLIENT_N_INFOS, 0);
        ClientInfo[] clientInfo = new ClientInfo[numOfClients];
        for (int k = 0; k < numOfClients; ++k) {
            clientInfo[k] = intent.getParcelableExtra(C.CLIENT_INFO + k);
        }
        serverController.updateClients(clientInfo);
        this.updateClientsGui();
    }

    private void updateThreadsInfo(Intent intent) {
        ThreadInfo threadInfo = intent.getParcelableExtra(C.THREAD_INFO);
        int nArgs = intent.getIntExtra(C.THREAD_N_ARGUMENTS, 0);
        Argument[] arguments = new Argument[nArgs];
        for (int i = 0; i < nArgs; i++) {
            arguments[i] = (Argument) intent.getSerializableExtra(C.THREAD_ARGUMENTS + i);
        }
        clientsController.updateThreadInfoAndArguments(threadInfo, arguments);
        this.updateClientsGui();
    }

    // Gui
    private void updateClientsGui() {
        int[] threadIDs = clientsController.getAllThreadIDs();
        int newIndex = table.getChildCount();
        for (int threadID : threadIDs) {
            if (!threadToView.containsKey(threadID)) {
                String title = clientsController.getThreadTitle(threadID);
                TableRow row = new TableRow(this);
                Button start = new Button(this);
                start.setText("Start " + title);
                start.setOnClickListener(getThreadStarter(threadID));
                Button stop = new Button(this);
                stop.setText("Stop " + title);
                stop.setOnClickListener(getThreadStopper(threadID));
                row.addView(start);
                row.addView(stop);
                threadToView.put(threadID, newIndex);
                table.addView(row, newIndex);
                newIndex++;
            }
        }
    }

    private View.OnClickListener getThreadStarter(final int threadID) {
        return new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                clientsController.startThread(threadID);
            }
        };
    }

    private View.OnClickListener getThreadStopper(final int threadID) {
        return new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                clientsController.stopThread(threadID);
            }
        };
    }

    //Interface
    public void startServer(View view) {
        String serverName = "";
        if (!serverController.isBufferServerServiceRunning()) {
            serverName = serverController.startServerService();
        }
    }

    public void startClients(View view) {
        String clientsName = "";
        if (!clientsController.isClientsServiceRunning()) {
            clientsName = clientsController.startClientsService();
        }
    }

    public void stopServer(View view) {
        boolean result = false;
        if (serverController.isBufferServerServiceRunning()) {
            result = serverController.stopServerService();
        }
    }

    public void stopClients(View view) {
        boolean result = false;
        if (clientsController.isClientsServiceRunning()) {
            result = clientsController.stopClientsService();
        }
    }
}

