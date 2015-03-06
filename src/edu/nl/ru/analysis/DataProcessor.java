package edu.nl.ru.analysis;

import nl.fcdonders.fieldtrip.bufferclient.BufferClientClock;
import nl.fcdonders.fieldtrip.bufferclient.BufferEvent;
import nl.fcdonders.fieldtrip.bufferclient.Header;
import nl.fcdonders.fieldtrip.bufferclient.SamplesEventsCount;
import org.apache.log4j.Logger;

import java.io.IOException;

/**
 * Created by Pieter on 23-2-2015.
 */
public class DataProcessor implements Runnable {

    String buffHost, endValue;
    NFEvent endType, predEventType, baselineEventType;
    int buffPort, baselineEnd, timeout_ms;
    double overlap;
    double[] lengthTrMs, lengthTrSample, stepMs, predFilt;
    boolean normLat, verbose;
    Object[] hdr;

    static Logger log = Logger.getLogger(DataProcessor.class);

    public DataProcessor() {
        log.setLevel(null);
        buffHost = "localhost";
        buffPort = 1972;
        hdr = new Object[]{};
        endType = NFEvent.STIMULUS_TEST;
        endValue = "end";
        verbose = false;
        predEventType = NFEvent.CLASSIFIER_PREDICTION;
        baselineEventType = NFEvent.STIMULUS_BASELINE;
        baselineEnd = 1000;
        lengthTrMs = new double[]{};
        lengthTrSample = new double[]{};
        overlap = .5;
        stepMs = new double[]{};
        predFilt = new double[]{};
        timeout_ms = 1000;
        normLat = true;
    }

    @Override
    public void run() {
        log.info("Running the dataprocessor");

        String[] args = new String[]{};

        String hostname = "localhost";
        int port = 1972;
        int timeout = 5000;

        if (args.length >= 1) {
            hostname = args[0];
        }
        if (args.length >= 2) {
            try {
                port = Integer.parseInt(args[1]);
            } catch (NumberFormatException e) {
                port = 0;
            }
            if (port <= 0) {
                log.error("Second parameter (" + args[1] + ") is not a valid port number.");
                System.exit(1);
            }
        }
        if (args.length >= 3) {
            try {
                port = Integer.parseInt(args[2]);
            } catch (NumberFormatException e) {
                timeout = 5000;
            }
        }

        BufferClientClock C = new BufferClientClock();

        Header hdr = null;
        while (hdr == null) {
            try {
                log.info("Connecting to " + hostname + ":" + port);
                C.connect(hostname, port);
                //C.setAutoReconnect(true);
                if (C.isConnected()) {
                    hdr = C.getHeader();
                }
            } catch (IOException e) {
                hdr = null;
            }
            if (hdr == null) {
                log.warn("Invalid Header... waiting");
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }
        //float[][] data = C.getFloatData(0,hdr.nSamples-1);
        log.info("#channels....: " + hdr.nChans);
        log.info("#samples.....: " + hdr.nSamples);
        log.info("#events......: " + hdr.nEvents);
        log.info("Sampling Freq: " + hdr.fSample);
        log.info("data type....: " + hdr.dataType);
        for (int n = 0; n < hdr.nChans; n++) {
            if (hdr.labels[n] != null) {
                log.debug("Ch. " + n + ": " + hdr.labels[n]);
            }
        }

        // Now do the echo-server
        int nEvents = hdr.nEvents;
        boolean endExpt = false;
        while (!endExpt) {
            SamplesEventsCount sec = null; // Block until there are new events
            try {
                sec = C.waitForEvents(nEvents, timeout);
            } catch (IOException e) {
                e.printStackTrace();
            }
            if (sec.nEvents > nEvents) {
                // get the new events
                BufferEvent[] evs = new BufferEvent[0];
                try {
                    evs = C.getEvents(nEvents, sec.nEvents - 1);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                nEvents = sec.nEvents;// update record of which events we've seen
                // filter for ones we want
                log.info("Got " + evs.length + " events");
                for (int ei = 0; ei < evs.length; ei++) {
                    BufferEvent evt = evs[ei];
                    String evttype = evt.getType().toString(); // N.B. to*S*tring, not upper case!
                    // only process if it's an event of a type we care about
                    // In our case, don't echo our own echo events....
                    if (!evttype.equals("echo")) {  // N.B. use equals, not == to compare string contents!
                        if (evttype.equals("exit")) { // check for a finish event
                            endExpt = true;
                        }
                        // Print the even to the console
                        log.info(ei + ") t:" + evt.getType().toString() + " v:" + evt.getValue().toString() + " s:" + evt.sample);
                        // Now create the echo event, with auto-completed sample number
                        // N.B. -1 for sample means auto-compute based on the real-time-clock
                        try {
                            C.putEvent(new BufferEvent("echo", evt.getValue().toString(), -1));
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                }
            } else { // timed out without new events
                log.debug("Timeout waiting for events");
            }
        }
        try {
            C.disconnect();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
