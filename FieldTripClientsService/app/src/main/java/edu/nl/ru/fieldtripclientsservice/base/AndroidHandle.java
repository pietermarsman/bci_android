package edu.nl.ru.fieldtripclientsservice.base;

import android.content.Context;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

public interface AndroidHandle {

	FileInputStream openReadFile(String path) throws IOException;

	FileOutputStream openWriteFile(String path) throws IOException;

	void toast(String message);

	void toastLong(String message);

	void updateStatus(String status);

    /**
     * Useful for revering to the environment of the android device
     * @return The android context
     */
    public Context getContext();
}
