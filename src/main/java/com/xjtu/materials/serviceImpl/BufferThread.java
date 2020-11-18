package com.xjtu.materials.serviceImpl;


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * @author:Scout
 * @data:2020/11/17
 * @description:
 */
public class BufferThread extends Thread {
    BufferedReader bf;
    String str;
    InputStream is;
    public BufferThread(InputStream is){
        this.is = is;
        this.bf = new BufferedReader(new InputStreamReader(is));
    }
    @Override
    public void run() {
        super.run();
        try {
            while ((str = bf.readLine()) != null)
                System.out.println(str);
            is.close();
        }catch (IOException e){
            e.printStackTrace();
        }
    }
}
