package com.xjtu.materials.service;

import com.sun.jna.Library;
import com.sun.jna.Native;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;

public interface TricliniCrystalGeneralcase {
    public interface JNATestDll extends Library {

        JNATestDll instanceDll  = (JNATestDll)Native.loadLibrary("TriclinicCrystalgeneralcase",JNATestDll.class);
        int triclinic_crystal_generalcase(String file0);
    }

    public static void main(String[] args) {
        String file = "C:\\Users\\Xiongkai\\Desktop\\Al_Elastic_Constants.txt";
        String name = "Al Elastic Constants.txt";
        if (name == "Al Elastic Constants.txt") {
            ArrayList<String> arr = new ArrayList<String>();
            try {
                FileInputStream fis = new FileInputStream(file);
                InputStreamReader isr = new InputStreamReader(fis);
                BufferedReader reader = new BufferedReader(isr);
                String line = "";
                int i = 0;
                while ((line=reader.readLine()) != null) {
                    if (line.contains("Elastic Stiffness Constants")) {
                        i = 1;
                    }
                    if (i >= 1 && i<10) {
                        i++;
                        if (i>=4) {
                            arr.add(line);
                        }
                    }
                }
                System.out.println(arr);
            } catch (Exception ex)  {
                ex.printStackTrace();
            }
            ArrayList matrix = new ArrayList();
            for (String s: arr) {
                String[] nums = s.split(" ");
                float[] row = new float[6];
                int i = 0;
                for (String num: nums) {
                    float n = Float.parseFloat(num);
                    row[i] = n;
                    i++;
                }
                matrix.add(row);
            }
            System.out.println(matrix);
        }
    }
}
