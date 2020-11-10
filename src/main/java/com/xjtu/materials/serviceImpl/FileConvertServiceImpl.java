package com.xjtu.materials.serviceImpl;

import com.xjtu.materials.service.FileConvertService;
import com.xjtu.materials.pojo.UpLoadMaterial;


import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;


public class FileConvertServiceImpl implements FileConvertService {
    @Override
    public boolean cell2cif(UpLoadMaterial material){
        return false;
    }
    @Override
    public boolean cif2cell(UpLoadMaterial material) throws Exception {
        String exe = "python";
        String command = "D:\\calculator_simple.py";
        String num1 = "1";
        String num2 = "2";
        String[] cmdArr = new String[] {exe, command, num1, num2};
        Process process = Runtime.getRuntime().exec(cmdArr);
        InputStream is = process.getInputStream();
        DataInputStream dis = new DataInputStream(is);
        String str = dis.readLine();
        process.waitFor();
        System.out.println(str);
        return  true;
    }
}
