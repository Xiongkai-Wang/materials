package com.xjtu.materials.serviceImpl;

import com.xjtu.materials.pojo.ThemoParam;
import com.xjtu.materials.service.ThermoService;

import java.io.*;

/**
 * @author:Scout
 * @data:2020/11/13
 * @description:
 */
public class ThermoServiceImpl implements ThermoService {

    private ThemoParam themoParam;
    private String PythonPro = "pythonpro/calculatefit.py";
    private String CppPro = "cpppro/ANVT_metal_advanced.exe";


    public ThermoServiceImpl(ThemoParam themoParam){
        this.themoParam = themoParam;
    }

    public ThermoServiceImpl(ThemoParam themoParam,String storedir){
        this.themoParam = themoParam;

    }


    @Override
    public void getData(String filename)  {
        String[] cmd = {CppPro,
                String.valueOf(themoParam.getFrames()), String.valueOf(themoParam.getSpecies()), String.valueOf(themoParam.getGibbs()), String.valueOf(themoParam.getNormaldos()), String.valueOf(themoParam.getN()), String.valueOf(themoParam.getFormlua()), String.valueOf(themoParam.getBare()),
                String.valueOf(themoParam.getTmin()), String.valueOf(themoParam.getTmax()), String.valueOf(themoParam.getStep()), String.valueOf(themoParam.getM()), String.valueOf(themoParam.getPoisson()), String.valueOf(themoParam.getB0()),
                themoParam.getMass(),themoParam.getNum(),filename};
        try {
            File file = new File("cpppro");
            Process process = Runtime.getRuntime().exec(cmd, null, file);
            InputStream is1 = process.getInputStream();
            InputStream is2 = process.getErrorStream();
            BufferThread bt1 = new BufferThread(is1);
            BufferThread bt2 = new BufferThread(is2);
            bt1.start();
            bt2.start();

            process.waitFor();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    @Override
    public  void fit(String filename,String storename) {
        filename = "cpppro/"+filename;
        String[] args1 = new String[] { "python", PythonPro, filename,storename};
        try{
            Process pro = Runtime.getRuntime().exec(args1);
            BufferThread bt = new BufferThread(pro.getInputStream());
            BufferThread bt2 = new BufferThread(pro.getErrorStream());
            bt.start();
            bt2.start();
            pro.waitFor();
        } catch (IOException  |InterruptedException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void run(String storename){
        String filename = "data.txt";
        getData(filename);
        fit(filename,storename);
    }



}

