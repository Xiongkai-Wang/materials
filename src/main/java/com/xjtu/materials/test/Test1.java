package com.xjtu.materials.test;

import com.xjtu.materials.pojo.ThemoParam;
import com.xjtu.materials.service.ThermoService;
import com.xjtu.materials.serviceImpl.ThermoServiceImpl;
import org.python.core.Py;
import org.python.core.PySystemState;
import org.python.util.PythonInterpreter;
import org.springframework.beans.factory.annotation.Autowired;

public class Test1 {


    private ThemoParam themoParam = new ThemoParam();


    public  void test() {
        String temp = "name.txt";
        ThermoService service = new ThermoServiceImpl(themoParam);
        service.run(temp);


    }
}
