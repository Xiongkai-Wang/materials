package com.xjtu.materials.factory;

import com.xjtu.materials.pojo.ThemoParam;
import com.xjtu.materials.service.ThermoService;
import com.xjtu.materials.serviceImpl.ThermoServiceImpl;

/**
 * @author:Scout
 * @data:2020/11/14
 * @description:
 */
public class ThermoFactory {

    public static ThermoService getThermoService(ThemoParam themoParam) {
        return new ThermoServiceImpl(themoParam);
    }
}
