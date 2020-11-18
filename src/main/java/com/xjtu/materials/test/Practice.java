package com.xjtu.materials.test;

import com.xjtu.materials.pojo.ThemoParam;

import java.lang.reflect.Field;

/**
 * @author:Scout
 * @data:2020/11/17
 * @description:
 */
public class Practice {
    public static void main(String[] args) {
        Field[] fields = new ThemoParam().getClass().getDeclaredFields();
        for(Field field:fields){
            String name =field.getName();
            System.out.printf("\'%s\':%s\n",name,name);
        }
    }
}
