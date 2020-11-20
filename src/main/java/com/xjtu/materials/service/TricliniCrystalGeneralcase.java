package com.xjtu.materials.service;

import com.sun.jna.Library;
import com.sun.jna.Native;

public interface TricliniCrystalGeneralcase {
    public interface JNATestDll extends Library {

        JNATestDll instanceDll  = (JNATestDll)Native.loadLibrary("TriclinicCrystalgeneralcase",JNATestDll.class);
        int triclinic_crystal_generalcase(String file0);
    }
}
