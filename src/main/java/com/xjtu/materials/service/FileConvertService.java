package com.xjtu.materials.service;

import com.xjtu.materials.pojo.UpLoadMaterial;

import java.io.IOException;

public interface FileConvertService {
    boolean cell2cif(UpLoadMaterial material);

    boolean cif2cell(UpLoadMaterial material) throws IOException, Exception;
}
