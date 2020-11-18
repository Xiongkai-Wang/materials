package com.xjtu.materials.pojo;

import org.springframework.stereotype.Component;

/**
 * @author:Scout
 * @data:2020/11/13
 * @description:
 */
@Component
public class ThemoParam {
    int frames = 10;
    int species = 2;
    int Gibbs = 1;
    int N = 4;
    int normaldos = 1;
    int formlua = 1;
    int bare = 1;
    double Tmin = 5;
    double Tmax = 2000;
    double step = 40;
    double M = 669.566;
    double Poisson = 0.4;
    double B0 = 500;
    String mass = "192.22 92.906";
    String num = "3 1";
    public ThemoParam(){

    }

    public int getFrames() {
        return frames;
    }

    @Override
    public String toString() {
        return "ThemoParam{" +
                "frames=" + frames +
                ", species=" + species +
                ", Gibbs=" + Gibbs +
                ", N=" + N +
                ", normaldos=" + normaldos +
                ", formlua=" + formlua +
                ", bare=" + bare +
                ", Tmin=" + Tmin +
                ", Tmax=" + Tmax +
                ", step=" + step +
                ", M=" + M +
                ", Poisson=" + Poisson +
                ", B0=" + B0 +
                ", mass='" + mass + '\'' +
                ", num='" + num + '\'' +
                '}';
    }

    public void setFrames(int frames) {
        this.frames = frames;
    }

    public int getSpecies() {
        return species;
    }

    public void setSpecies(int species) {
        this.species = species;
    }

    public int getGibbs() {
        return Gibbs;
    }

    public void setGibbs(int gibbs) {
        Gibbs = gibbs;
    }

    public int getN() {
        return N;
    }

    public void setN(int n) {
        N = n;
    }

    public int getNormaldos() {
        return normaldos;
    }

    public void setNormaldos(int normaldos) {
        this.normaldos = normaldos;
    }

    public int getFormlua() {
        return formlua;
    }

    public void setFormlua(int formlua) {
        this.formlua = formlua;
    }

    public int getBare() {
        return bare;
    }

    public void setBare(int bare) {
        this.bare = bare;
    }

    public double getTmin() {
        return Tmin;
    }

    public void setTmin(double tmin) {
        Tmin = tmin;
    }

    public double getTmax() {
        return Tmax;
    }

    public void setTmax(double tmax) {
        Tmax = tmax;
    }

    public double getStep() {
        return step;
    }

    public void setStep(double step) {
        this.step = step;
    }

    public double getM() {
        return M;
    }

    public void setM(double m) {
        M = m;
    }

    public double getPoisson() {
        return Poisson;
    }

    public void setPoisson(double poisson) {
        Poisson = poisson;
    }

    public double getB0() {
        return B0;
    }

    public void setB0(double b0) {
        B0 = b0;
    }

    public String getMass() {
        return mass;
    }

    public void setMass(String mass) {
        this.mass = mass;
    }

    public String getNum() {
        return num;
    }

    public void setNum(String num) {
        this.num = num;
    }
}
