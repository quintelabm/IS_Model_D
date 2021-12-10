#ifndef _IS_Model_H_
#define _IS_Model_H_

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

//each 10 space equals 1 cm (1 × 10^−2 m)
const int    Xspace   = 10.0; // 10 mm = 1 cm
const int    Yspace   = 10.0; // 10 mm = 1 cm
const int    Zspace   = 10.0; // 10 mm = 1 cm
const int    SPACE    = Xspace * Yspace * Zspace; // 1 cm^3
const int    IC_SPACE = 0.4*(Xspace * Yspace * Zspace); //initial condition
const int    source   = 100*pow(10,0);
const double SCALE    = pow(10,-3);
const double MOL      = 6.02*pow(10,23);
const int    buffer   = 2;

class IS_Model{

  private:

    double A [buffer][Xspace][Yspace][Zspace];
    double CH [buffer][Xspace][Yspace][Zspace];
    double CA [buffer][Xspace][Yspace][Zspace];
    double MR [buffer][Xspace][Yspace][Zspace];
    double MA [buffer][Xspace][Yspace][Zspace];
    double F [buffer][Xspace][Yspace][Zspace];
    double D [buffer][Xspace][Yspace][Zspace];

    int simCase;
    int days;
    int points;
    double deltaT;
    double iterPerDay;
    double deltaX, deltaY, deltaZ;

    int lnv;
    int bv;
    double tol;
    double source_mr;
    double migration_ma;
    double migration_f;

    double m0;
    double a0;
    double th0;
    double b0;
    double p0;
    double f0;
    double t_estrela;
    double b_estrela;
    double p_estrela;
    double f_estrela;
    double m_estrela;
    double MA_T, MA_L, MR_T;
    double Th, B, P;
    double F_T, F_L, A_T;
    double CH_T, CA_T, D_T;
    double d_a;
    double d_mr;
    double d_ma;
    double beta_A;
    double k_A;
    double k_F;
    double m_A;
    double m_Mr;
    double m_Ma;
    double ca0;
    double chinf;
    double cainf;
    double d_ch;
    double d_ca;
    double m_ch;
    double m_ca;
    double b_ch_ma;
    double b_ca_ma;
    double x_Mr;
    double x_Ma;
    double x_F;
    double alpha_Ma;
    double gamma_ma;
    double lambda_mr;
    double lambda_ma;
    double lambda_afmr;
    double lambda_afma;
    double b_th;
    double b_p;
    double b_pb;
    double b_pp;
    double ro_t;
    double ro_b;
    double ro_p;
    double ro_f;
    double alpha_t;
    double alpha_b;
    double alpha_p;
    double alpha_f;
    double alpha_mr;
    double alpha_d;
    double theta_d;
    double d_f;

    int saveFiles;
    char *dir;
    FILE* datamatlabA;
    FILE* datamatlabCH;
    FILE* datamatlabCA;
    FILE* datamatlabMr;
    FILE* datamatlabMa;
    FILE* datamatlabT;
    FILE* datamatlabB;
    FILE* datamatlabP;
    FILE* datamatlabF;
    FILE* datamatlabL;
    FILE* datamatlabD;

    std::string Header();
    std::string Footer(long int t);
    int checkFile(FILE* theFile);    
    int calcIntegral(double vec[][Xspace][Yspace][Zspace], double *V);
    int calcIntegral_lv(double vec[][Xspace][Yspace][Zspace], double *V);
    int calcIntegral_bv(double vec[][Xspace][Yspace][Zspace], double *V);
    void initialize();
    void update(double vec[][Xspace][Yspace][Zspace]);
    double laplacian(double vec[][Xspace][Yspace][Zspace], int x, int y, int z);
    double chemotaxis(double vec[][Xspace][Yspace][Zspace], int x, int y, int z, int ite);     
    int is_bvase(int x, int y, int z);
    int is_lnvase(int x, int y, int z);

  public:
    IS_Model();
    //IS_Model(int sCase, int sFile);
    //~IS_Model();
    void setSaveFiles(int sf);
    void setSimulationCase(int sc);
    int solve();

};

#endif
