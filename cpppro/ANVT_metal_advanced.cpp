
//C++ code written for tabulated data//
//************************************//
//Code is wirtten by Bing Xiao//
// Emial: bxiao1@tulane.edu//
//Deparment of Physics and Eingineering Physics//
//Tulane University, Louisiana, 70118, USA//
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <string>
#define kb 1.3806505E-23                // Boltzman constant
#define H 1.60217646E-19                    // convert the eV to J
#define h 1.05457168E-34             // reduced Planck's constant
#define w 1.0E12                    // convert THz to Hz
#define vmin 0.0001                // lower bound for Debye integral
#define vmin_1 0.0001
#define max_in 50000                // Total number of intervals for Debye intergration
#define N_a 6.02E23                  //Avgavdlo constant
#define PI 3.14159265083
#define e_c 1.60217646E-19         // Unit C, elementary charge
#define weightt 1.660538921E-27     // Atomic mass to kg
#define nterms 4                    // number of fitting parameters in EOS
#define maxnpts 50                  // Maximum number of volumes, sufficient for QHA calculation
#define maxiterations 20    // Maximum iteration in the fitting loop
#define boundary 10.0 // boundary between dense T grid and corase one
/*------------------------------------------------------------------------------*/
   int param, iteration, nloops, n, cycle, nfree;
   int npts; // number of pairs in the fitting data
   double x[maxnpts], y[maxnpts], sigmay[maxnpts]; // data Volume~energy~uncertainty
   double z[maxnpts]; // This is an auxiliary quantity
   double weight[maxnpts];
   double yfit[maxnpts];
   double a[nterms], sigmaa[nterms], b[nterms];
   double beta[nterms], c[nterms];
   double finala[nterms], lastsigmaa[nterms];
   double alpha[nterms][nterms], arry[nterms][nterms];
   double aug[nterms][nterms*2];  // Matrix inversion
   double deriv[maxnpts][nterms]; // Derivatives
   double flambda; // Proporition of gradient search (0.001)
   double chisq;
   double chisql, fchisq, sy;
   char errorchoice;
   char filename3[128], answer[128];
   FILE *fp; // Input data
// Subroutine functions
//   void readdata(void), unweightedinput(void);
//   void weightedinput(void);
   void chisquare(void);
   void calcderivative(void), matrixinvert(void);
   void curvefit(int npoints);
//   void display(void);
   void uncertainties(void), jackknifedata(char *filename3, int k);
//   void print_data(void);
/*------------------------------------------------------------------------------*/
/* Subroutine for data interpolation*/
/* Lagrangian Interpolation Method */
double inter(double x2, double *xdata, double *ydata, int firstpoint, int numpoints)
{
       int i,j;
       double sum=0.0;
       double term;
       for(i=firstpoint;i<firstpoint+numpoints; i++){
                                               term=ydata[i];
                                               for(j=firstpoint; j<firstpoint+numpoints;j++){
                                                                 if(i!=j){
                                                                          term=term*(x2-xdata[j])/(xdata[i]-xdata[j]);
                                                                          }
                                                                 }
                                               sum=sum+term;
                                               }
                                               return sum;
}
/*Local Linear Interpolation */
//double inter1
/*------------------------------------------------------------------------------*/

void inttoString(int* key,int value){
    *key = value;
    printf("%d\n",value);
}

void doubletoString(double* key,double value){
    *key = value;
    printf("%lf\n",value);
}

void floattoString(float* key,float value){
    *key = value;
    printf("%f\n",value);
}
int main(int argc,char *argv[])
{
//   printf("/*###################################################################*/\n");
//   printf("/*#  Non-Linear Least Square Fitting Program For Free Energy Curve  #*/\n");
//   printf("/*#                        VERSION 1.01                             #*/\n");
//   printf("/*#    Equation of State Form: Birch-Murnaghan Third Order EOS      #*/\n");
//   printf("/*###################################################################*/\n");
//   printf("/*            THIS PROGRAM IS WRITTEN BY BING XIAO.                  */\n");
//   printf("/*             FIRST VERSION RELEASED: 18/03/2015                    */\n");
//   printf("/*              THE LATEST UPDATE: 10/09/2019.                       */\n");
//   printf("/*           CONTACT INFORMATION: b.xiao@ucl.ac.uk.                  */\n");
//   printf("/*             Department of Earth Sciences, UCL.                    */\n");
//   printf("/*-------------------------------------------------------------------*/\n");
//   printf("/*--------------------  Vibrational Part-----------------------------*/\n");
//   printf("/* This program evluates the thermodynamic properties of metallic    */\n");
//   printf("/* structures by three methods: Phonon density of states (1)         */\n");
//   printf("/* Gamma-point only phonon frequencies (2)  Debye Model (3)          */\n");
//   printf("/*                     Quantities Computed                           */\n");
//   printf("/*       Internal Energy(U), Zero Point Energy(ZPE), Entropy(S);     */\n");
//   printf("/*          Helmoltz Free Energy(A), Specific Heat (Cv)              */\n");
//   printf("/*----------------------Electron Subsystem---------------------------*/\n");
//   printf("/* Electron excitation-part is calculated using electron TDOS.       */\n");
//   printf("/*                     Quantities Computed                           */\n");
//   printf("/*       Internal Energy(Ue), Configurational Entropy(Se)            */\n");
//   printf("/*                   Electron Hemoltz Free Energy(Ae)                */\n");
//   printf("/*     Electron Specific Heat:Low Temperature Approximation(Cel)     */\n");
//   printf("/*     Electron Specific Heat:Exact Analytical Expression(Cel_f)     */\n");
//   printf("/*-------------------------------------------------------------------*/\n");
//   printf("/*----------------------Electric Resistivity-------------------------*/\n");
//   printf("/*     Electric resistivity is calculated by Bloch-Gruneisen model.  */\n");
//   printf("/*                  Electron-Electron Scattering(n=2)                */\n");
//   printf("/*            s-d Electron Scattering(Transition Metals)(n=3)        */\n");
//   printf("/*               Electron-Phonon Scattering(n=5)(Most Metals)        */\n");
//   printf("/*             Thermal Excited Electrons Near Fermi Level            */\n");
//   printf("/* Reference: R.G. Goetsch,et.al.,Physical Review B,85,054517(2012)  */\n");
//   printf("/*------------------ Electron Thermal Conductivity-------------------*/\n");
//   printf("/*   Electronic thermal conductivity is estimated empirically by:    */\n");
//   printf("/*                         Wiedemann-Franz Law                       */\n");
//   printf("/*   The total electric conductivity is used in the calculation.     */\n");
//   printf("/*-------------------------------------------------------------------*/\n");
//   printf("/*----------------------IMPORTANT REMINDER---------------------------*/\n");
//   printf("/* You must prepare the following files in order to run the program. */\n");
//   printf(" phonon_Vx.dat (phonon frequencies or DOS are stored for the volume V).\n");
//   printf(" EV_dft.dat (DFT static calculations for energy versus volume curve). \n");
//   printf(" tdos_Vx.dat (DFT computed total electron density of states at Vx).   \n");
//   printf("/*----------------------FOR S-D SCATTERING---------------------------*/\n");
//   printf("/*    You also need provide the partial DOS for s+d states.          */\n");
//   printf("/*                        sddos_Vx.dat                               */\n");
//   printf("/*-------------------------------------------------------------------*/\n");
//   printf("/*###################################################################*/\n");
//   printf("/*#     The follow paper should be cited in your work, thanks!      #*/\n");
//   printf("/*#     Yingchun Ding, Bing Xiao, RSC Advances, 5, 18391 (2015).    #*/\n");
//   printf("/*###################################################################*/\n");
   int i;
   int L;
   int P_1, P_2, P_3, P_4;
   int rows, rown, rowns, rownp;
   int bare;
   double T;                      // temperature
   double A_lat;                    //Helmoltz free energy
   double entropy_e;                  // Entropy of electrons at finite temperature
   double entropy_l;
   double entropy_t; // S_e+S_l
   double zpe;                   // Zero point energy
   double Ezpe;
   double E_1, E_2;
   double A_e;                  // electron free energy
   double A_t;   // The addition of free energies of lattice and electrons
   float *E_x, *V_x; // DFT computed energy and volume
   float B0;   // Bulk modulus at equilibrium
   float M; // Molecule weight
   float m_f; // Mass per formula unit
   float inv_m;
   float V; // Cell volume
   float V_se; // Seitz radius
   float Rou_5, R_5; // Electric conductivity: Phonon scattering
   float Rou_3, R_3; // Electric resistivity: s-d scattering
   float Rou_2, R_2;  //electron-electron correlation
   double U_e;       // Total internal energy of valenece electrons
   double U_l;       // Total internal energy of lattice vibration
   double U_t; // U_e+U_l+zpe+E_x[i]
   double Cv;     // Specific heat of lattice vib
   double C;     // ratio of T and Debye temperature
   double Cvt;  // total specific heat Cvt = Cve+Cel
   double density; // Density
   double N_e; // Total number of conducting eletrons
   double N_sd; // Total numer of conducting sd electrons
   double p_1, p_2, p_3;
   float Tmax, Tmin;
   float step;
   double vmax;
   double vmax_1;
   double Debye_1, Debye_2;  // Intermediate parameters
   float Poisson; // Poisson's ratio
   float Debye;       //Debye temperature
   float mass_A, mass_B, mass_C, mass_D, mass_E, mass_F;
   float num_A, num_B, num_C, num_D, num_E, num_F;
   double Lwf; // Wiedamman-Franz constant
   double lamada; // Thermal conductivity
   int species;
   char file1[128],file2[128],file3[128],file4[128];
   char files1[128],files2[128],files3[128],files4[128];
   char files5[128],files6[128],files7[128];
   char files8[128], files9[128], files10[128];
   char filename[128];
   char filename1[128], filename2[128]; // output11, and output12;
   int SMAX, PMAX;
   int tmp, tmp_1;
   int N;         // Total number of atoms per cell
   float N_f;      // Total number of atom per formula
   int formula; //	number of molecules in the cell
   int j,l;
   int Gibbs;
   int normaldos;
   int frames; // Total number of volumes
   int u;  // Simpson's rule
   int rowtemp; // total number of different temperatures;
   int rowtotal; // total number of rows in combined file
   float dummy1, dummy2, dummy3, dummy4;
   int dummy;
   double x1,y1,z1;
   double *xdata,*ydata; // Frequency and dos
   int *gdata; // Degeneracy of phonon modes at Gamma point
   double *dos,*edata;
   double *pdos_e,*pdos; // s-d density of states
   float f(double, double, double);  // Fermi-Dirac distribution function
   float devf(double, double, double); // Derivative of f
   float g(double);
   float f1(double);
   float f2(double);
   float Cut_1(double);
   float Cut_2(double);
   float Cut_3(double);
   float De(double);
   float simpson(int no, float min, float max);
   float simpson_1(int no, float min, float max);
   float simpson_2(int no, float min, float max);
   float simpson_3(int no, float min, float max);
// Output files
   FILE *input; // phonon_Vx.dat
   FILE *input1; // tods_Vx.dat
   FILE *input2;  // sddos_Vx.dat
   FILE *input3;  // EV_dft.dat
   FILE *output;
   FILE *output1;
   FILE *output2;
   FILE *output3;
   FILE *output4;
   FILE *output5;
   FILE *output6;
   FILE *output7;
   FILE *output8;
   FILE *output9;
   FILE *output10;
   FILE *output11;  // Data combined file
   FILE *output12; // EV data at each temperature
   FILE *output13;
// Input parameters

   printf("/*------------------Please Provide Required Parameters---------------*/\n");
   printf("Total number of volume-files in the calculations:(MUST BE AN INTEGER)\n");
   frames = atoi(argv[1]);
   printf("frames=%d\n",frames);
   //scanf("%d",&frames);
   printf("Number of different species in structure:(MUST BE AN INTEGER)\n");
   species = atoi(argv[2]);
   printf("species=%d\n",species);
  // scanf("%d",&species);
   printf("Choose the method to calculate the thermodynamic quantities of lattice:\n");
   printf("Phonon density of states(1);Gamma-point only frequency(2);Debye model(3)\n");
   printf("Enter your option below:(MUST BE AN INTEGER)\n");
   Gibbs = atoi(argv[3]);
   printf("Gibbs=%d\n",Gibbs);
  // scanf("%d",&Gibbs);
   printf("Enter number of atoms per cell:(MUST BE AN INTEGER)\n");
   N = atoi(argv[4]);
   printf("N=4\n",N);
  // scanf("%d",&N);
   printf("The input phonon dos is normalized to: 3N(1) or 1(2)?\n");
   printf("Enter the nromalization option:(1 or 2)\n");
   inttoString(&normaldos,atoi(argv[5]));
  // scanf("%d",&normaldos);
   printf("Enter total number of formula units in the cell:(MUST BE AN INTEGER)\n");
   inttoString(&formula,atoi(argv[6]));
  // scanf("%d", &formula);
   printf("Output options for lattice and electronic free energies:\n");
   printf("1:Lattice or electron energy plus E_DFT\n");
   printf("2:Lattice or electron contribution only\n");
   printf("Enter your option below:\n");
   inttoString(&bare,atoi(argv[7]));
  // scanf("%d",&bare);
   printf("#####################################################################\n");
   printf("#Number of isothermal E-V curves is defined by the following inputs.#\n");
   printf("#          User should think carefully before typing!               #\n");
   printf("#All other temperature dependent thermal physical properties are    #\n");
   printf("#anticipated to change moderately (or slightly).                    #\n");
   printf("#        Mode one: Dense Temperature Grid, T_step <= 5K             #\n");
   printf("#             Mode two: Corase T-grid, T_step > 5K                  #\n");
   printf("#####################################################################\n");
   printf("Enter temperatures(K):Starting, Final and Step\n");
   printf("Format: 10.0 1000.0 1.0\n");
   Tmin = atof(argv[8]);
   Tmax = atof(argv[9]);
   step = atof(argv[10]);
   printf("%f %f %f\n",Tmin,Tmax,step);
 //  scanf("%f %f %f",&Tmin,&Tmax,&step);
   printf("Molecule weight(g/mol.cell)\n");
   floattoString(&M,atof(argv[11]));
 //  scanf("%f",&M);
   printf("The following parameters are defined for equilibrium geometry:\n");
   printf("Enter poisson's ratio:\n");
   floattoString(&Poisson,atof(argv[12]));
 //  scanf("%f", &Poisson);
   printf("Enter bulk modulus(GPa):\n");
   floattoString(&B0,atof(argv[13]));
 //  scanf("%f", &B0);
   printf("Enter the parameter for each element in the cell.\n");
   /*
   if(species==1){
                  printf("Element A:\n");
                  printf("Enter the atomic weight:(Float number)\n");
   //                scanf("%f",&mass_A);
                  printf("Number of atoms per cell:(MUST BE AN INTEGER)\n");
                  printf("Example: 2 \n");
     //             scanf("%d",&num_A);
                  inv_m = (float)(1.0*num_A*(1.0/weightt)*(N_a/N)/mass_A);
                  }
   else if(species==2){
                       printf("Element A and B:\n");
                       printf("Enter the atomic weight:(24.305 15.999)\n");
                       mass_A = 192.22;
                       mass_B = 92.906;
                       printf("%f %f\n",mass_A,mass_B);
                       //scanf("%f %f",&mass_A,&mass_B);
                       printf("Number of atoms per cell:(MUST BE TWO INTEGERS)\n");
                       printf("Example: 2 4\n");
                       //scanf("%d %d",&num_A,&num_B);
                       floattoString(&num_A,3);
                       floattoString(&num_B,1);
                       inv_m = (float)((N_a/N)*(1.0/weightt)*(1.0*num_A/mass_A+1.0*num_B/mass_B));
                       }
   else if(species==3){
                       printf("Element A, B, C:\n");
                       printf("Enter the atomic weight:(24.305 15.999 12.000)\n");
                       scanf("%f %f %f",&mass_A,&mass_B,&mass_C);
                       printf("Number of atoms per cell:(MUST BE THREE INTEGERS)\n");
                       printf("Example: 2 4 5\n");
                       scanf("%d %d %d",&num_A,&num_B,&num_C);
                       inv_m = (float)((N_a/N)*(1.0/weightt)*(1.0*num_A/mass_A+1.0*num_B/mass_B+1.0*num_C/mass_C));
                       }
   else if(species==4){
                       printf("Element A, B, C, D:\n");
                       printf("Enter the atomic weight:(24.305 15.999 12.000 56.065)\n");
                       scanf("%f %f %f %f",&mass_A,&mass_B,&mass_C,&mass_D);
                       printf("Number of atoms per cell:(MUST BE FOUR INTEGERS)\n");
                       printf("Example: 2 4 5 8\n");
                       scanf("%d %d %d %d",&num_A,&num_B,&num_C,&num_D);
                       inv_m = (float)((N_a/N)*(1.0/weightt)*(1.0*num_A/mass_A+1.0*num_B/mass_B+1.0*num_C/mass_C+1.0*num_D/mass_D));
                       }
   else if(species==5){
                       printf("Element A, B, C, D, E:\n");
                       printf("Enter the atomic weight:(1.000 24.305 15.999 12.000 56.065)\n");
                       scanf("%f %f %f %f %f",&mass_A,&mass_B,&mass_C,&mass_D,&mass_E);
                       printf("Number of atoms per cell:(MUST BE FIVE INTEGERS)\n");
                       printf("Example: 2 4 10 3 5\n");
                       scanf("%d %d %d %d %d",&num_A,&num_B,&num_C,&num_D,&num_E);
                       inv_m = (float)((N_a/N)*(1.0/weightt)*(1.0*num_A/mass_A+1.0*num_B/mass_B+1.0*num_C/mass_C+1.0*num_D/mass_D+1.0*num_E/mass_E));
                       }
                       */
   printf("please input %d numbers",species);
   // ∂¡»°mess;
   char *p;
   p = strtok(argv[14], " ");
   mass_A = atof(p);
   p = strtok(NULL, " ");
   mass_B = atof(p);
   p = strtok(NULL, " ");
   mass_C = atof(p);
   p = strtok(NULL, " ");
   mass_D = atof(p);
   p = strtok(NULL, " ");
   mass_E = atof(p);
   printf("mass : %f,%f,%f,%f,%f",mass_A,mass_B,mass_C,mass_D,mass_E);

   p = strtok(argv[15], " ");
   num_A = atof(p);
   p = strtok(NULL, " ");
   num_B = atof(p);
   p = strtok(NULL, " ");
   num_C = atof(p);
   p = strtok(NULL, " ");
   num_D = atof(p);
   p = strtok(NULL, " ");
   num_E = atof(p);
   printf("num : %f,%f,%f,%f,%f",num_A,num_B,num_C,num_D,num_E);

   printf("/*--------------------END OF INPUTS----------------------------------*/\n");
   if(Gibbs==1||Gibbs==2){
                          printf("The phonon frequencies must be angular frequency(THz).\n");
                          printf("Please Check Your Numbers Carefully.\n");
                          }
// Read the DFT static energy and volume data
/*------------------------------------------------------------------------------*/
/* determine the structure of EV_dft file */
   sprintf(file1,"EV_dft.dat");
   input3=fopen(file1,"r");
   if(input3==NULL){
                      perror("Error while opening the file.\n");
                      printf("%s does not exist.\n",file1);
                      getchar();
                      exit(EXIT_FAILURE);
                      }
      else {
           printf("File is loaded\n");
           }
   P_1=0;
   while(!feof(input3)){
                       fscanf(input3,"%f %f",&dummy1, &dummy2);
                       P_1++;
                       }
   rows=P_1-1;
   fclose(input3);
//   printf("Total number of rows in file '%s' is %d\n",file1,rows);
   if(rows!=frames) {
                     perror("INCONSISTENCY IS DETECTED, PLEASE CHECK YOUR DATA.\n");
                     printf("Total number of volumes is %d\n",frames);
                     printf("Total number of volumes in %s is %d\n",file1,rows);
                     getchar();
                     exit(EXIT_FAILURE);
                     }
// Allocate the memory to the EV curve
   E_x=(float*)malloc(rows*sizeof(float));
   V_x=(float*)malloc(rows*sizeof(float));
   if(E_x==NULL||V_x==NULL){
                            printf("Error-could not allocate an array.\n");
                            exit(EXIT_FAILURE);
                            }
   input3=fopen(file1,"r");
   for(i=0;i<rows;i++){
                       fscanf(input3,"%f %f",&(V_x[i]),&(E_x[i]));
                       }
   fclose(input3);
/*------------------------------------------------------------------------------*/
   sprintf(files9,"V_ZPE.dat");
   sprintf(filename1,"Data_combined.dat");
   output8=fopen(files9,"w");
   output11=fopen(filename1,"w");
   fprintf(output8,"V\tE(V)(eV)\tZPE(eV/cell)\n");
   rowtotal=0;
   for(i=1;i<=frames;i++){
                          if(Gibbs==1){
                                       sprintf(file2,"phonon_V%d.dat",i);
                                       input=fopen(file2,"r");
                                       if(input==NULL){
                                                       perror("Error while opening the file.\n");
                                                       printf("%s does not exist.\n",file2);
                                                       getchar();
                                                       exit(EXIT_FAILURE);
                                                       }
                                       P_2=0;
                                       while(!feof(input)){
                                                           fscanf(input,"%f %f",&dummy1, &dummy2);
                                                           P_2++;
                                                           }
                                       rown=P_2-1;
                                       fclose(input);
                                       xdata=(double*)malloc(rown*sizeof(double));
                                       ydata=(double*)malloc(rown*sizeof(double));
                                       if(xdata==NULL||ydata==NULL){
                                                                    printf("Error-could not allocate an array.\n");
                                                                    exit(EXIT_FAILURE);
                                                                    }
                                       input=fopen(file2,"r");
                                       for (j=0;j<rown;j++){
                                                            fscanf(input,"%lf %lf",&(xdata[j]),&(ydata[j]));
                                                            }
                                       fclose(input);
                                       }
                            else if(Gibbs==2){
                                              sprintf(file2,"phonon_V%d.dat",i);
                                              input=fopen(file2,"r");
                                              if(input==NULL){
                                                              perror("Error while opening the file.\n");
                                                              printf("%s does not exist.\n",file2);
                                                              getchar();
                                                              exit(EXIT_FAILURE);
                                                              }
                                              P_2=0;
                                              while(!feof(input)){
                                                                  fscanf(input,"%d %f",&dummy, &dummy2);
                                                                  P_2++;
                                                                  }
                                              rown=P_2-1;
                                              fclose(input);
                                              gdata=(int*)malloc(rown*sizeof(int));
                                              ydata=(double*)malloc(rown*sizeof(double));
                                              if(gdata==NULL||ydata==NULL){
                                                                           printf("Error-could not allocate an array.\n");
                                                                           exit(EXIT_FAILURE);
                                                                           }
                                              input=fopen(file2,"r");
                                              for (j=0;j<rown;j++){
                                                                   fscanf(input,"%d %lf",&(gdata[j]),&(ydata[j]));
                                                                   }
                                              fclose(input);
                                              }
// Read electronic dos and partial dos
                  sprintf(file3,"tdos_V%d.dat",i);
                  sprintf(file4,"sddos_V%d.dat",i);
                  input1=fopen(file3,"r");
                  input2=fopen(file4,"r");
                  if(input1==NULL){
                                  perror("Error while opening the file.\n");
                                  printf("Error: %s\n", strerror(errno));
                                  printf("%s does not exist.\n",file3);
                                  getchar();
                                  exit(EXIT_FAILURE);
                                  }
                  if(input2==NULL){
                                  perror("Error while opening the file.\n");
                                  printf("%s does not exist.\n",file4);
                                  getchar();
                                  exit(EXIT_FAILURE);
                                  }
                  P_3=0;
                  while(!feof(input1)){
                                      fscanf(input1,"%f %f",&dummy1, &dummy2);
                                      P_3++;
                                      }
                  rowns=P_3-1;
                  fclose(input1);
                  P_4=0;
                  while(!feof(input2)){
                                      fscanf(input2,"%f %f",&dummy3, &dummy4);
                                      P_4++;
                                      }
                  rownp=P_4-1;
                  fclose(input2);
// Allocate memory to the arrays
                  edata=(double*)malloc(rowns*sizeof(double));
                  dos=(double*)malloc(rowns*sizeof(double));
                  pdos_e=(double*)malloc(rownp*sizeof(double));
                  pdos=(double*)malloc(rownp*sizeof(double));
                  if(edata==NULL||dos==NULL||pdos_e==NULL||pdos==NULL){
                                                                       printf("Error-could not allocate an array.\n");
                                                                       exit(EXIT_FAILURE);
                                                                       }
                  input1=fopen(file3,"r");
                  for(j=0;j<rowns;j++){
                                       fscanf(input1,"%lf %lf",&(edata[j]),&(dos[j]));
                                      }
                  fclose(input1);
                  input2=fopen(file4,"r");
                  for(j=0;j<rownp;j++){
                                       fscanf(input2,"%lf %lf",&(pdos_e[j]),&(pdos[j]));
                                       }
                  fclose(input2);
// ALL DATA ARE IMPORTED
/*------------------------------------------------------------------------------*/
                  sprintf(files1,"ANVT_lattice_V%d.dat",i);
                  sprintf(files2,"ANVT_electron_V%d.dat",i);
                  sprintf(files3,"Resistivity_normal_V%d.dat",i);
                  sprintf(files4,"Resistivity_raw_V%d.dat",i);
                  sprintf(files5,"ANVT_le_V%d.dat",i);
                  sprintf(files8,"Thermal_el_V%d.dat",i);
                  sprintf(files10,"Thermal_condt_V%d.dat",i);
                  sprintf(filename,"Electron_Cv%d.dat",i);
                  output=fopen(files1,"w");
                  output2=fopen(files2,"w");
                  output3=fopen(files3,"w");
                  output5=fopen(files4,"w");
                  output6=fopen(files5,"w");
                  output7=fopen(files8,"w");
                  output9=fopen(files10,"w");
                  output10=fopen(filename,"w");
                  printf("T(K)\tA_l(eV/cell)\tA_e(eV/cell)\tS_e(eV/K.cell)\tS_l(eV/cell)\tCv(eV/cell)\tCv(J/K.cell)\tRows_e\n");
                  fprintf(output,"T(K)\tA_lat(eV/cell)\tS_lat(eV/K.cell)\tU_lat(eV/cell)\tCv(eV/K.cell)\n");
                  fprintf(output2,"T(K)\tA_el(eV/cell)\tS_el(eV/K.cell)\tU_el(eV/cell)\n");
                  fprintf(output3,"T(K)\tN_ce(e/cell)\tN_sd(e/cell)\tR_pe(Om.m)\tR_sd(Om.m)\tR_ee(Om.m)\n");
                  fprintf(output5,"T(K)\tN_ce(e/cell)\tN_sd(e/cell)\tR_pe(Om.m)\tR_sd(Om.m)\tR_ee(Om.m)\n");
                  fprintf(output6,"T(K)\tA(NVT)(eV/cell)\tU(NVT)(eV/cell)\tS(NVT)(eV/cell)\tCv(T)(eV/cell)\tCv(J/K.mol)\n");
                  fprintf(output7,"T(K)\tA_lat(eV/cell)\tA_e(eV/cell)\tU_lat(eV/cell)\tU_e(eV/cell)\tS_lat(eV/cell)\tS_e(eV/cell)\tRows\n");
                  fprintf(output9,"T(K)\tR_pe(Om.m)\tlamada(W/K.m)\n");
                  fprintf(output10,"T(K)\tCv(J/mol.K.cell)\tCv_full(J/mol.K.cell)\tDOS(states/eV)\n");
                  rowtemp=0;
                  for(T=Tmin;T<=Tmax;T=T+step){
// Resort the electron density of states files according to temperature
                                               rowtemp++;  // total number of different temperatures
                                               rowtotal++;  // total number of rows in combined file
                                               sprintf(files6,"tdos_new.dat");
                                               output1=fopen(files6,"w");
                                               for(j=0;j<rowns;j++){
                                                                   if(edata[j]>-10.0*(1.38/1.6)*pow(10.0,-4.0)*T&&edata[j]<10.0*(1.38/1.6)*pow(10.0,-4.0)*T){
                                                                                                           fprintf(output1,"%lf\t%lf\n",edata[j],dos[j]);
                                                                                                                                                     }
                                                                  }
                                               fclose(output1);
                                               sprintf(files7,"pdos_new.dat");
                                               output4=fopen(files7,"w");
                                               for(j=0;j<rownp;j++){
                                                                   if(pdos_e[j]>-10.0*(1.38/1.6)*pow(10.0,-4.0)*T&&pdos_e[j]<10.0*(1.38/1.6)*pow(10.0,-4.0)*T){
                                                                                                           fprintf(output4,"%lf\t%lf\n",pdos_e[j],pdos[j]);
                                                                                                                                                      }
                                                                   }
                                               fclose(output4);
/*------------------------------------------------------------------------------*/
                                               /*Properties of electron subsystem*/
                                               int SMAX, PMAX;
                                               int tmp, tmp_1;
                                               double *sdata,*sdos;
                                               double *sd,*sdo;
                                               double R;      // A constant, covert eV to J
                                               double F_1, F_11;
                                               double F_2, F_22;
                                               double F_3, F_4, F_5, F_6, F_7, F_8; // Variables
                                               float degen;
                                               double g1,g2,g11;
                                               double V_f;
                                               N_f=(float)(N/formula);
                                               m_f=(float)(M/formula)*weightt;
                                               V_f=V_x[i-1]/(float)formula;
                                               g1=3.0/(f1(Poisson)+f2(Poisson));
                                               g2=6.0*PI*PI*sqrt(V_f*pow(10.0,-30.0))*N_f;
                                               g11=pow(g1,1.0/3.0);
                                               Debye_1=(h/kb)*pow(g2,1.0/3.0)*g11*sqrt(B0*pow(10.0,9.0)/m_f);
                                               output1=fopen(files6,"r");
                                               tmp=0;
                                               while(!feof(output1)){
                                                                    fscanf(output1,"%f %f",&dummy3, &dummy4);
                                                                    tmp++;
                                                                    }
                                               SMAX=tmp-1;
                                               fclose(output1);
                                               sdata=(double*)malloc(SMAX*sizeof(double));
                                               sdos=(double*)malloc(SMAX*sizeof(double));
                                               output1=fopen(files6,"r");
                                               for(j=0;j<SMAX;j++){
                                                                   fscanf(output1,"%lf %lf",&(sdata[j]),&(sdos[j]));
                                                                   }
                                               fclose(output1);
                                               //-------------------------------
                                               output4=fopen(files7,"r");
                                               tmp_1=0;
                                               while(!feof(output4)){
                                                                     fscanf(output4,"%f %f",&dummy1, &dummy2);
                                                                     tmp_1++;
                                                                     }
                                               PMAX=tmp_1-1;
                                               fclose(output4);
                                               sd=(double*)malloc(PMAX*sizeof(double));
                                               sdo=(double*)malloc(PMAX*sizeof(double));
                                               output4=fopen(files7,"r");
                                               for(j=0;j<PMAX;j++){
                                                                   fscanf(output4,"%lf %lf",&(sd[j]),&(sdo[j]));
                                                                   }
                                               fclose(output4);
                                               //Conduction electrons
                                               R=(1.6021764/1.3806505)*pow(10.0,4.0);
                                               N_e=0.0;
                                               for(l=1;l<SMAX;l++){
                                                                   N_e+=f(sdata[l],R,T)*sdos[l]*fabs(sdata[l]-sdata[l-1]);
                                                                   }
                                               N_sd=0.0;
                                               for(l=1;l<PMAX;l++){
                                                                   N_sd+=f(sd[l],R,T)*sdo[l]*fabs(sd[l]-sd[l-1]);
                                                                   }
                                              //Electronic entropy Se(V,T)
                                              entropy_e=0.00;
                                              for(l=1;l<SMAX;l++){
                                                                  F_1=1.0-f(sdata[l],R,T);
                                                                  F_2=f(sdata[l],R,T);
                                                                  entropy_e+=(-(1.3806505/1.6)*pow(10.0,-4)*((F_1)*log(F_1)+F_2*log(F_2))*sdos[l]*fabs(sdata[l]-sdata[l-1]));
                                                                  }
                                             //Electron internal energy
                                             E_1=0.0;
                                             E_2=0.0;
                                             for(l=1;l<rowns;l++){
                                                                 F_11=1.0-f(edata[l],R,T);
                                                                 F_22=f(edata[l],R,T);
                                                                 E_1+=dos[l]*fabs(edata[l]-edata[l-1])*edata[l]*F_22;
                                                                 if(edata[l]<=0.0){
                                                                                   E_2+=dos[l]*fabs(edata[l]-edata[l-1])*edata[l];
                                                                                  }
                                                                 }
                                             A_e=(E_1-E_2)-T*entropy_e;
                                             U_e=E_1-E_2;
                                             //Electric resistivity
                                             V_se=pow((V_x[i-1]/N),0.33333);   // Seitz radius
                                             F_5=pow(N_e/N,0.66667);
                                             F_7=pow(N_sd/N,0.66667);
                                             F_6=4.77928*pow(10.0,-7);
                                             F_4=(1.0*inv_m)/((V_se)*(F_5)*Debye_1);
                                             F_8=(1.0*inv_m)/((V_se)*(F_7)*Debye_1);
                                             vmax_1=Debye_1/T;
                                             Debye_2=T/Debye_1;
                                             Rou_5=4.0*(F_6)*(F_4)*pow(Debye_2,5)*simpson_1(max_in, vmin_1, vmax_1);
                                             Rou_3=4.0*(F_6)*(F_8)*pow(Debye_2,3)*simpson_2(max_in, vmin_1, vmax_1);
                                             Rou_2=4.0*(F_6)*(F_4)*pow(Debye_2,2)*simpson_3(max_in, vmin_1, vmax_1);
                                             R_5=Rou_5/(0.9464635*(F_6)*(F_4));
                                             R_3=Rou_3/(0.9464635*(F_6)*(F_8));
                                             R_2=Rou_2/(0.9464635*(F_6)*(F_4));
                                             // Thermal conductivity
                                             Lwf=2.45*pow(10.0,-8); // W.Om/K^2
                                             lamada=(1.0/Rou_5)*Lwf*T;
                                             if(bare==1){
                                             	        fprintf(output2,"%f\t%lf\t%le\t%lf\n",T,A_e+E_x[i-1],entropy_e,U_e+E_x[i-1]);
                                                        }
                                             else if(bare==2){
                                             	        fprintf(output2,"%f\t%lf\t%le\t%lf\n",T,A_e,entropy_e,U_e);
                                                        }
                                             fprintf(output3,"%f\t%lf\t%lf\t%le\t%le\t%le\n",T,N_e,N_sd,R_5,R_3,R_2);
                                             fprintf(output5,"%f\t%lf\t%lf\t%le\t%le\t%le\n",T,N_e,N_sd,Rou_5,Rou_3,Rou_2);
                                             fprintf(output9,"%f\t%le\t%lf\n",T,Rou_5,lamada);
                                             //---------2016-02-12----------------------------
                                             // Specific heat of electron (thermal excitations)
                                             // Consider the low temperature expression first
                                             // Find the density of states at Fermi level
                                             double dos_f; // DOS at Fermi level
                                             double Cel; // electron specific heat caculated from low temperature expression
                                             double Cel_f;
                                             for(l=1;l<SMAX;l++){
                                                                 if(sdata[l]>0.000001&&sdata[l-1]<-0.000001){
                                                                                                            dos_f=0.5*(sdos[l]+sdos[l-1]);
                                                                                                            }
                                                                 }
                                             Cel=(2.0*PI*PI/3.0)*N_a*(dos_f/H)*kb*kb*T; //The unit for dos is states/eV.cell
                                             // The full expression for electron specific heat
                                             Cel_f=0.000000;
                                             for(l=1;l<rowns;l++){
                                                                 Cel_f += (edata[l]*R/T)*(edata[l]*R/T)*(dos[l]/H)*(-1.0)*devf(edata[l],R,T)*(2.0*kb*kb*T*N_a)*fabs(edata[l]-edata[l-1]);
                                                                 }
                                             fprintf(output10,"%f\t%le\t%le\t%lf\n",T,Cel,Cel_f,dos_f);
//--------------------------------------------------------------------------------------------------------------------------
//------------------------------Lattice Dynamics--------------------------------
// Calculate thermodynamical properties from Phonon density of states
                                             if(Gibbs==1){
                                                          A_lat=0.0;
                                                          zpe=0.0;
                                                          entropy_l=0.0;
                                                          U_l=0.0;
                                                          Cv=0.0;
                                                          for(j=1;j<rown;j++){
                                                                              if(xdata[j]>=0.0){
                                                                                                if(normaldos==2){
                                                                                                                 p_1=1-exp(-h*xdata[j]*w/(kb*T));
                                                                                                                 p_2=1.0/(exp(h*w*xdata[j]/(kb*T))-1);
                                                                                                                 p_3=1.0/(exp(-h*w*xdata[j]/(kb*T))-1);
                                                                                                                 zpe += (3.0*N/2.0)*w*h*(1.0/H)*(ydata[j]*xdata[j]*(xdata[j]-xdata[j-1]));
                                                                                                                 A_lat += (3.0*N*kb*T/H)*ydata[j]*log(p_1)*(xdata[j]-xdata[j-1]);
                                                                                                                 entropy_l += (3.0*N/H)*((-kb*ydata[j]*log(p_1)*(xdata[j]-xdata[j-1]))+(h*w*xdata[j]*ydata[j]/T)*p_2*(xdata[j]-xdata[j-1]));
                                                                                                                 U_l += (3.0*N/H)*(h*w*xdata[j]*ydata[j])*(p_2)*(xdata[j]-xdata[j-1]);
                                                                                                                 Cv += (3.0*N/H)*kb*ydata[j]*(xdata[j]*h*w/(kb*T))*(xdata[j]*h*w/(kb*T))*exp(-h*w*xdata[j]/(kb*T))*(pow(p_3,2.0)*(xdata[j]-xdata[j-1]));
                                                                                                                 }
                                                                                                else if(normaldos==1){
                                                                                                                      p_1=1-exp(-h*xdata[j]*w/(kb*T));
                                                                                                                      p_2=1.0/(exp(h*w*xdata[j]/(kb*T))-1);
                                                                                                                      p_3=1.0/(exp(-h*w*xdata[j]/(kb*T))-1);
                                                                                                                      zpe += (1.0/2.0)*w*h*(1.0/H)*(ydata[j]*xdata[j]*(xdata[j]-xdata[j-1]));
                                                                                                                      A_lat += (1.0*kb*T/H)*ydata[j]*log(p_1)*(xdata[j]-xdata[j-1]);
                                                                                                                      entropy_l += (1.0/H)*((-kb*ydata[j]*log(p_1)*(xdata[j]-xdata[j-1]))+(h*w*xdata[j]*ydata[j]/T)*p_2*(xdata[j]-xdata[j-1]));
                                                                                                                      U_l += (1.0/H)*(h*w*xdata[j]*ydata[j])*(p_2)*(xdata[j]-xdata[j-1]);
                                                                                                                      Cv += (1.0/H)*kb*ydata[j]*(xdata[j]*h*w/(kb*T))*(xdata[j]*h*w/(kb*T))*exp(-h*w*xdata[j]/(kb*T))*(pow(p_3,2.0)*(xdata[j]-xdata[j-1]));
                                                                                                                     }
                                                                                                }
                                                                               }
                                                           }
                                        else if(Gibbs==2){
                                                          A_lat=0.0;
                                                          zpe=0.0;
                                                          entropy_l=0.0;
                                                          Cv=0.0;
                                                          U_l=0.0;
                                                          for(j=0;j<rown;j++){
                                                                              if(ydata[j]>=0.0){
                                                                                                degen=(float)gdata[j];
                                                                                                p_1=1.0-exp(-ydata[i]*h*pow(10.0,12)/(kb*T));
                                                                                                p_2=1.0/(exp(h*ydata[i]*pow(10.0,12)/(kb*T))-1.0);
                                                                                                p_3=(h*ydata[i]*pow(10.0,12)/(kb*T))*(h*ydata[i]*pow(10.0,12)/(kb*T));
                                                                                                zpe += 0.5*(h*degen*ydata[i]*pow(10.0,12)*pow(10.0,19))/1.6;
                                                                                                U_l += h*degen*ydata[i]*pow(10.0,12)*p_2*pow(10.0,19)/1.6;
                                                                                                A_lat += kb*T*degen*log(p_1)*pow(10.0,19)/1.6;
                                                                                                Cv += kb*degen*p_3*exp(h*ydata[i]*pow(10.0,12)/(kb*T))*p_2*p_2*pow(10.0,19)/1.60;
                                                                                                entropy_l += kb*degen*(log(p_1)+(h*ydata[i]*pow(10.0,12)/(kb*T))*exp(-ydata[i]*h*pow(10.0,12)/(kb*T))*(1.0/p_1))*pow(10.0,19)/1.6;
                                                                                                }
                                                                             }
                                                          }
                                        else if(Gibbs==3){
                                                          density=(M*weightt)/(V_x[i-1]*pow(10.0,-30.0));
                                                          Debye=(h/kb)*pow(g2,1.0/3.0)*g11*sqrt(B0*pow(10.0,9.0)/m_f);
                                                          vmax=Debye/T;
                                                          Cv= float(formula)*N_f*3.0*N_a*kb*(4.0*pow(T/Debye,3.0)*simpson(max_in, vmin, vmax)-3.0*(Debye/T)/(exp(Debye/T)-1.0))*(1.0/(9.632*10000.0));
                                                          zpe=float(formula)*(9.0/8)*N_f*kb*Debye/H;
                                                          entropy_l=float(formula)*(-3.0*N_f*kb*log(1.0-exp(-vmax))+4.0*N_f*kb*pow(T/Debye,3.0)*(simpson(max_in, vmin, vmax)))/H;
                                                          U_l=float(formula)*3.0*N_f*kb*T*pow(T/Debye,3.0)*(simpson(max_in, vmin, vmax))/H;
                                                          A_lat=float(formula)*N_f*kb*T*(3.0*log(1.0-exp(-vmax))-pow(T/Debye,3.0)*simpson(max_in, vmin, vmax))/H;
                                                          }
                                       if(bare==1){
                                       	           fprintf(output,"%f\t%lf\t%lf\t%lf\t%lf\n",T,A_lat+zpe+E_x[i-1],entropy_l,U_l+zpe,Cv);
                                                  }
                                       else if(bare==2){
                                       	               fprintf(output,"%f\t%lf\t%lf\t%lf\t%lf\n",T,A_lat+zpe,entropy_l,U_l+zpe,Cv);
                                                      }

                                       A_t=A_lat+A_e+zpe+E_x[i-1];
                                       U_t=U_l+U_e+zpe+E_x[i-1];
                                       entropy_t=entropy_e+entropy_l;
                                       Cvt=Cv+Cel/(23.061*4.189*1000.0);
                                       Ezpe=zpe;
                                       fprintf(output6,"%f\t%lf\t%lf\t%lf\t%lf\t%f\n",T,A_t,U_t,entropy_t,Cv,Cv*9.632*10000.0);
                                       fprintf(output7,"%f\t%lf\t%le\t%lf\t%le\t%lf\t%le\t%d\n",T,A_lat,A_e,U_l,U_e,entropy_l,entropy_e,SMAX);
                                       printf("%f\t%lf\t%10.9f\t%10.9f\t%10.9f\t%10.9f\t%10.9f\t%d\n",T,A_lat,A_e,entropy_e,entropy_l,Cv,Cv*9.632*10000.0,SMAX);
                                       free(sdata); free(sdos); free(sd); free(sdo);
                                     //  fprintf(output8,"%f\t%f\t%lf\n",V_x[i-1],E_x[i-1],zpe);
                                       fprintf(output11,"%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\n",T,V_x[i-1],E_x[i-1],(E_x[i-1]+Ezpe),U_t,entropy_t,Cvt,A_t);
                                       }
                                      fprintf(output8,"%f\t%f\t%lf\n",V_x[i-1],E_x[i-1],Ezpe);
                                      printf("Zero_point_energy(eV/cell)\n");
                                      printf("%lf\n",zpe);
                                      printf("Heat_capacity from 3NR at high temperature(J/K.cell)\n");
                                      printf("%lf\n",3.0*N*8.314);
                                      printf("Total number of different temperatures:%d\n",rowtemp);
                                      fclose(output);
                                      fclose(output1);
                                      fclose(output2);
                                      fclose(output3);
                                      fclose(output5); fclose(output9);
                                      fclose(output6); fclose(output7);
                                      fclose(output10);
                                      free(xdata); free(ydata);
                                      if(Gibbs==2){
                                                   free(gdata);
                                                   }
                                      free(edata); free(dos); free(pdos_e); free(pdos);
                                      printf("Calculations are finished for volume %d.\n",i);
                       }
   fclose(output8); fclose(output11);
// Resort the data, and ouput the EV data at each temperature
   int evoption;
   char filename4[128];
   printf("Do you want to write the E-V curve at each temperature:Yes(1) or No(2)\n");
   printf("Enter your option below:(MUST BE AN INTEGER)\n");
   evoption = 1;
  // scanf("%d",&evoption);

   if(evoption==1){
                   int fileorder;
                   float *tdata, *vdata, *eedata;
                   double *ezpe, *euzpe, *shang, *spheat, *sumall;
                  // FILE *output4;
                  // char filename3[128];
                   printf("Total number of rows in the file %s:%d\n",filename,rowtotal);
                   printf("Number of different cell volumes:%d\n",frames);
                   printf("Number of different temperatures:%d\n",rowtemp);
                   printf("Number of rows computed from the latter two parameters:%d\n",frames*rowtemp);
                   if(rowtotal!=(frames*rowtemp)){
                                                  printf("Inconsistent data information is found.\n");
                                                  printf("Program is terminating.\n");
                                                  printf("Please go back to check your input parameters.\n");
                                                  exit(EXIT_FAILURE);
                                                  }
                   tdata=(float*)malloc(rowtotal*sizeof(float));
                   vdata=(float*)malloc(rowtotal*sizeof(float));
                   eedata=(float*)malloc(rowtotal*sizeof(float));
                   ezpe=(double*)malloc(rowtotal*sizeof(double));
                   euzpe=(double*)malloc(rowtotal*sizeof(double));
                   shang=(double*)malloc(rowtotal*sizeof(double));
                   spheat=(double*)malloc(rowtotal*sizeof(double));
                   sumall=(double*)malloc(rowtotal*sizeof(double));
                   output11=fopen(filename1,"r");
                   for(i=0;i<rowtotal;i++){
                                           fscanf(output11,"%f %f %f %lf %lf %lf %lf %lf",&(tdata[i]),&(vdata[i]),&(eedata[i]),&(ezpe[i]),&(euzpe[i]),&(shang[i]),&(spheat[i]),&(sumall[i]));
                                           }
                   fclose(output11);
                   fileorder=0;
                   for(T=Tmin;T<=Tmax;T=T+step){
                                                fileorder++;
                                                sprintf(filename2,"EV_data_%d.dat",fileorder);
                                                sprintf(filename4,"EV_fit_data_%d.dat",fileorder);
                                                output12=fopen(filename2,"w");
                                                output13=fopen(filename4,"w");
                                                fprintf(output12,"Volume(A^3)\tE0(eV/cell)\tE0+zpe(eV/cell)\tU(eV/cell)\tS(eV/cell)\tCv(eV/cell)\tF(eV/Cell)\n");
                                                fprintf(output12,"Temperature(K):%f\n",T);
                                                for(i=0;i<frames;i++){
                                                                      for(j=0;j<rowtotal;j++){
                                                                                              if(tdata[j]==T&&vdata[j]==V_x[i]){
                                                                                                                                fprintf(output12,"%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\n",vdata[j],eedata[j],ezpe[j],euzpe[j],shang[j],spheat[j],sumall[j]);
                                                                                                                                fprintf(output13,"%f\t%lf\t%lf\n",vdata[j],sumall[j],spheat[j]);
                                                                                                                                }
                                                                                              }
                                                                      }
                                                fclose(output12); fclose(output13);
                                                }
                   free(tdata); free(vdata); free(eedata);
                   free(ezpe); free(euzpe); free(shang); free(spheat); free(sumall);
                   }
   else if(evoption==2){
                        printf("You can exit the program now.\n");
                        exit(0);
                        }
   free(E_x); free(V_x);
   printf("ALL CALCULATIONS ARE DONE FOR FREE ENERGY.\n");
// Nonlinear Least Squares Curve Fitting Program
/* Marquardt algorithm from P.R. Bevington, "Data Reduction and Error Analysis for the Physical Sciences,*/
/* McGraw-Hill, 1969; Adapted by Wayne Weimer & David Harris. Jackknife error algorithm of M.S. Caceci, */
/* Anal. Chem. 1989, 61, 2324. Translated to ANSI C by Dogulas Harris & Tim Seufert 7-94 */
   printf("/*-------------------------------------------------------------------*/\n");
   printf("/*          Non-Linear Least Square Fitting Program                  */\n");
   printf("/*  The algorithm is implemented by Bing Xiao in current program     */\n");
   printf("/*  The only EOS available in the current version: BM-3th order      */\n");
   printf("/*-------------------------------------------------------------------*/\n");
// define parameters and arrays
//   int maxiterations; // maximum iterations
   double param_1[rowtemp], param_2[rowtemp], param_3[rowtemp], param_4[rowtemp]; // fitting parameters
//   double fitp_1[2], fitp_2[2], fitp_3[2], fitp_4[2];
   char filename5[128];
   FILE *output14;
//   double cov; // convergence critia
//   double cov_1, cov_2, cov_3, cov_4;  // Convergence ciriteron
// Star the main program and read the input data at each temperature
   printf("/--------------------------------------------------------------------/\n");
   printf("Least Square Curve Fitting is Initialized.\n");
   printf("'nterms' defines the total number of fitting parameters in expression.\n");
   printf("'func' defines the expression in the fitting.\n");
   printf("Both places must be modified for the new problem.\n");
   printf("/--------------------------------------------------------------------\n");
// Reading the input
//   printf("Enter the maximum iterations: (400)\n");
//   scanf("%d",&maxiterations);
   for(L=1;L<=rowtemp;L++){
                           param_1[L-1]=0.000000; param_2[L-1]=0.000000;
                           param_3[L-1]=0.000000; param_4[L-1]=0.000000;
                           }
   sprintf(filename5,"EOS_fitting_param.dat");
   output14=fopen(filename5,"w");
   FILE *fpp = NULL ;
   char storename[128];
   sprintf(storename,argv[16]);
   fpp = fopen(storename,"w") ;
   for(L=1;L<=rowtemp;L++){
                           npts=frames;
                           sprintf(filename3,"EV_fit_data_%d.dat",L);
                           fp=fopen(filename3,"r");
                           if(fp==NULL){
                                        printf("Unable to open the input file %s\n",filename3);
                                        printf("Program is terminating.\n");
                                        exit(1);
                                        }
                           for(i=0;i<frames;i++){
                                                 fscanf(fp,"%lf %lf %lf",&(x[i]),&(y[i]),&(z[i]));
                                                 }
                           // set all uncertainity of input data as 1.0
                           for(i=0;i<frames;i++){
                                                 sigmay[i]=1.000000;
                                                 weight[i]=1.0/(sigmay[i]*sigmay[i]);
                                                 }
                           fclose(fp);
                           // Printf the input data for checking
                           printf("Data set:%d\n",L);
                           printf("Temperature(K):%f\n",(float)(Tmin+(L-1)*step));
                           printf("Volume(A^3)\tF(ev/cell)\t(+-)Sigma(eV/cell)\tWeight\n");
                           for(i=0;i<frames;i++){
                                                 fprintf(fpp,"%lf\t\t%lf\n",x[i],y[i]);
                                                 printf("%lf\t%lf\t%lf\t%lf\n",x[i],y[i],sigmay[i],weight[i]);
                                                 }

                       /*    // Enter the initial guess for fitting parameters
                           if(L==1){
                                    printf("For the first run, the fitting parameters must be given.\n");
                                    printf("All fitting parameters can not be exact zero (0.000000).\n");
                                    printf("If the value is actually small, please input it as 0.00001.\n");
                                    printf("Parameter 1: E0--->Equilibrium Free Energy.\n");
                                    printf("Parameter 2: V0--->Equilibrium Cell Volume.\n");
                                    printf("Parameter 3: B0--->Iosthermal Bulk Modulus.\n");
                                    printf("Parameter 4: B'--->Derivative of B respect to P.\n");
                                    printf("/##################################################/\n");
                                    printf("/#Please provide the best initial values for them.#/\n");
                                    printf("/##################################################/\n");
                                    printf("Now enter the initial guess for all fitting parameters:\n");
                                    int S;
                                    while((S=getchar())!=EOF&&S!='\n');
                                    for(i=0;i<nterms;i++){
                                                          do{
                                                             printf("Parameter #%d = ", i+1);
                                                             gets(answer);
                                                             }
                                                          while((a[i]=atof(answer))==0.0);
                                                          }
                                    }
                           else if(L>1){
                                        a[0]=param_1[L-2];
                                        a[1]=param_2[L-2];
                                        a[2]=param_3[L-2];
                                        a[3]=param_4[L-2];
                                        printf("Initial guess for fitting parameters(Previous Data Set):\n");
                                        printf("Parameter a[0]=%lf\n",a[0]);
                                        printf("Parameter a[1]=%lf\n",a[1]);
                                        printf("Parameter a[2]=%lf\n",a[2]);
                                        printf("Parameter a[3]=%lf\n",a[3]);
                                        }
                           // Initializing the parameters in fitting
                           flambda=0.000000; iteration=0; cycle=0;
                           do{
                              curvefit(npts);
                              iteration++;
                              // Display the fitting parameters during each iteration
                              printf("\nIteration #%d\n",iteration);
                              for(i=0;i<nterms;i++){
                                                    printf("A[%d]=%lf\n",i,a[i]);
                                                    finala[i]=a[i];
                                                    }
                              printf("Sum of squares of residuals=%lf\n",fchisq*nfree);
                              sy=sqrt(fchisq);
                              // For corase T-grid, we  check the convergence for every isothermal curve
                              if(step>boundary){
                                              if(iteration<2){
                                                              answer[0]='Y';
                                                              }
                                              else{
                                                   printf("When the following critia is fulfilled,\n");
                                                   printf("you should stop the iteration.\n");
                                                   printf("Sum of squares of residuals <= 0.000050.\n");
                                                   printf("Do you want have another iteration:(Y/N)?\n");
                                                   gets(answer);
                                                   }
                                              }
                              // For dense T-gird, we only check the first three volumes
                              else if(step<=boundary){
                                                 if(L<=3){
                                                          if(iteration<1){
                                                                          answer[0]='Y';
                                                                          }
                                                          else{
                                                               printf("When the following critia is fulfilled,\n");
                                                               printf("you should stop the iteration.\n");
                                                               printf("Sum of squares of residuals <= 0.000050.\n");
                                                               printf("Do you want have another iteration:(Y/N)?\n");
                                                               gets(answer);
                                                               }
                                                          }
                                                 else if(L>3){
                                                              // One iteration if sufficient if dense grid is used in the calculation
                                                              if(iteration<1){
                                                                              answer[0]='Y';
                                                                              }
                                                              else{
                                                                   answer[0]='N';
                                                                   }
                                                              }
                                                 }
                              } while(answer[0]!='N');
                            param_1[L-1]=a[0]; param_2[L-1]=a[1];
                            param_3[L-1]=a[2]; param_4[L-1]=a[3];
                            printf("Fitting Parameters in final iteration for data set %d are:\n",L);
                            printf("%lf\t%lf\n%lf\t%lf\n",param_1[L-1],param_2[L-1],param_3[L-1],param_4[L-1]);
                            fprintf(output14,"%f\t%lf\t%lf\t%lf\t%lf\n",(float)(Tmin+(L-1)*step),a[0],a[1],a[2],a[3]);
                            if(step>boundary){
                                         printf("\nDo you want to calculate the uncertainities in fitting parameters:\n");
                                         printf("\nEnter you choice:(Y/N)?\n");
                                         gets(answer);
                                         }
                            else if(step<=boundary){
                                               answer[0]='N';
                                               }
                            if(answer[0]=='Y') uncertainties();*/
                           }
   fclose(fpp);
   fclose(output14);
   return 0;

//------------------------------------------------------------------------------
/*Volumetric Thermal Expansion Coefficient*/
   printf("/*-------------------------------------------------------------------*/\n");
   printf("/*Computing Thermal Expansion Coefficient by Finite Difference Method*/\n");
   printf("/*-------------------------------------------------------------------*/\n");
   double refV; // Reference volume
   float *tmpdata;
   double *eqenergy, *eqvolume, *eqmodulus, *eqderivative;
   double vtec;
   double vtec1;
   double *difcvp; // Cp-Cv terms
   double *difcvv; // Cv in array
   char filename6[128];
   FILE *output15;
   printf("Enter the equilibrium volume at 0K:(A^3/cell)\n");
   scanf("%lf",&refV);
   tmpdata=(float*)malloc(rowtemp*sizeof(float));
   eqenergy=(double*)malloc(rowtemp*sizeof(double));
   eqvolume=(double*)malloc(rowtemp*sizeof(double));
   eqmodulus=(double*)malloc(rowtemp*sizeof(double));
   eqderivative=(double*)malloc(rowtemp*sizeof(double));
   difcvp=(double*)malloc(rowtemp*sizeof(double));
   difcvv=(double*)malloc(rowtemp*sizeof(double));
   output14=fopen(filename5,"r");
   for(i=0;i<rowtemp;i++){
                          fscanf(output14,"%f %lf %lf %lf %lf",&(tmpdata[i]),&(eqenergy[i]),&(eqvolume[i]),&(eqmodulus[i]),&(eqderivative[i]));
                          }
   fclose(output14);
   for(i=0;i<rowtemp;i++){
                          difcvp[i]=0.000000;
                          difcvv[i]=0.000000;
                          }
   sprintf(filename6,"TEC.dat");
   output15=fopen(filename6,"w");
   printf("Temperature(K)\tTEC(K^-1)\n");
   fprintf(output15,"T(K)\tTEC_V0(K^-1)\tTEC_V0(T0)(K^-1)\n");
   for(i=1;i<rowtemp;i++){
                          vtec=(1.0/refV)*(eqvolume[i]-eqvolume[i-1])/(tmpdata[i]-tmpdata[i-1]); // alpha = (1/V0)*dV/dT
                          vtec1=(1.0/eqvolume[i-1])*(eqvolume[i]-eqvolume[i-1])/(tmpdata[i]-tmpdata[i-1]);
                          difcvp[i]=(vtec1*vtec1)*tmpdata[i]*(eqvolume[i]*pow(10.0,-30.0))*(eqmodulus[i]*160.0*pow(10.0,9.0))/(1.60*pow(10.0,-19.0)); // ev/cell.K
                          difcvv[i]=vtec1;
                          printf("%f\t%le\t%le\n",tmpdata[i],vtec,vtec1);
                          fprintf(output15,"%f\t%le\t%le\n",tmpdata[i],vtec,vtec1);
                          }
   fclose(output15);
//------------------------------------------------------------------------------
/*  Calculating the specific heat at constant pressure (Cp)*/
   printf("/*-------------------------------------------------------------------*/\n");
   printf("/*      Computing the Specific Heat at Constant Pressure (Cp)        */\n");
   printf("/*              Adiabatic Bulk Modulus form B_t(T)                   */\n");
   printf("/*           Thermodynamic Gruneisen Parameter(Gamma)                */\n");
   printf("/*-------------------------------------------------------------------*/\n");
   printf("######################################################################\n");
   printf("#The Volume versus Cv data will be interpolated to determine the Cv  #\n");
   printf("#at equilibrium cell volume at each temperature, and then the Cp is  #\n");
   printf("#calculated from: C_p(T)=C_v(T)+alpha^2*V(T)*B(T)*T.                 #\n");
   printf("######################################################################\n");
   double *itpov, *itpofree, *itpocv;  // Import the E-V data at each temperature
   char filename10[128];
   char filename7[128];
   char filename8[128], filename9[128];
   int firstpoint, numpoints;
   int np, kp;
   int intporows;
   double x_s, y_s;
   double itpstep;
   double vleft, vright;
   double cvleft, cvright;
   double cvref;        // Cv
   double cpref;        // Cp
   double modulusref;   // Isobaric bulk modulus
   double gruneisen;
   double *tpdatax, *tpdatay;
   printf("Enter the increment of cell volume in the interpolation:(1.0)\n");
   scanf("%lf",&itpstep);
   FILE *ftp;
   FILE *fftp;
   FILE *output16;
   FILE *output17;
   sprintf(filename9,"Thermal_physicaldata_I.dat");
   sprintf(filename10,"Thermal_physicaldata_II.dat");
   output16=fopen(filename9,"w");
   output17=fopen(filename10,"w");
   printf("Temperature(K)\tVolume(A^3)\tTEC(K^-1)\tCv(eV/cell.K)\tCp(eV/cell.K)\tGamma\n");
   fprintf(output16,"Temperature(K)\tVolume(A^3)\tTEC(K^-1)\tCv(eV/cell.K)\tCp(eV/cell.K)\tGamma\tB_t(eV/A^2)\tB_s(eV/A^2)\n");
   fprintf(output17,"Temperature(K)\tVolume(A^3)\tTEC(K^-1)\tCv(J/mol.K.cell)\tCp(J/mol.K.cell)\tGamma\tB_t(GPa)\tB_s(GPa)\n");
   for(i=1;i<rowtemp;i++){
                          // Import data
                          itpov=(double*)malloc(rowtemp*sizeof(double));
                          itpofree=(double*)malloc(rowtemp*sizeof(double));
                          itpocv=(double*)malloc(rowtemp*sizeof(double));
                          sprintf(filename7,"EV_fit_data_%d.dat",i+1);
                          ftp=fopen(filename7,"r");
                          for(j=0;j<frames;j++){
                                                fscanf(ftp,"%lf %lf %lf",&(itpov[j]),&(itpofree[j]),&(itpocv[j]));
                                                }
                          fclose(ftp);
                          // Interpolating data
                          /*for(kp=0;kp<=np-2;kp++){
                                                 numpoints=frames;
                                                 if(kp==0){
                                                           firstpoint=0;
                                                           }
                                                 else if(kp==np-2){
                                                                  firstpoint=kp-1;
                                                                  }
                                                 else{
                                                      firstpoint=kp-1;
                                                      }
                          */
                          firstpoint=0; numpoints=frames;
                          sprintf(filename8,"V_Cv_interpolation_%d.dat",i+1);
                          fftp=fopen(filename8,"w");
                          intporows=0;
                          for(x_s=itpov[0];x_s<=itpov[frames-1];x_s+=itpstep){
                                                                              intporows++;
                                                                              y_s=inter(x_s, itpov, itpocv, firstpoint, numpoints);
                                                                              fprintf(fftp,"%lf\t%lf\n",x_s,y_s);
                                                                              }
                         fclose(fftp);
                         // find volumes on the left and right hand side of equilibrium one
                         tpdatax=(double*)malloc(intporows*sizeof(double));
                         tpdatay=(double*)malloc(intporows*sizeof(double));
                         fftp=fopen(filename8,"r");
                         for(j=0;j<intporows;j++){
                                                  fscanf(fftp,"%lf %lf",&(tpdatax[j]),&(tpdatay[j]));
                                                  }
                         fclose(fftp);
                         for(j=0;j<intporows;j++){
                                                  double vdiff_1=(tpdatax[j]-eqvolume[i]);
                                                  if(vdiff_1>0.0){
                                                                  vright=tpdatax[j];
                                                                  cvright=tpdatay[j];
                                                                  break;
                                                                  }
                                                  }
                         for(j=intporows-1;j>0;j--){
                                                    double vdiff_2=(tpdatax[j]-eqvolume[i]);
                                                    if(vdiff_2<0.0){
                                                                    vleft=tpdatax[j];
                                                                    cvleft=tpdatay[j];
                                                                    break;
                                                                    }
                                                    }
                         cvref=(cvleft+cvright)/2.0;
                         cpref=cvref+difcvp[i];
                         modulusref = eqmodulus[i]*(cpref/cvref);
                         gruneisen=difcvv[i]*(eqmodulus[i]*160.0*pow(10.0,9.0))*(eqvolume[i]*pow(10.0,-30.0))/(cvref*1.60*pow(10.0,-19.0));
                         printf("%lf\t%lf\t%le\t%lf\t%lf\t%lf\n",tmpdata[i],eqvolume[i],difcvv[i],cvref,cpref,gruneisen);
                         fprintf(output16,"%lf\t%lf\t%le\t%lf\t%lf\t%lf\t%lf\t%lf\n",tmpdata[i],eqvolume[i],difcvv[i],cvref,cpref,gruneisen,eqmodulus[i],modulusref);
                         fprintf(output17,"%lf\t%lf\t%le\t%lf\t%lf\t%lf\t%lf\t%lf\n",tmpdata[i],eqvolume[i],difcvv[i],(cvref*23.061*4.189*1000.0),(cpref*23.061*4.189*1000.0),gruneisen,eqmodulus[i]*160.0,modulusref*160.0);
                         free(tpdatax); free(tpdatay);
                         free(itpov); free(itpofree); free(itpocv);
                         }
   fclose(output16); fclose(output17);
   free(difcvp); free(difcvv);
   printf("Calculations are finsihed for thermal-physical properties.\n");
   printf("Results are saved in two files:\n");
   printf("File 1:%s (unit in eV, eV/A^2)\n",filename9);
   printf("File 2:%s (unit in J, GPa)\n",filename10);
// Clean up some files
   printf("Good time to clean up some adundant files in the folder.\n");
   system("PAUSE");
   int status;
   int status1;
   int status2;
   char delfile[128];
   char delfile1[128];
   char delfile2[128];
   for(L=1;L<=rowtemp;L++){
                           sprintf(delfile,"EV_fit_data_%d.dat",L);
                           status=remove(delfile);
                           if(status==0){
                                         printf("%s is removed.\n",delfile);
                                         }
                           else{
                                printf("Could not find %s in current folder.\n",delfile);
                                }
                           }
   for(L=2;L<=rowtemp;L++){
                           sprintf(delfile1,"V_Cv_interpolation_%d.dat",L);
                           status1=remove(delfile1);
                           if(status1==0){
                                         printf("%s is removed.\n",delfile1);
                                         }
                           else{
                                printf("Could not find %s in current folder.\n",delfile1);
                                }
                           }
   printf("You may also delete all files:EV_data_xxx.dat\n");
   int S1;
   while((S1=getchar())!=EOF&&S1!='\n');
   printf("Do you want to delete them:(Y/N)?\n");
   gets(answer);
   if(answer[0]=='Y'){
                      for(L=1;L<=rowtemp;L++){
                                              sprintf(delfile2,"EV_data_%d.dat",L);
                                              status2=remove(delfile2);
                                              if(status2==0){
                                                             printf("%s is deleted.\n",delfile2);
                                                             }
                                               else{
                                                    printf("Could not find %s is current folder!\n",delfile2);
                                                    }
                                               }
                      }
   free(tmpdata); free(eqenergy); free(eqvolume); free(eqmodulus); free(eqderivative);
   printf("All calculations are done!\n");
   printf("Send email to 'b.xiao@ucl.ac.uk' for questions\n");
   printf("Press 'Enter' key to exit program.\n");
   getchar();getchar();
}
/*------------------------------------------------------------------------------*/
/* Subroutines for Free energy calculations */
float f(double x1, double y1, double z1)
{
      return 1.0/(exp((x1*y1)/z1)+1.0);
      }
float devf(double x1, double y1, double z1)
{
      double K1;
      double K2;
      double K3;
      K1=exp((x1*y1)/z1);
      K2=(1.0/(exp((x1*y1)/z1)+1.0))*(1.0/(exp((x1*y1)/z1)+1.0));
      K3=-1.0*y1/z1;
      return K1*K3*K2;
      }
float De(double x1)
{
   return 3.0*pow(x1,3.0)/(exp(x1)-1.0);
      }
float g(double y1)
{
      return 3.0/(2.0*pow((2.0/3)*(1.0+y1)/(1.0-2.0*y1),3.0/2)+pow((1.0/3)*(1.0+y1)/(1.0-y1),3.0/2));
      }
float f1(double x1)
{
      return 2.0*pow((2.0/3.0)*(1.0+x1)/(1.0-2.0*x1),3.0/2.0);
      }
float f2(double x1)
{
      return pow((1.0/3.0)*(1.0+x1)/(1.0-x1),3.0/2.0);
      }
// Phonon-electron sacttering
float Cut_1(double x1)
{
      return pow(x1,5)/((exp(x1)-1)*(1-exp(-x1)));
      }
// s-d electron scattering for transition metal
float Cut_2(double x1)
{
      return pow(x1,3)/((exp(x1)-1)*(1-exp(-x1)));
      }
// Electron-electron interaction
float Cut_3(double x1)
{
      return pow(x1,2)/((exp(x1)-1)*(1-exp(-x1)));
      }

// Simpson's rule for Debye integration //
float simpson(int no, float min, float max)
{
   int n;
   float interval, sum_1=0., t;
   interval = ((max-min)/(no-1));

   for (n=2; n<no; n+=2)                /* loop for even points  */
   {
       t = interval*(n-1);
       sum_1 += 4*De(t);
   }
   for (n=3; n<no; n+=2)                /* loop for odd points  */
   {
      t = interval*(n-1);
      sum_1 += 2*De(t);
   }
   sum_1 +=  De(min) + De(max);	 	/* add first and last value */
   sum_1 *= interval/3.;        		/* then multilpy by interval*/
   return (sum_1);
}
// Phonon-electron scattering
float simpson_1(int no, float min, float max)
{
   int n;
   float interval, sum_2=0., t;
   interval = ((max-min)/(no-1));

   for (n=2; n<no; n+=2)                /* loop for even points  */
   {
       t = interval*(n-1);
       sum_2 += 4*Cut_1(t);
   }
   for (n=3; n<no; n+=2)                /* loop for odd points  */
   {
      t = interval*(n-1);
      sum_2 += 2*Cut_1(t);
   }
   sum_2 +=  Cut_1(min) + Cut_1(max);	 	/* add first and last value */
   sum_2 *= interval/3.;        		/* then multilpy by interval*/
   return (sum_2);
}
float simpson_2(int no, float min, float max)
{
   int n;
   float interval, sum_3=0., t;
   interval = ((max-min)/(no-1));

   for (n=2; n<no; n+=2)
   {
       t = interval*(n-1);
       sum_3 += 4*Cut_2(t);
   }
   for (n=3; n<no; n+=2)
   {
      t = interval*(n-1);
      sum_3 += 2*Cut_2(t);
   }
   sum_3 +=  Cut_2(min) + Cut_2(max);
   sum_3 *= interval/3.;
   return (sum_3);
}
float simpson_3(int no, float min, float max)
{
   int n;
   float interval, sum_4=0., t;
   interval = ((max-min)/(no-1));

   for (n=2; n<no; n+=2)
   {
       t = interval*(n-1);
       sum_4 += 4*Cut_3(t);
   }
   for (n=3; n<no; n+=2)
   {
      t = interval*(n-1);
      sum_4 += 2*Cut_3(t);
   }
   sum_4 +=  Cut_3(min) + Cut_3(max);
   sum_4 *= interval/3.;
   return (sum_4);
}
/*------------------------------------------------------------------------------*/
/* Suroutines for Non-linear Least Square Fitting to EOS */
double func(int i)
{
 int loop;
 double value;
 if(param==1){
              for(loop=0;loop<nterms;loop++) c[loop]=b[loop];
              }
 else{
      for(loop=0;loop<nterms;loop++) c[loop]=a[loop];
      }
 // c[0]=E0, c[1]=V0, c[2]=B0, c[3]=B'
 double ratio=c[1]/x[i];
 double r1=pow(ratio,2.0/3.0);
 double r2=r1-1.0;
 double r3=6.0-4.0*r1;
 value=c[0]+(9.0*c[1]*c[2]/16.0)*(pow(r2,3.0)*c[3]+pow(r2,2.0)*r3);
 return value;
}
//------------------------------------------------------------------------------
/* Sum of square of differences between meausred and calculated values */
void chisquare(void)
{
     int i;
     fchisq=0;
     for(i=0;i<npts;i++){
                         fchisq += weight[i]*(y[i]-yfit[i])*(y[i]-yfit[i]);
                         }
     fchisq /= nfree;
}
//------------------------------------------------------------------------------
/*  Numberical Derivatives */
void calcderivative(void)
{
     int i, m;
     double atemp, delta;
     for(m=0;m<nterms;m++){
                           atemp = a[m];
                           delta = fabs(a[m]/100000.0);
                           a[m] = atemp + delta;
                           for(i=0;i<npts;i++) deriv[i][m] = (func(i)-yfit[i])/delta;
                           a[m] = atemp;
                           }
}
//------------------------------------------------------------------------------
/* Invert the Matrix array */
void matrixinvert(void)
{
     int i, j, k, ik[nterms], jk[nterms];
     double rsave, amax;
     for(k=0;k<nterms;k++){
                           amax = 0.0;
                           for(i=k;i<nterms;i++){
                                                 for(j=k;j<nterms;j++){
                                                                       if(fabs(amax)<=fabs(arry[i][j])){
                                                                                                      amax = arry[i][j];
                                                                                                      ik[k] = i;
                                                                                                      jk[k] = j;
                                                                                                      }
                                                                       }
                                                 }
                           i = ik[k];
                           if(i>k){
                                   for(j=0;j<nterms;j++){
                                                         rsave = arry[k][j];
                                                         arry[k][j] = arry[i][j];
                                                         arry[i][j] = -1.0*rsave;
                                                         }
                                   }
                           j = jk[k];
                           if(j>k){
                                   for(i=0;i<nterms;i++){
                                                         rsave = arry[i][k];
                                                         arry[i][k] = arry[i][j];
                                                         arry[i][j] = -1.0*rsave;
                                                         }
                                   }
                           for(i=0;i<nterms;i++){
                                                 if(i!=k){
                                                          arry[i][k] = -1.0*arry[i][k]/amax;
                                                          }
                                                 }
                           for(i=0;i<nterms;i++){
                                                 for(j=0;j<nterms;j++){
                                                                       if(j!=k&&i!=k){
                                                                                      arry[i][j] = arry[i][j] + arry[i][k]*arry[k][j];
                                                                                      }
                                                                       }
                                                 }
                           for(j=0;j<nterms;j++){
                                                 if(j!=k){
                                                          arry[k][j] = arry[k][j]/amax;
                                                          }
                                                 }
                           arry[k][k] = 1.0/amax;
                           }
     for(k=nterms-1;k>-1;k--){
                              j=ik[k];
                              if(j>k){
                                      for(i=0;i<nterms;i++){
                                                            rsave = arry[i][k];
                                                            arry[i][k] = -1.0*arry[i][j];
                                                            arry[i][j] = rsave;
                                                            }
                                      }
                              i=jk[k];
                              if(i>k){
                                      for(j=0;j<nterms;j++){
                                                            rsave = arry[k][j];
                                                            arry[k][j] = -1.0*arry[i][j];
                                                            arry[i][j] = rsave;
                                                            }
                                      }
                              }
}
//------------------------------------------------------------------------------
/* Curve fitting Algorithm */
void curvefit(int npoints)
{
     int i, j, k;
     nfree = npoints-nterms;
     for(j=0;j<nterms;j++){
                           b[j] = beta[j] = 0;
                           for(k=0;k<=j;k++) alpha[j][k] = 0;
                           }
     param=0;
     for(i=0;i<npoints;i++) yfit[i] = func(i);
     chisquare();
     chisql = fchisq;
     calcderivative();
     for(i=0;i<npoints;i++){
                            for(j=0;j<nterms;j++){
                                                  beta[j] += weight[i]*(y[i]-yfit[i])*deriv[i][j];
                                                  for(k=0;k<=j;k++)alpha[j][k] += (weight[i]*deriv[i][j]*deriv[i][k]);
                                                  }
                            }
     for(j=0;j<nterms;j++){
                           for(k=0;k<=j;k++) alpha[k][j] = alpha[j][k];
                           }
     nloops = 0;
     do{
        param = 1;
        for(j=0;j<nterms;j++){
                              for(k=0;k<nterms;k++) arry[j][k] = alpha[j][k]/sqrt(alpha[j][j]*alpha[k][k]);
                              arry[j][j]=flambda + 1;
                              }
        matrixinvert();
        for(j=0;j<nterms;j++){
                              b[j] = a[j];
                              for(k=0;k<nterms;k++) b[j] += beta[k]*arry[j][k]/sqrt(alpha[j][j]*alpha[k][k]);
                              }
        for(i=0;i<npoints;i++) yfit[i] = func(i);
        chisquare();
        if((chisql-fchisq)<0) flambda *= 10.0;
        nloops++;
        } while(fchisq > chisql);
     for(j=0;j<nterms;j++){
                           a[j] = b[j];
                           sigmaa[j] = sqrt(arry[j][j]/alpha[j][j]);
                           }
     flambda /= 10.0;
     }
//------------------------------------------------------------------------------
/* Uncertainities */
void uncertainties(void)
{
     int i, k;
     double ajack[nterms][maxnpts], avajack[nterms];
     do{
            cycle++;
            printf("Starting to calculate uncertainities:\n");
            for(i=0;i<npts;i++){
                                jackknifedata(filename3,i);
                                for(k=0;k<=iteration;k++) curvefit(npts-1);
                                printf("Now playing with the data point %d\n",i+1);
                                for(k=0;k<nterms;k++) ajack[k][i] = a[k];
                                }
             printf("\n\n");
             for(k=0;k<nterms;k++) avajack[k] =0;
             for(k=0;k<nterms;k++){
                                   for(i=0;i<npts;i++) avajack[k] += ajack[k][i];
                                   avajack[k] = avajack[k]/npts;
                                   }
             for(k=0;k<nterms;k++) sigmaa[k]=0;
             for(k=0;k<nterms;k++){
                                   for(i=0;i<npts;i++)
                                   sigmaa[k] += (ajack[k][i]-avajack[k])*(ajack[k][i]-avajack[k]);
                                   sigmaa[k] = sqrt((npts-1)*sigmaa[k]/npts);
                                   printf("Parameter[%d]=%lf +- %lf\n",k,finala[k],sigmaa[k]);
                                   if(cycle>1) printf("\t(Previous uncertainty = %lf)\n\n",lastsigmaa[k]);
                                   lastsigmaa[k] = sigmaa[k];
                                   }
             printf("Standard deviation of y = %lf\n",sy);
             printf("Above result is based %d iterations.\n",iteration);
             iteration+=5;
             printf("Iteration will now be increased to %d to see if the estimates of uncertainity change.\n",iteration);
             printf("If the new values are similar to previous results.\n");
             printf("It is the time to stop the calculation.\n");
             printf("\tDo you want to try another cycle now:(Y/N)?\n");
             gets(answer);
             } while(answer[0]=='Y');
}
//------------------------------------------------------------------------------
void jackknifedata(char *filename3, int k)
{
     int n=0;
     double dumy1, dumy2, dumy3;
     fp = fopen(filename3,"rb");
     while(!feof(fp)){
                      fscanf(fp,"%lf %lf %lf",&dumy1,&dumy2,&dumy3);
                      //fread(&x[n],sizeof(double),1,fp);
                      //fread(&y[n],sizeof(double),1,fp);
                      x[n]=dumy1;
                      y[n]=dumy2;
                      z[n]=dumy3;
                      sigmay[n] = 1.0;
                      weight[n] = 1.0/(sigmay[n]*sigmay[n]);
                      n++;
                      }
     npts = n-1;
     fclose(fp);
     for(n=0;n<(npts-1);n++){
                             if(n>=k){
                                      x[n] = x[n+1]; y[n] = y[n+1];
                                      weight[n] = weight[n+1];
                                      }
                             }
     }
//------------------------------------------------------------------------------
