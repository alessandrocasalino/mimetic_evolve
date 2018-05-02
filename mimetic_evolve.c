// Alessandro Casalino
//
// Paper with equations: arXiv:1803.02620
//
// Compile with gcc-7 -O2 mimetic_evolve.c -o mimetic_evolve.exe
// Run with ./mimetic_evolve.exe
//
//
// NOTE: The scale factor is defined as F, and the derivative of the scale factor as A
// NOTE: Units of measurement with 16 \pi G = 1, as in the paper
//
// NOTE: inline use for functions
// Speed up computation but need to compile gcc-7 (Homebrew on macOS), or alternatively use -O2 with clang (might fail to compile though)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>




// PHYSICAL PARAMETERS VALUES

// Action parameters
double a=0.0;
double b=0.0;
double c=0.0;

// INITIAL CONDITIONS
// Initial value of the scale factor
double F_init = 1e-5;

// Values of fractional density for the cosmological matter today
double OmegaR_0 = 8e-5;
double OmegaB_0 = 0.0486;
double OmegaLambda_0 = 0.6911;
double OmegaCDM_0 = 0.2589;

// Physical constants for output conversions
double _G_ = 6.67428e-11;                 /**< Newton constant in m^3/Kg/s^2 */
double _MPc_over_m_ = 3.085677581282e22;  // Conversion factor from Mpc to m
double _Gyr_over_Mpc_ = 3.06601394e2;     // Conversion factor from Gyr to Mpc
double _c_ = 2.99792458e8;                // Speed of light in m/s
double _H0_ = 67.74;                      // H0 of LCDM in Km/s/Mpc

// MODEL SELECTION - POTENTIAL SELECTION (write the correspoding int)
//  0: No potential (this means only rhob, rhor and eventually rhoCDM can be modelled)
//  1: LCMD (dark fluid exactly as dark matter + dark energy)
//  2: QUINTESSENCE (i.e. 1/\phi^2 potential)
//  3: PHANTOM
int model = 1;
// Potential constants (makes something only for model>1)
double V_const = 2e3;
double ts = 1e10;


// COMPUTATION PARAMETERS VALUES

// Number of points used in the computation
int points = (int) 1e6;

// Number of c values to explore in the physical range -2/3<c<=0 (only c=0 if 1)
int c_values = 1;

// Raise this value to make the csv file smaller, but decreasing resolution
// The value inserted is the ratio between the number of values written in a full resolution file / decreased resolution file
int csv_resolution = 10;

double delta_bisection = 1e-2;
double rhof_init_min = 1e-10;
double rhof_init_max = 1e16; //1e16 max for LCDM (model 1)


// PERTURBATION FACTOR COMPUTATION
// Provides the factor of \delta\phi of equation (22) (arXiv:1803.02620) as output in .csv file
int PERT_FACT_MODE = 1;

// TEST MODE
// Provides some informations on the terminal and the .csv file during the computation (1: on, others: off)
int TEST_MODE = 0;


// Definition of the POTENTIAL as function of t (remember t = \phi due to mimetic constraint)
// For a list of the models see above
inline double V(double t) {

  double result;

  if (model==1) {
    result = 6. * OmegaLambda_0;
  }
  else if (model==2) {
    result = V_const / t / t;
  }
  else if (model==3) {
    result = V_const / (ts-t) / (ts-t);
  }
  else {
    result = 0.;
  }

  return result;

}

// Definition of the velocity of scalar perturbations
inline double cs2(double a, double b, double c){

  return 2.*(b-c)*(a-1.)/(2.*a-b-2.)/(4.-4.*a-b+3.*c);

}

// Definitions of DIFFERENTIAL EQUATION system
// Need to have first order equations to use Runge Kutta
//
// This is the derivative of the scale factor F defined as A
inline double D1(double t, double F, double A, double rhof, double rhor, double rhob, double a, double b, double c) {

    return A;

}

inline double D2(double t, double F, double A, double rhof, double rhor, double rhob, double a, double b, double c) {

    return -1./F/2.*A*A-F/(4.-b+3.*c-4.*a)*(1./3.*rhor-V(t));

}

inline double D3(double t, double F, double A, double rhof, double rhor, double rhob, double a, double b, double c) {

    return -3. * A/F * (rhof-V(t));

}

inline double D4(double t, double F, double A, double rhof, double rhor, double rhob, double a, double b, double c) {

    return -4. * A/F * rhor;

}

inline double D5(double t, double F, double A, double rhof, double rhor, double rhob, double a, double b, double c) {

    return -3. * A/F * rhob;

}


// Function used to print results stored in vectors as a csv file
void csv(double * t, double * F, double * A, double * rhof, double * rhor, double * rhob, double a, double b, double c, char * filename) {

    int i = 0;

    FILE *fp;
    fp = fopen (filename, "w+");
    fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,%s", "t [Gyr]", "a(t)", "H(t) [Km/s/Mpc]", "H_prime(t) [(Km/s/Mpc)^2]", "Omega_df", "Omega_r", "Omega_b", "omega_df");
    if(PERT_FACT_MODE==1) fprintf(fp, ",%s", "Pert. factor [Km^2/s^2/Mpc^2]");
    if(TEST_MODE==1) fprintf(fp, ",%s", "H_check [Km/s/Mpc]");
    fprintf(fp, "\n");

    double MPl = (4.-b+3.*c-4.*a)/4.;
    double CS2 = cs2(a,b,c);
    double k_H0 = 1.;
    double H_c = 0.,pert_factor = 0.;

    // Time conversion factor
    double tcf = 1./_H0_/(60.*60.*24.*365.*1e9)*_MPc_over_m_/1000.;

    while( F[i] <= 1.0){

      double H = sqrt(1./MPl)*sqrt(rhor[i]/6.+rhof[i]/6.+rhob[i]/6.);
      double H_prime = - 3./2. * H * H - rhor[i]/MPl/3./4. + V(t[i])/MPl/4.;

      fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e", t[i] * tcf, F[i], H * _H0_, H_prime * _H0_ * _H0_, rhof[i]/MPl/H/H/6., rhor[i]/MPl/H/H/6., rhob[i]/MPl/H/H/6., -V(t[i])/rhof[i]);

      if(PERT_FACT_MODE==1){

        pert_factor = H_prime * _H0_ * _H0_ + 2.*(a-1.)/(2.*a-b-2.)/(4.-4.*a-b+3.*c)*(4./3.*rhor[i]+rhob[i]) * _H0_ * _H0_ + k_H0 * k_H0 * CS2/F[i]/F[i];

        fprintf(fp, ",%e", pert_factor);

      }

      if(TEST_MODE==1){

        H_c = A[i]/F[i];

        fprintf(fp, ",%e", H_c * _H0_);

      }

      fprintf(fp, "\n");

      if(F[i+csv_resolution]<= 1.0){
        i=i+csv_resolution;
      }
      else{
        i++;
      }

    }

    fclose(fp);

}

int scan_for_F0 (double * F) {

  int i=0;

  while( F[i] <= 1.0){
    i=i+1;
  }

  return i;

}

// This divides an interval in a logarithmic scale of basis 10, and store results in input vector v
void logscale10 (double * v, double A, double B, int points){

    int i;

    double a = log10(A);
    double b = log10(B);

    double h = (b - a) / (points-1.0);

    for(i=0;i<points;i++){
        v[i] = pow(10.0, a + i * h);
    }

}

// This function is a temporal step of the Runge-Kutta 4 method:
// evolves the system for a time equal to DeltaT (h in the program)
void rk4_step(double * f, double a, double b, double c, double h) {

    double k_F[4], k_A[4], k_rhof[4], k_rhor[4], k_rhob[4];
    double t = f[0];
    double F = f[1];
    double A = f[2];
    double rhof = f[3];
    double rhor = f[4];
    double rhob = f[5];

    k_F[0] = D1(t, F, A, rhof, rhor, rhob, a, b, c);
    k_A[0] = D2(t, F, A, rhof, rhor, rhob, a, b, c);
    k_rhof[0] = D3(t, F, A, rhof, rhor, rhob, a, b, c);
    k_rhor[0] = D4(t, F, A, rhof, rhor, rhob, a, b, c);
    k_rhob[0] = D5(t, F, A, rhof, rhor, rhob, a, b, c);

    k_F[1] = D1(t + 0.5 * h, F + 0.5 * h * k_F[0], A + 0.5 * h * k_A[0], rhof + 0.5 * h * k_rhof[0], rhor + 0.5 * h * k_rhor[0], rhob + 0.5 * h * k_rhob[0], a, b, c);
    k_A[1] = D2(t + 0.5 * h, F + 0.5 * h * k_F[0], A + 0.5 * h * k_A[0], rhof + 0.5 * h * k_rhof[0], rhor + 0.5 * h * k_rhor[0], rhob + 0.5 * h * k_rhob[0], a, b, c);
    k_rhof[1] = D3(t + 0.5 * h, F + 0.5 * h * k_F[0], A + 0.5 * h * k_A[0], rhof + 0.5 * h * k_rhof[0], rhor + 0.5 * h * k_rhor[0], rhob + 0.5 * h * k_rhob[0], a, b, c);
    k_rhor[1] = D4(t + 0.5 * h, F + 0.5 * h * k_F[0], A + 0.5 * h * k_A[0], rhof + 0.5 * h * k_rhof[0], rhor + 0.5 * h * k_rhor[0], rhob + 0.5 * h * k_rhob[0], a, b, c);
    k_rhob[1] = D5(t + 0.5 * h, F + 0.5 * h * k_F[0], A + 0.5 * h * k_A[0], rhof + 0.5 * h * k_rhof[0], rhor + 0.5 * h * k_rhor[0], rhob + 0.5 * h * k_rhob[0], a, b, c);

    k_F[2] = D1(t + 0.5 * h, F + 0.5 * h * k_F[1], A + 0.5 * h * k_A[1], rhof + 0.5 * h * k_rhof[1], rhor + 0.5 * h * k_rhor[1], rhob + 0.5 * h * k_rhob[1], a, b, c);
    k_A[2] = D2(t + 0.5 * h, F + 0.5 * h * k_F[1], A + 0.5 * h * k_A[1], rhof + 0.5 * h * k_rhof[1], rhor + 0.5 * h * k_rhor[1], rhob + 0.5 * h * k_rhob[1], a, b, c);
    k_rhof[2] = D3(t + 0.5 * h, F + 0.5 * h * k_F[1], A + 0.5 * h * k_A[1], rhof + 0.5 * h * k_rhof[1], rhor + 0.5 * h * k_rhor[1], rhob + 0.5 * h * k_rhob[1], a, b, c);
    k_rhor[2] = D4(t + 0.5 * h, F + 0.5 * h * k_F[1], A + 0.5 * h * k_A[1], rhof + 0.5 * h * k_rhof[1], rhor + 0.5 * h * k_rhor[1], rhob + 0.5 * h * k_rhob[1], a, b, c);
    k_rhob[2] = D5(t + 0.5 * h, F + 0.5 * h * k_F[1], A + 0.5 * h * k_A[1], rhof + 0.5 * h * k_rhof[1], rhor + 0.5 * h * k_rhor[1], rhob + 0.5 * h * k_rhob[1], a, b, c);

    k_F[3] = D1(t + h, F + h * k_F[2], A + h * k_A[2], rhof + h * k_rhof[2], rhor + h * k_rhor[2], rhob + h * k_rhob[2], a, b, c);
    k_A[3] = D2(t + h, F + h * k_F[2], A + h * k_A[2], rhof + h * k_rhof[2], rhor + h * k_rhor[2], rhob + h * k_rhob[2], a, b, c);
    k_rhof[3] = D3(t + h, F + h * k_F[2], A + h * k_A[2], rhof + h * k_rhof[2], rhor + h * k_rhor[2], rhob + h * k_rhob[2], a, b, c);
    k_rhor[3] = D4(t + h, F + h * k_F[2], A + h * k_A[2], rhof + h * k_rhof[2], rhor + h * k_rhor[2], rhob + h * k_rhob[2], a, b, c);
    k_rhob[3] = D5(t + h, F + h * k_F[2], A + h * k_A[2], rhof + h * k_rhof[2], rhor + h * k_rhor[2], rhob + h * k_rhob[2], a, b, c);

    f[1] = f[1] + 1.0 / 6.0 * h * ( k_F[0] + 2.0 * k_F[1] + 2.0 * k_F[2] + k_F[3]);
    f[2] = f[2] + 1.0 / 6.0 * h * ( k_A[0] + 2.0 * k_A[1] + 2.0 * k_A[2] + k_A[3]);
    f[3] = f[3] + 1.0 / 6.0 * h * ( k_rhof[0] + 2.0 * k_rhof[1] + 2.0 * k_rhof[2] + k_rhof[3]);
    f[4] = f[4] + 1.0 / 6.0 * h * ( k_rhor[0] + 2.0 * k_rhor[1] + 2.0 * k_rhor[2] + k_rhor[3]);
    f[5] = f[5] + 1.0 / 6.0 * h * ( k_rhob[0] + 2.0 * k_rhob[1] + 2.0 * k_rhob[2] + k_rhob[3]);
    f[0] = f[0] + h;

}

// This function evolves the system with the Runge-Kutta 4 method until t_stop
double rk4(double * t, double * F, double * A, double * rhof, double * rhor, double * rhob, double a, double b, double c, double OmegaDf_0) {

    double f[6];
    int i=0;

    int j = 0;

    double MPl = (4.-b+3.*c-4.*a)/4.;

    // Call the function for a logarithmic scale of time
    // We don't start from t = 0 to avoid problems with quintessence potentials
    logscale10(t,1e-10,3.,points);

    while( i < points-1){

      f[0] = t[i];
      f[1] = F[i];
      f[2] = A[i];
      f[3] = rhof[i];
      f[4] = rhor[i];
      f[5] = rhob[i];

      rk4_step(f, a, b, c, t[i+1]-t[i]);

      F[i+1] = f[1];
      A[i+1] = f[2];
      rhof[i+1] = f[3];
      rhor[i+1] = f[4];
      rhob[i+1] = f[5];
      i++;

      if (F[i]<=1.){
        j=i;
      }

    }

    double H = sqrt(1./MPl) * sqrt(rhor[j]/6.+rhof[j]/6.+rhob[j]/6.);
    return rhof[j]/MPl/H/H/6.-OmegaDf_0;

}


double bisection (double min, double max, double * t, double * F, double * A, double * rhof, double * rhor, double * rhob, double a, double b, double c, double OmegaDf_0) {

  double C = (min+max) / 2.;

  double MPl = (4.-b+3.*c-4.*a)/4.;

  while(fabs((max-min)/min)>delta_bisection){

    double rhor_init = 6. * OmegaR_0/pow(F_init,4.);
    double rhob_init = 6. * OmegaB_0/pow(F_init,3.);

    A[0] = F_init * sqrt(1./MPl)*sqrt(OmegaR_0/pow(F_init,4.) + min/6. + OmegaB_0/pow(F_init,3.));
    rhof[0] = min;
    double rk4_min = rk4(t, F, A, rhof, rhor, rhob, a, b, c, OmegaDf_0);
    A[0] = F_init * sqrt(1./MPl)*sqrt(OmegaR_0/pow(F_init,4.) + C/6. + OmegaB_0/pow(F_init,3.));
    rhof[0] = C;
    double rk4_C = rk4(t, F, A, rhof, rhor, rhob, a, b, c, OmegaDf_0);

    if(rk4_min*rk4_C>=0) {
      min=C;
    }
    else {
      max=C;
    }

    C=(max+min)/2.;

    if (TEST_MODE == 1) printf("TEST_MODE ON - min: %e , max: %e, C: %e, rk4_min: %e , rk4_C: %e \n", min, max, C, rk4_min, rk4_C);

  }

  double result = (max+min)/2.;

  if(TEST_MODE==1) printf("\n");
  printf("\t-> Result of bisection method is rho_df: %e (internal units).\n", result);
  if(TEST_MODE==1) printf("\t--> Confront with LCDM value: %e (internal units).\n", 6. * OmegaLambda_0 + 6. * OmegaCDM_0 / pow(F_init,3.));

  return result;

}

int main() {

    double rhor_init = 6. * OmegaR_0/pow(F_init,4.);
    double rhob_init = 6. * OmegaB_0/pow(F_init,3.);

    double OmegaDf_0 = 1. - OmegaR_0 - OmegaB_0;

    double MPl = (4.-b+3.*c-4.*a)/4.;

    int i = 0;
    while(i<c_values){

        printf("\n\t\t----------------------------------------\n\n");

        // Definition of the vector needed for the evolution functions
        double * t; double * F; double * A; double * rhof; double * rhor; double * rhob;
        t = (double *) malloc(sizeof(double) * points);
        F = (double *) malloc(sizeof(double) * points);
        A = (double *) malloc(sizeof(double) * points);
        rhof = (double *) malloc(sizeof(double) * points);
        rhor = (double *) malloc(sizeof(double) * points);
        rhob = (double *) malloc(sizeof(double) * points);

        // Initial conditions (t=0)
        F[0] = F_init;
        rhor[0] = rhor_init;
        rhob[0] = rhob_init;

        if(model!=0){
          printf(" Searching for best initial value for the dark fluid ...\n \n");
          rhof[0] = bisection (rhof_init_min, rhof_init_max, t, F, A, rhof, rhor, rhob, a, b, c, OmegaDf_0);
          A[0] = F_init * sqrt(1./MPl) * sqrt(rhor_init/6.+rhof[0]/6.+rhob_init/6.);
        }
        else{
          printf(" Model with no potential. Consider CDM without dark energy. \n");
          rhof[0] = OmegaCDM_0/pow(F_init,3.);
          A[0] = F_init * sqrt(1./MPl) * sqrt(rhor_init/6.+rhof[0]/6.+rhob_init/6.);
        }

        printf("\n\n Evolving the system with Runge-Kutta 4 ...\n");
        printf(" Action parameters used are:\n");
        printf("  a: %2.3f\n", a);
        printf("  b: %2.3f\n", b);
        printf("  c: %2.3f\n", c);

        rk4(t, F, A, rhof, rhor, rhob, a, b, c, OmegaDf_0);

        char filename[50];
        sprintf (filename, "RK4_c_%0.2f_m_%d.csv", c, model);
        int last_int = scan_for_F0(F);

        printf("\n RESULTS:\n");
        printf("\t-> H0: %f km/s/Mpc\n", A[last_int]/F[last_int] * _H0_);
        printf("\t-> Age of the Universe: %f Gyr\n", t[last_int] /_H0_/(60.*60.*24.*365.*1e9)*_MPc_over_m_/1000.);
        printf("\t-> Velocity fo scalar perturbations (cs2): %f \n", cs2(a,b,c));

        csv(t, F, A, rhof, rhor, rhob, a, b, c, filename);

        printf("\n The results are saved in '.csv' files. The name is labelled with the value of c, and the model (m).\n");

        if(TEST_MODE==1) printf("\n TEST_MODE ON: check the values of H in .csv file. They must be equal!");

        printf("\n\t\t----------------------------------------\n");

        free(t);free(F);free(A);free(rhof);free(rhor);free(rhob);

        i++;
        c=c-0.1;

    }

    printf("\n");

    exit(0);

}
