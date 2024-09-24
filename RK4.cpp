#include<fftw3.h> // for using fftw library
#include<iostream> // for using cerr or cout
#include<algorithm> // for using max function
#include<cstdlib> // for using atof function
#include<cmath> // for using math functions
#include<limits> // for using numeric limits
#include<vector> // for using vectors
#include<fstream> //for using files
#include<sstream>
using namespace std;

int N = 128;
int init_iters = 1000;
int iter_write = 1000;
int iter_max = 100000;

double L = 100;
double vA = 1.0;
double vB = 1.0;
double vS = 1.0;
double XAB = -2.0;
double XAS = 3.0;
double XBS = 1.0;
double zA = 1.0;
double zB = -1.0;
double wA = 1.0;
double wB = 1.0;
double wS = 1.0;
double ep = 0.11429; //ep = ep_r/4*Pi*ell_B with ep_r = 80 and ell_B = 55.7
double gammaAS = 1;
double gammaBS = 1;

double dt1 = 0.0001;
double dt2 = 0.01;

double kaA = 1.0;
double kaB = 1.0;

double psi_0 = 0.7;
double psi_L = 0.0;

double dx = L/(double)N;

double *x, *cA, *cA_tmp, *cA_diff, *cB, *cB_tmp, *cB_diff, *psi, *psihat, *muA, *muB, *jA, *jB, *kA, *kB;

fftw_plan pfpsi, pbpsi;

void initialize_structures()
{
    x = (double*) fftw_malloc(sizeof(double) * N);
    cA = (double*) fftw_malloc(sizeof(double) * N);
    cA_tmp = (double*) fftw_malloc(sizeof(double) * N);
    cA_diff = (double*) fftw_malloc(sizeof(double) * N);
    cB = (double*) fftw_malloc(sizeof(double) * N);
    cB_tmp = (double*) fftw_malloc(sizeof(double) * N);
    cB_diff = (double*) fftw_malloc(sizeof(double) * N);
    psi = (double*) fftw_malloc(sizeof(double) * N);
    psihat = (double*) fftw_malloc(sizeof(double) * N);
    muA = (double*) fftw_malloc(sizeof(double) * N);
    muB = (double*) fftw_malloc(sizeof(double) * N);
    jA = (double*) fftw_malloc(sizeof(double) * (N+1));
    jB = (double*) fftw_malloc(sizeof(double) * (N+1));
    kA = (double*) fftw_malloc(sizeof(double) * N);
    kB = (double*) fftw_malloc(sizeof(double) * N);

    pfpsi = fftw_plan_r2r_1d(N, psi, psihat, FFTW_REDFT10, FFTW_ESTIMATE);
    pbpsi = fftw_plan_r2r_1d(N, psihat, psi, FFTW_REDFT01, FFTW_ESTIMATE);
}

void initialize_profile()
{
    for (int j = 0; j < N; j++)
    {
        x[j] = (j + 0.5) * dx - L/2;

        if (j < N/2)
        {
            cA[j] = 0.43;
            cB[j] = 0.43;
        }
        else
        {
            cA[j] = 0.07;
            cB[j] = 0.07;
        }
    }
}

void calculate_norm(double* psi)
{
    for (int j = 0; j < N; j++)
    {
        psi[j] = psi[j]/(2.*(double)N);
    }
}

void calculate_psi(double* cA, double* cB)
{
    for (int j = 0; j < N; j++)
    {
        psi[j] = - zA*cA[j] - zB*cB[j];
    }

    fftw_execute(pfpsi);

    for (int j = 1; j < N; j++)
    {
        double q = M_PI*((double)j)/L;
        psihat[j] /= - q*q*ep;
    }

    fftw_execute(pbpsi);
    calculate_norm(psi);

    double a = (psi[0] - psi[N-1] - psi_0 + psi_L) / ((double)N-1);
    double b = psi_0 - psi[0];

    for (int j = 0; j < N; j++)
    {
        psi[j] += a*j + b;
    }
}

void calculate_mu(double* cA, double* cB)
{
    for (int j = 0; j < N; j++)
    {
        double lap_cA, lap_cB;
        if(j == 0)
        {
            lap_cA = (cA[0] - 2.*cA[1] + cA[2])/(dx*dx);
            lap_cB = (cB[0] - 2.*cB[1] + cB[2])/(dx*dx);
        }
        else if (j == N-1)
        {
            lap_cA = (cA[N-1] - 2.*cA[N-2] + cA[N-3])/(dx*dx);
            lap_cB = (cB[N-1] - 2.*cB[N-2] + cB[N-3])/(dx*dx);
        }
        else
        {
            lap_cA = (cA[j+1] -2.*cA[j] +cA[j-1])/(dx*dx);
            lap_cB = (cB[j+1] -2.*cB[j] +cB[j-1])/(dx*dx);
        }
        muA[j] = log(vA*cA[j]) + 1. - (vA/vS)*log(1 - vA*cA[j] - vB*cB[j]) - (vA/vS) + wA - wS*(vA/vS) + XAB*cB[j] + XAS*(1 - 2.*vA*cA[j] - vB*cB[j])/vS - XBS*cB[j]*(vA/vS) - kaA*lap_cA + zA*psi[j];
        muB[j] = log(vB*cB[j]) + 1. - (vB/vS)*log(1 - vA*cA[j] - vB*cB[j]) - (vB/vS) + wB - wS*(vB/vS) + XAB*cA[j] + XBS*(1 - vA*cA[j] - 2.*vB*cB[j])/vS - XAS*cA[j]*(vB/vS) - kaB*lap_cB + zB*psi[j];
    }
}

void calculate_j(double* cA, double* cB)
{
    jA[0] = 0.0;
    jA[N] = 0.0;
    jB[0] = 0.0;
    jB[N] = 0.0;
    for (int j = 1; j < N; j++)
    {
        double cA_avg = 0.5*(cA[j-1]+cA[j]);
        double cB_avg = 0.5*(cB[j-1]+cB[j]);
        jA[j] = gammaAS*vA*cA_avg*(1-vA*cA_avg-vB*cB_avg)*(muA[j-1]-muA[j])/dx;
        jB[j] = gammaBS*vB*cB_avg*(1-vA*cA_avg-vB*cB_avg)*(muB[j-1]-muB[j])/dx;
    }
}

void f_RK(double* cA, double* cB)
{
    calculate_psi(cA, cB);

    calculate_mu(cA, cB);
    calculate_j(cA, cB);

    for (int j = 0; j < N; j++)
    {
        kA[j] = - (jA[j+1] - jA[j]) / dx;
        kB[j] = - (jB[j+1] - jB[j]) / dx;
    }
}

void run_numerics()
{
    double t = 0;

    for (int iter = 1; iter <= iter_max; iter++) //time iteration
    {
        double dt;

        if (iter <= init_iters)
            dt = dt1;
        else
            dt = dt2;

        //start of RK4 scheme

        f_RK(cA, cB); //k's are filled with k1

        for (int j = 0; j < N; j++)
        {
            cA_diff[j] = dt*kA[j]/6.;
            cB_diff[j] = dt*kB[j]/6.;
        }

        for (int j = 0; j < N; j++)
        {
            cA_tmp[j] = cA[j] + 0.5*dt*kA[j];
            cB_tmp[j] = cB[j] + 0.5*dt*kB[j];
        }

        f_RK(cA_tmp, cB_tmp); //k's are filled with k2

        for (int j = 0; j < N; j++)
        {
            cA_diff[j] += dt*kA[j]/3.;
            cB_diff[j] += dt*kB[j]/3.;
        }

        for (int j = 0; j < N; j++)
        {
            cA_tmp[j] = cA[j] + 0.5*dt*kA[j];
            cB_tmp[j] = cB[j] + 0.5*dt*kB[j];
        }

        f_RK(cA_tmp, cB_tmp); //k's are filled with k3

        for (int j = 0; j < N; j++)
        {
            cA_diff[j] += dt*kA[j]/3.;
            cB_diff[j] += dt*kB[j]/3.;
        }

        for (int j = 0; j < N; j++)
        {
            cA_tmp[j] = cA[j] + dt*kA[j];
            cB_tmp[j] = cB[j] + dt*kB[j];
        }

        f_RK(cA_tmp, cB_tmp); //k's are filled with k4

        for (int j = 0; j < N; j++)
        {
            cA[j] += cA_diff[j] + dt*kA[j]/6.;
            cB[j] += cB_diff[j] + dt*kB[j]/6.;
        }

        //end of RK4 scheme

        t += dt;

        if (iter % iter_write == 0)
        {
            calculate_psi(cA, cB);

            stringstream name;
            name << "data_" << iter/iter_write << ".txt";
            ofstream file (name.str().c_str());
            file.precision(numeric_limits<double>::digits10);
            file << "#vA = " << vA << endl;
            file << "#vB = " << vB << endl;
            file << "#vS = " << vS << endl;
            file << "#XAB = " << XAB << endl;
            file << "#XAS = " << XAS << endl;
            file << "#XBS = " << XBS << endl;
            file << "#wA = " << wA << endl;
            file << "#wB = " << wB << endl;
            file << "#wS = " << wS << endl;
            file << "#zA = " << zA << endl;
            file << "#zB = " << zB << endl;
            file << "#kaA = " << kaA << endl;
            file << "#kaB = " << kaB << endl;
            file << "#ep = " << ep << endl;
            file << "#dt1 = " << dt1 << endl;
            file << "#dt2 = " << dt2 << endl;
            file << "#init_iters = " << init_iters << endl;
            file << "#iter_write = " << iter_write << endl;
            file << "#iter_max = " << iter_max << endl;
            file << "#L = " << L << endl;
            file << "#N = " << N << endl;
            file << "#dx = " << dx << endl;
            file << "#t = " << t << endl;
            file << "#" << endl;
            file << "# x  cA(x)  cB(x)  psi(x)" << endl;

            for(int j = 0; j < N; ++j)
            {
                file << x[j] << " " << cA[j] << " " << cB[j] << " " << psi[j] << endl;
            }
            file.clear();
            file.close();
        }
    }
}

void clean_up()
{
    fftw_destroy_plan(pfpsi);
    fftw_destroy_plan(pbpsi);

    fftw_free(x);
    fftw_free(cA);
    fftw_free(cA_tmp);
    fftw_free(cA_diff);
    fftw_free(cB);
    fftw_free(cB_tmp);
    fftw_free(cB_diff);
    fftw_free(psi);
    fftw_free(psihat);
    fftw_free(muA);
    fftw_free(muB);
    fftw_free(jA);
    fftw_free(jB);
    fftw_free(kA);
    fftw_free(kB);
}

int main(int argc, char *argv[])
{
    initialize_structures();
    initialize_profile();
    run_numerics();

    clean_up();

    return 0;
}
