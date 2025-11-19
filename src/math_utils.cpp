// src/math_utils.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include "math_utils.h"

void physicalMesh(const std::vector<double>& Xc, 
                  const std::vector<double>& Yc,
                  double k, double ymax,
                  std::vector<std::vector<double>>& Xc_p,
                  std::vector<std::vector<double>>& Yc_p) {
    
    int Nx = Xc.size();
    int Ny = Yc.size();
    
    // 重置并重新分配内存
    Xc_p.assign(Nx, std::vector<double>(Ny, 0.0));
    Yc_p.assign(Nx, std::vector<double>(Ny, 0.0));
    
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Xc_p[i][j] = Xc[i];
            
            if (Xc[i] < 1.0) {
                Yc_p[i][j] = Yc[j];
            } else {
                Yc_p[i][j] = Yc[j] + (1.0 - Yc[j] / ymax) * k * (Xc[i] - 1.0);
            }
        }
    }
}

void uInitialFunc(const std::vector<std::vector<double>>& Xc_p,
                  const std::vector<std::vector<double>>& Yc_p,
                  std::vector<std::vector<std::vector<double>>>& U,
                  double gam){

    int Nx = Xc_p.size();
    int Ny = Yc_p[0].size();

    std::vector<std::vector<double>> rho(Nx,std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> u(Nx,std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> v(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> p(Nx,std::vector<double>(Ny,0.0));
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            rho[i][j]=99719/(287.14*293.15);
            u[i][j]=686.47; 
            v[i][j]=0.0;
            p[i][j]=99719;
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            U[i][j][0]=rho[i][j];
            U[i][j][1]=rho[i][j]*u[i][j];
            U[i][j][2]=rho[i][j]*v[i][j];
            U[i][j][3]=p[i][j]/(gam-1.0)+0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]);
        }
    }
}

void partialdiff_jacobi(const std::vector<double>& Xc_b,
                           const std::vector<double>& Yc_b,
                           double k,double l,double ymax,
                        std::vector<std::vector<double>>&px_x,
                        std::vector<std::vector<double>>&px_y,
                        std::vector<std::vector<double>>&py_x,
                        std::vector<std::vector<double>>&py_y,
                        std::vector<std::vector<double>>&x_px,
                        std::vector<std::vector<double>>&x_py,
                        std::vector<std::vector<double>>&y_px,
                        std::vector<std::vector<double>>&y_py,
                        std::vector<std::vector<double>>&Jinv){
    for(int i=0;i<Xc_b.size();i++){
        for(int j=0;j<Yc_b.size();j++){
            px_x[i][j]=1.0;
            px_y[i][j]=0.0;
            py_x[i][j]=(1.0-Yc_b[j]/ymax)*k;
            py_y[i][j]=1.0 - (k/ymax)*(Xc_b[i]-1.0);

            x_px[i][j]=1.0;
            x_py[i][j]=0.0;
            y_px[i][j]=-py_x[i][j]/py_y[i][j];
            y_py[i][j]=-1.0/py_y[i][j];
        }
    }
    int Nl=Xc_b.size();
    for(int i=Nl-1;i>=0;i--){
        if(Xc_b[i]<=1.0){
            Nl=i;
            break;
        }
    }
    for(int i=0;i<=Nl;i++){
        for(int j=0;j<Yc_b.size();j++){
            px_x[i][j]=1.0;
            px_y[i][j]=0.0;
            py_x[i][j]=0.0;
            py_y[i][j]=1.0;
            x_px[i][j]=1.0;
            x_py[i][j]=0.0;
            y_px[i][j]=0.0;
            y_py[i][j]=1.0;
        }
    }

    for(int i=0;i<Xc_b.size();i++){
        for(int j=0;j<Yc_b.size();j++){
            Jinv[i][j]=y_py[i][j];
        }
    }
}

double  MaxNormEigenvalues_Euler2d(const std::vector<std::vector<std::vector<double>>>& U,
                                   const std::vector<std::vector<double>>& Jinv,
                                   const std::vector<std::vector<double>>& x_px,
                                   const std::vector<std::vector<double>>& x_py,
                                   const std::vector<std::vector<double>>& y_px,
                                   const std::vector<std::vector<double>>& y_py,
                                   double gam, char direction){

    int Nx=U.size();
    int Ny=U[0].size();
    std::vector<std::vector<std::vector<double>>> U_p(Nx,std::vector<std::vector<double>>(Ny,std::vector<double>(4,0.0)));
    std::vector<std::vector<double>> rho(Nx,std::vector<double>(Ny,0.0)); 
    std::vector<std::vector<double>> m(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> w(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> E(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> u(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> v(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> p(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> c(Nx,std::vector<double>(Ny,0.0));


    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int m=0;m<4;m++){
                U_p[i][j][m] = U[i][j][m]/Jinv[i+3][j+3];
            }
            rho[i][j]=U_p[i][j][0];
            m[i][j]=U_p[i][j][1];
            w[i][j]=U_p[i][j][2];
            E[i][j]=U_p[i][j][3];

            u[i][j]=m[i][j]/rho[i][j];
            v[i][j]=w[i][j]/rho[i][j];
            p[i][j]=(gam-1.0)*(E[i][j]-0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
            c[i][j]=sqrt(gam*p[i][j]/rho[i][j]);
        }
    }


    double lambda=0.0;

    if(direction=='x'){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                const double alpha=x_px[i][j]*Jinv[i][j];
                const double beta=x_py[i][j]*Jinv[i][j];
                const double eigenvalue1=fabs(alpha*u[i][j]+beta*v[i][j]+c[i][j]*sqrt(alpha*alpha+beta*beta));
                const double eigenvalue2=fabs(alpha*u[i][j]+beta*v[i][j]-c[i][j]*sqrt(alpha*alpha+beta*beta));
                lambda =fmax(lambda,fmax(eigenvalue1,eigenvalue2));
                }
            }
        }

    if(direction=='y'){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                const double alpha=y_px[i][j]*Jinv[i][j];
                const double beta=y_py[i][j]*Jinv[i][j];
                const double eigenvalue1=fabs(alpha*u[i][j]+beta*v[i][j]+c[i][j]*sqrt(alpha*alpha+beta*beta));
                const double eigenvalue2=fabs(alpha*u[i][j]+beta*v[i][j]-c[i][j]*sqrt(alpha*alpha+beta*beta));
                lambda =fmax(lambda,fmax(eigenvalue1,eigenvalue2));
                }
            }
        }
    return lambda; 
}

                                   