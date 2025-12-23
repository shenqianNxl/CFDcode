// src/math_utils.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include "math_utils.h"

//调试代码
// 写二维矩阵到CSV（行: j, 列: i）
void write2DCSV(const std::string& filename,
                       const std::vector<std::vector<double>>& A) {
    std::ofstream f(filename);
    f.setf(std::ios::scientific);
    f << std::setprecision(16);
    const size_t Nx = A.size();
    const size_t Ny = Nx ? A[0].size() : 0;
    for (size_t j = 0; j < Ny; ++j) {
        for (size_t i = 0; i < Nx; ++i) {
            f << A[i][j];
            if (i + 1 < Nx) f << ",";
        }
        f << "\n";
    }
}
// 新增：将 U 或 U_p 的四个分量分别输出为二维网格 CSV
void dumpUComponents(const std::string& outDir,
                            const std::string& prefix,
                            const std::vector<std::vector<std::vector<double>>>& U3) {
    using std::vector;
    const size_t Nx = U3.size();
    const size_t Ny = Nx ? U3[0].size() : 0;
    vector<vector<double>> c0(Nx, vector<double>(Ny));
    vector<vector<double>> c1(Nx, vector<double>(Ny));
    vector<vector<double>> c2(Nx, vector<double>(Ny));
    vector<vector<double>> c3(Nx, vector<double>(Ny));
    for (size_t i = 0; i < Nx; ++i) {
        for (size_t j = 0; j < Ny; ++j) {
            c0[i][j] = U3[i][j][0]; // rho
            c1[i][j] = U3[i][j][1]; // rho*u
            c2[i][j] = U3[i][j][2]; // rho*v
            c3[i][j] = U3[i][j][3]; // E
        }
    }
    std::filesystem::create_directories(outDir);
    write2DCSV(outDir + "/" + prefix + "_rho.csv", c0);
    write2DCSV(outDir + "/" + prefix + "_rhou.csv", c1);
    write2DCSV(outDir + "/" + prefix + "_rhov.csv", c2);
    write2DCSV(outDir + "/" + prefix + "_E.csv",   c3);
}

void physicalMesh(const std::vector<double>& Xc, 
                  const std::vector<double>& Yc,
                  double k, double ymax,double l,
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
            
            if (Xc[i] < l) {
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
    //Xc_b.size()=Nx+6,Yc_b.size()=Ny+6
    for(size_t i=0;i<Xc_b.size();i++){
        for(size_t j=0;j<Yc_b.size();j++){
            px_x[i][j]=1.0;
            px_y[i][j]=0.0;
            py_x[i][j]=(1.0-Yc_b[j]/ymax)*k;
            py_y[i][j]=1.0 - (k/ymax)*(Xc_b[i]-1.0);

            x_px[i][j]=1.0;
            x_py[i][j]=0.0;
            y_px[i][j]=-py_x[i][j]/py_y[i][j];
            y_py[i][j]=1.0/py_y[i][j];
        }
    }
    int Nl=Xc_b.size();
    for(int i=Nl-1;i>=0;i--){
        if(Xc_b[i]<=l){
            Nl=i;
            break;
        }
    }
    for(int i=0;i<=Nl;i++){
        for(size_t j=0;j<Yc_b.size();j++){
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

    for(size_t i=0;i<Xc_b.size();i++){
        for(size_t j=0;j<Yc_b.size();j++){
            Jinv[i][j]=py_y[i][j];
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
    /*U的尺寸是Nx*Ny*4*/
    int Nx=U.size();
    int Ny=U[0].size();
    std::vector<std::vector<std::vector<double>>> U_p(Nx,std::vector<std::vector<double>>(Ny,std::vector<double>(4,0.0)));
    std::vector<std::vector<double>> rho(Nx,std::vector<double>(Ny,0.0)); 
    std::vector<std::vector<double>> m(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> w(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> E(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> u(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> v(Nx,std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> P(Nx,std::vector<double>(Ny,0.0));
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
            P[i][j]=(gam-1.0)*(E[i][j]-0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
            c[i][j]=sqrt(gam*P[i][j]/rho[i][j]);
        }
    }


    double lambda=0.0;

    if(direction=='x'){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                const double alpha=x_px[i+3][j+3]*Jinv[i+3][j+3];
                const double beta=x_py[i+3][j+3]*Jinv[i+3][j+3];
                const double eigenvalue1=fabs(alpha*u[i][j]+beta*v[i][j]+c[i][j]*sqrt(alpha*alpha+beta*beta));
                const double eigenvalue2=fabs(alpha*u[i][j]+beta*v[i][j]-c[i][j]*sqrt(alpha*alpha+beta*beta));
                lambda =fmax(lambda,fmax(eigenvalue1,eigenvalue2));
                }
            }
        }

    if(direction=='y'){
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny;j++){
                const double alpha=y_px[i+3][j+3]*Jinv[i+3][j+3];
                const double beta=y_py[i+3][j+3]*Jinv[i+3][j+3];
                const double eigenvalue1=fabs(alpha*u[i][j]+beta*v[i][j]+c[i][j]*sqrt(alpha*alpha+beta*beta));
                const double eigenvalue2=fabs(alpha*u[i][j]+beta*v[i][j]-c[i][j]*sqrt(alpha*alpha+beta*beta));
                lambda =fmax(lambda,fmax(eigenvalue1,eigenvalue2));
                }
            }
        }
    return lambda; 
}

                                   