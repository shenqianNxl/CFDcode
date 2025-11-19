#include "FluxFunctions_Euler2d.h"

void FluxFunctions_Euler2d(const std::vector<std::vector<std::vector<double>>>& U,
                    std::vector<std::vector<std::vector<double>>>& F,
                    std::vector<std::vector<std::vector<double>>>& G,
                    double gam){
    int Nx=U.size();
    int Ny=U[0].size();

    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            double rho=U[i][j][0];
            double m=U[i][j][1];
            double w=U[i][j][2];
            double E=U[i][j][3];

            double u=m/rho;
            double v=w/rho;
            double P=(gam-1.0)*(E-0.5*rho*(u*u+v*v));

            F[i][j][0]=m;
            F[i][j][1]=rho*u*u + P;
            F[i][j][2]=m*v;
            F[i][j][3]=(E+P)*u;

            G[i][j][0]=w;
            G[i][j][1]=w*u;
            G[i][j][2]=w*v + P;
            G[i][j][3]=(E+P)*v;
        }
    }
}