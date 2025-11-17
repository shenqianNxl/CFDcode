#include "RK3_WENO_Euler2d.h"


void RK3_WENO_Euler2d(std::vector<std::vector<std::vector<double>>>& Unew,
                      double t, double dx, double dy, double dt,
                      const std::vector<std::vector<std::vector<double>>>& U,
                      const std::vector<std::vector<double>>& x_px,
                      const std::vector<std::vector<double>>& x_py,
                      const std::vector<std::vector<double>>& y_px,
                      const std::vector<std::vector<double>>& y_py,
                      const std::vector<std::vector<double>>& Jinv,
                      double gam,
                      const std::string& method_splitflux,
                      const std::string& method_WENO){

    int Nx=U.size();
    int Ny=U[0].size();

    std::vector<std::vector<std::vector<double>>> Fhat(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> Ghat(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> F1hat(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> G1hat(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> F2hat(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> G2hat(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> U1(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<std::vector<double>>> U2(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));


    //第一阶段
    Flux_LFsplitBased_Euler2d(U,x_px,x_py,y_px,y_py,Jinv,gam,Fhat,Ghat,method_splitflux,method_WENO);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int m=0;m<4;m++){
                U1[i][j][m]=U[i][j][m]-dt/dx*(Fhat[i+1][j][m]-Fhat[i][j][m])-dt/dy*(Ghat[i][j+1][m]-Ghat[i][j][m]);
            }
        }
    }
    //第二阶段
    Flux_LFsplitBased_Euler2d(U1,x_px,x_py,y_px,y_py,Jinv,gam,F1hat,G1hat,method_splitflux,method_WENO);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int m=0;m<4;m++){
                U2[i][j][m]=0.75*U[i][j][m]+0.25*(U1[i][j][m]-dt/dx*(F1hat[i+1][j][m]-F1hat[i][j][m])-dt/dy*(G1hat[i][j+1][m]-G1hat[i][j][m]));
            }
        }
    }
    //第三阶段
    Flux_LFsplitBased_Euler2d(U2,x_px,x_py,y_px,y_py,Jinv,gam,F2hat,G2hat,method_splitflux,method_WENO);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int m=0;m<4;m++){
                Unew[i][j][m]=(1.0/3.0)*U[i][j][m]+(2.0/3.0)*(U2[i][j][m]-dt/dx*(F2hat[i+1][j][m]-F2hat[i][j][m])-dt/dy*(G2hat[i][j+1][m]-G2hat[i][j][m]));
            }
        }
    }
}

void  Flux_LFsplitBased_Euler2d(const std::vector<std::vector<std::vector<double>>>& U,
                             const std::vector<std::vector<double>>& x_px,
                             const std::vector<std::vector<double>>& x_py,
                             const std::vector<std::vector<double>>& y_px,
                             const std::vector<std::vector<double>>& y_py,
                             const std::vector<std::vector<double>>& Jinv,
                             double gam,
                             std::vector<std::vector<std::vector<double>>>& Fhat,
                             std::vector<std::vector<std::vector<double>>>& Ghat,
                             const std::string& method_splitflux,
                             const std::string& method_WENO){
    int Nx=U.size();
    int Ny=U[0].size();      
    std::vector<std::vector<std::vector<double>>> U_p(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));      

    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int m=0;m<4;m++){
                U_p[i][j][m]=U[i][j][m]/Jinv[i+3][j+3];
            }
        }
    }
}