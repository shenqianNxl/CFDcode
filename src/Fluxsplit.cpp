//目前仅仅实现全局L-F格式，后续格式待添加
#include "Fluxsplit.h"

void Fluxsplit_x(//输入
                 const std::vector<std::vector<std::vector<double>>>& Ub_ax_p,
                 const std::vector<std::vector<std::vector<double>>>& F_p_x,
                 const std::vector<std::vector<std::vector<double>>>& G_p_x,
                 //输出
                 std::vector<std::vector<std::vector<double>>>& Fpos_ax,
                 std::vector<std::vector<std::vector<double>>>& Fneg_ax,
                 //其他参数
                 const std::vector<std::vector<double>>& Jinv,
                 const std::vector<std::vector<double>>& x_px,
                 const std::vector<std::vector<double>>& x_py,
                 const std::string& method_splitflux){
    double gam=1.4; //暂时写死，后续可以改成参数传递
    int Nx=Ub_ax_p.size()-6;
    int Ny=Ub_ax_p[0].size();
    std::vector<std::vector<double>> alpha(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> beta(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> theta(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> rho(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> m(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> w(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> E(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> u(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> v(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> P(Nx+6, std:: vector<double>(Ny,0.0));
    std::vector<std::vector<double>> c(Nx+6, std:: vector<double>(Ny,0.0));
    
    //预处理，计算必要的参数
    for(int i=0;i<Nx+6;i++){
        for(int j=0;j<Ny;j++){
            alpha[i][j]=Jinv[i][j+3]*x_px[i][j+3];
            beta[i][j]=Jinv[i][j+3]*x_py[i][j+3];
            theta[i][j]=sqrt(alpha[i][j]*alpha[i][j]+beta[i][j]*beta[i][j]);

            rho[i][j]=Ub_ax_p[i][j][0];
            m[i][j]=Ub_ax_p[i][j][1];
            w[i][j]=Ub_ax_p[i][j][2];
            E[i][j]=Ub_ax_p[i][j][3];
            u[i][j]=m[i][j]/rho[i][j];
            v[i][j]=w[i][j]/rho[i][j];
            P[i][j]=(gam-1)*(E[i][j]-0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
            c[i][j]=sqrt(gam*P[i][j]/rho[i][j]);
        }
    }

    //不同的通量分裂格式
    if(method_splitflux=="LF"){
        double lambda_max=0.0;
        //计算最大特征值
        for(int i=0;i<Nx+6;i++){
            for(int j=0;j<Ny;j++){
                lambda_max=std::max(lambda_max, std::max(fabs(alpha[i][j]*u[i][j]+beta[i][j]*v[i][j]+theta[i][j]*c[i][j]), 
                                            fabs(alpha[i][j]*u[i][j]+beta[i][j]*v[i][j]-theta[i][j]*c[i][j])));
            }
        }

            //计算正负通量
        for(int i=0;i<Nx+6;i++){
            for(int j=0;j<Ny;j++){
                for(int k=0;k<4;k++){
                    double F=Jinv[i][j+3]*(x_px[i][j+3]*F_p_x[i][j][k]+x_py[i][j+3]*G_p_x[i][j][k]);
                    Fpos_ax[i][j][k]=0.5*(F+lambda_max*Ub_ax_p[i][j][k]);
                    Fneg_ax[i][j][k]=0.5*(F -lambda_max*Ub_ax_p[i][j][k]);
                }
            }             
        }
    }
}

void Fluxsplit_y(//输入
                 const std::vector<std::vector<std::vector<double>>>& Ub_ay_p,
                 const std::vector<std::vector<std::vector<double>>>& F_p_y,
                 const std::vector<std::vector<std::vector<double>>>& G_p_y,
                 //输出
                 std::vector<std::vector<std::vector<double>>>& Gpos_ay,
                 std::vector<std::vector<std::vector<double>>>& Gneg_ay,
                 //其他参数
                 const std::vector<std::vector<double>>& Jinv,
                 const std::vector<std::vector<double>>& y_px,
                 const std::vector<std::vector<double>>& y_py,
                 const std::string& method_splitflux){
                double gam=1.4; //暂时写死，后续可以改成参数传递
    int Nx=Ub_ay_p.size();
    int Ny=Ub_ay_p[0].size()-6;
    std::vector<std::vector<double>> alpha(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> beta(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> theta(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> rho(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> m(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> w(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> E(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> u(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> v(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> P(Nx, std:: vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>> c(Nx, std:: vector<double>(Ny+6,0.0));
    
    //预处理，计算必要的参数
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny+6;j++){
            alpha[i][j]=Jinv[i+3][j]*y_px[i+3][j];
            beta[i][j]=Jinv[i+3][j]*y_py[i+3][j];
            theta[i][j]=sqrt(alpha[i][j]*alpha[i][j]+beta[i][j]*beta[i][j]);

            rho[i][j]=Ub_ay_p[i][j][0];
            m[i][j]=Ub_ay_p[i][j][1];
            w[i][j]=Ub_ay_p[i][j][2];
            E[i][j]=Ub_ay_p[i][j][3];
            u[i][j]=m[i][j]/rho[i][j];
            v[i][j]=w[i][j]/rho[i][j];
            P[i][j]=(gam-1)*(E[i][j]-0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
            c[i][j]=sqrt(gam*P[i][j]/rho[i][j]);
        }
    }

    //不同的通量分裂格式
    if(method_splitflux=="LF"){
        double lambda_max=0.0;
        //计算最大特征值
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny+6;j++){
                lambda_max=std::max(lambda_max, std::max(fabs(alpha[i][j]*u[i][j]+beta[i][j]*v[i][j]+theta[i][j]*c[i][j]), 
                                            fabs(alpha[i][j]*u[i][j]+beta[i][j]*v[i][j]-theta[i][j]*c[i][j])));
            }
        }

        //计算正负通量
        for(int i=0;i<Nx;i++){
            for(int j=0;j<Ny+6;j++){
                for(int k=0;k<4;k++){
                    double G=Jinv[i+3][j]*(y_px[i+3][j]*F_p_y[i][j][k]+y_py[i+3][j]*G_p_y[i][j][k]);
                    Gpos_ay[i][j][k]=0.5*(G+lambda_max*Ub_ay_p[i][j][k]);
                    Gneg_ay[i][j][k]=0.5*(G -lambda_max*Ub_ay_p[i][j][k]);
                }
            }
        }
    }
}