#include "TreatBCs_Euler2d.h"


void TreatBCs_Euler2d(std::vector<std::vector<std::vector<double>>>& U,
                        std::vector<std::vector<std::vector<double>>>& Ub_ax,
                        std::vector<std::vector<std::vector<double>>>& Ub_ay){
    int Nx=U.size();
    int Ny=U[0].size();

    std::string lefttype = "Dirichlet";
    std::string righttype = "supersonicout";
    std::string bottomtype = "Reflect";
    std::string uppertype = "subsonicout";
    TreatsBCsLeftRight_Euler(Ub_ax,U,lefttype,righttype);
    TreatsBCsBottomUpper_Euler(Ub_ay,U,bottomtype,uppertype);
}

void TreatsBCsLeftRight_Euler(std::vector<std::vector<std::vector<double>>>& Ub_ax,
                              const std::vector<std::vector<std::vector<double>>>& U,
                              const std::string& lefttype,
                              const std::string& righttype){
    //这里的Nx和Ny是计算网格的尺寸，不包含ghost cell
    int Nx=U.size();
    int Ny=U[0].size();
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            Ub_ax[i+3][j]=U[i][j];
        }
    }

    //左边界
    if(lefttype=="Dirichlet"){
        for(int i=0;i<3;i++){
            for(int j=0;j<Ny;j++){
                Ub_ax[i][j][0]=99719/(287.14*293.15); //rho
                Ub_ax[i][j][1]= 99719/(287.14*293.15)*686.47; //rho*u
                Ub_ax[i][j][2]=0.0; //rho*v
                Ub_ax[i][j][3]=99719/0.4 + 0.5*99719/(287.14*293.15)*686.47*686.47;//E

            }
            
        }
    }

    //右边界
    if(righttype=="supersonicout"){
        for(int i=0;i<3;i++){
            for(int j=0;j<Ny;j++){
                Ub_ax[Nx+3+i][j]=U[Nx-1][j];
            }
        }
    }
}

void TreatsBCsBottomUpper_Euler(std::vector<std::vector<std::vector<double>>>& Ub_ay,
                               const std::vector<std::vector<std::vector<double>>>& U,
                               const std::string& bottomtype,
                               const std::string& uppertype){
    int Nx=U.size();
    int Ny=U[0].size();
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            Ub_ay[i][j+3]=U[i][j];
        }
    }


    //下边界,折角前后分情况处理
    int bend_index=Nx/4; //假设折角位置在y方向的1/3处
    if(bottomtype=="Reflect"){
        //折角前部分
        for(int i=0;i<bend_index;i++){
            for(int j=0;j<3;j++){
                Ub_ay[i][j][0]=U[i][2-j][0]; //rho
                Ub_ay[i][j][1]=U[i][2-j][1]; //rho*u
                Ub_ay[i][j][2]=-U[i][2-j][2]; //rho*v
                Ub_ay[i][j][3]=U[i][2-j][3];//E
            }
        }
        //折角后部分
        for(int i=bend_index;i<Nx;i++){
            for(int j=0;j<3;j++){
                Ub_ay[i][j][0]=U[i][2-j][0]; //rho
                Ub_ay[i][j][1]=(sqrt(3.0)/2)*U[i][2-j][1] +0.5*U[i][2-j][2]; //rho*u
                Ub_ay[i][j][2]=0.5*U[i][2-j][1] -(sqrt(3.0)/2)*U[i][2-j][2]; //rho*v
                Ub_ay[i][j][3]=U[i][2-j][3];//E
            }
        }
    }
    //上边界
    if(uppertype=="subsonicout"){
        for(int i=0;i<Nx;i++){
            double rho_i=U[i][Ny-1][0];
            double m_i=U[i][Ny-1][1];
            double w_i=U[i][Ny-1][2];
            double E_i=U[i][Ny-1][3];

            double u_i=m_i/rho_i;
            double v_i=w_i/rho_i;
            double P_i=(0.4)*(E_i -0.5*rho_i*(u_i*u_i+v_i*v_i));
            double c_i=sqrt(1.4*P_i/rho_i);
            double s_i=P_i/(pow(rho_i,1.4));

            double rho_inf=99719/(287.14*293.15);
            double u_inf=686.47;
            double v_inf=0.0;
            double P_inf=99719;
            double c_inf=sqrt(1.4*P_inf/rho_inf);

            double v_e=(v_inf + v_i)/2 + (c_i - c_inf)/2.0/0.4;
            double u_e=u_i;
            double c_e=(1.4-1.0)*0.25*(v_i-v_inf)+0.5*(c_i+c_inf);
            double rho_e=pow(c_e*c_e/ (1.4*s_i), 2.5);
            double P_e=rho_e*c_e*c_e/1.4;
            for(int j=0;j<3;j++){
                Ub_ay[i][Ny+3+j][0]=rho_e; //rho
                Ub_ay[i][Ny+3+j][1]=rho_e*u_e; //rho*u
                Ub_ay[i][Ny+3+j][2]=rho_e*v_e; //rho*v
                Ub_ay[i][Ny+3+j][3]=P_e/0.4 + 0.5*rho_e*(u_e*u_e+v_e*v_e);//E
            }
        }
    }
}