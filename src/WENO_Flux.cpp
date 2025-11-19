#include "WENO_Flux.h"
void WENO_FluxPos1d(//输入
                    const std::vector<std::vector<std::vector<double>>> &Fpos_ax,
                    const std::string& method_WENO,
                    const std::string& recon_direction,
                    //输出
                    std::vector<std::vector<std::vector<double>>> &Fhat_pos){
    if(method_WENO=="WENO-JS"){
        double rec_d0=0.1;
        double rec_d1=0.6;
        double rec_d2=0.3;
        double epsilon=1.0e-10;

        int Nx=Fpos_ax.size();
        int Ny=Fpos_ax[0].size();

        if(recon_direction=="x"){
            for(int i=0;i<Nx-5;i++){
                for(int j=0;j<Ny;j++){
                    for(int k=0;k<4;k++){
                        ///计算三个子模板的光滑指标beta0,beta1,beta2
                    double beta0=13.0/12.0*(Fpos_ax[i][j][k]-2.0*Fpos_ax[i+1][j][k]+Fpos_ax[i+2][j][k])*
                                 (Fpos_ax[i][j][k]-2.0*Fpos_ax[i+1][j][k]+Fpos_ax[i+2][j][k])+
                                 1.0/4.0*(Fpos_ax[i][j][k]-4.0*Fpos_ax[i+1][j][k]+3.0*Fpos_ax[i+2][j][k])*
                                 (Fpos_ax[i][j][k]-4.0*Fpos_ax[i+1][j][k]+3.0*Fpos_ax[i+2][j][k]);

                    double beta1=13.0/12.0*(Fpos_ax[i+1][j][k]-2.0*Fpos_ax[i+2][j][k]+Fpos_ax[i+3][j][k])*
                                 (Fpos_ax[i+1][j][k]-2.0*Fpos_ax[i+2][j][k]+Fpos_ax[i+3][j][k])+
                                 1.0/4.0*(Fpos_ax[i+1][j][k]-Fpos_ax[i+3][j][k])*
                                 (Fpos_ax[i+1][j][k]-Fpos_ax[i+3][j][k]);

                    double beta2=13.0/12.0*(Fpos_ax[i+2][j][k]-2.0*Fpos_ax[i+3][j][k]+Fpos_ax[i+4][j][k])*
                                 (Fpos_ax[i+2][j][k]-2.0*Fpos_ax[i+3][j][k]+Fpos_ax[i+4][j][k])+
                                 1.0/4.0*(3.0*Fpos_ax[i+2][j][k]-4.0*Fpos_ax[i+3][j][k]+Fpos_ax[i+4][j][k])*
                                 (3.0*Fpos_ax[i+2][j][k]-4.0*Fpos_ax[i+3][j][k]+Fpos_ax[i+4][j][k]);

                    //计算非线性权重
                    double alpha0=rec_d0/( (epsilon+beta0)*(epsilon+beta0) );
                    double alpha1=rec_d1/( (epsilon+beta1)*(epsilon+beta1) );
                    double alpha2=rec_d2/( (epsilon+beta2)*(epsilon+beta2) );
                    double alpha_sum=alpha0+alpha1+alpha2;
                    double w0=alpha0/alpha_sum;
                    double w1=alpha1/alpha_sum;
                    double w2=alpha2/alpha_sum;
                    //计算重构值
                    Fhat_pos[i][j][k]=w0*( (2.0*Fpos_ax[i][j][k]-7.0*Fpos_ax[i+1][j][k]+11.0*Fpos_ax[i+2][j][k])/6.0 )+
                                 w1*( (-Fpos_ax[i+1][j][k]+5.0*Fpos_ax[i+2][j][k]+2.0*Fpos_ax[i+3][j][k])/6.0 )+
                                 w2*( (2.0*Fpos_ax[i+2][j][k]+5.0*Fpos_ax[i+3][j][k]-Fpos_ax[i+4][j][k])/6.0 );
                    }
                }
            }   
        }
        else if(recon_direction=="y"){
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny-5;j++){
                    for(int k=0;k<4;k++){
                    ///计算三个子模板的光滑指标beta0,beta1,beta2
                    double beta0=13.0/12.0*(Fpos_ax[i][j][k]-2.0*Fpos_ax[i][j+1][k]+Fpos_ax[i][j+2][k])*
                                 (Fpos_ax[i][j][k]-2.0*Fpos_ax[i][j+1][k]+Fpos_ax[i][j+2][k])+
                                 1.0/4.0*(Fpos_ax[i][j][k]-4.0*Fpos_ax[i][j+1][k]+3.0*Fpos_ax[i][j+2][k])*
                                 (Fpos_ax[i][j][k]-4.0*Fpos_ax[i][j+1][k]+3.0*Fpos_ax[i][j+2][k]);

                    double beta1=13.0/12.0*(Fpos_ax[i][j+1][k]-2.0*Fpos_ax[i][j+2][k]+Fpos_ax[i][j+3][k])*
                                 (Fpos_ax[i][j+1][k]-2.0*Fpos_ax[i][j+2][k]+Fpos_ax[i][j+3][k])+
                                    1.0/4.0*(Fpos_ax[i][j+1][k]-Fpos_ax[i][j+3][k])*
                                    (Fpos_ax[i][j+1][k]-Fpos_ax[i][j+3][k]);

                    double beta2=13.0/12.0*(Fpos_ax[i][j+2][k]-2.0*Fpos_ax[i][j+3][k]+Fpos_ax[i][j+4][k])*
                                 (Fpos_ax[i][j+2][k]-2.0*Fpos_ax[i][j+3][k]+Fpos_ax[i][j+4][k])+
                                 1.0/4.0*(3.0*Fpos_ax[i][j+2][k]-4.0*Fpos_ax[i][j+3][k]+Fpos_ax[i][j+4][k])*
                                 (3.0*Fpos_ax[i][j+2][k]-4.0*Fpos_ax[i][j+3][k]+Fpos_ax[i][j+4][k]);

                    //计算非线性权重
                    double alpha0=rec_d0/( (epsilon+beta0)*(epsilon+beta0) );
                    double alpha1=rec_d1/( (epsilon+beta1)*(epsilon+beta1) );
                    double alpha2=rec_d2/( (epsilon+beta2)*(epsilon+beta2) );
                    double alpha_sum=alpha0+alpha1+alpha2;
                    double w0=alpha0/alpha_sum;
                    double w1=alpha1/alpha_sum;
                    double w2=alpha2/alpha_sum;
                    //计算重构值
                    Fhat_pos[i][j][k]=w0*( (2.0*Fpos_ax[i][j][k]-7.0*Fpos_ax[i][j+1][k]+11.0*Fpos_ax[i][j+2][k])/6.0 )+
                                 w1*( (-Fpos_ax[i][j+1][k]+5.0*Fpos_ax[i][j+2][k]+2.0*Fpos_ax[i][j+3][k])/6.0 )+
                                 w2*( (2.0*Fpos_ax[i][j+2][k]+5.0*Fpos_ax[i][j+3][k]-Fpos_ax[i][j+4][k])/6.0 );
                                
                    }
                }
            }   
        }
        else{
            std::cout<<"Invalid recon_direction in WENO_FluxPos1d. Use 'x' or 'y'."<<std::endl;
        }
    }
}


void WENO_FluxNeg1d(//输入
                    const std::vector<std::vector<std::vector<double>>> &Fneg_ax,
                    const std::string& method_WENO,
                    const std::string& recon_direction,
                    //输出
                    std::vector<std::vector<std::vector<double>>> &Fhat_neg){
                        if(method_WENO=="WENO-JS"){
        double rec_d0=0.1;
        double rec_d1=0.6;
        double rec_d2=0.3;
        double epsilon=1.0e-10;

        int Nx=Fneg_ax.size();
        int Ny=Fneg_ax[0].size();

        if(recon_direction=="x"){
            for(int i=0;i<Nx-5;i++){
                for(int j=0;j<Ny;j++){
                    for(int k=0;k<4;k++){
                        ///计算三个子模板的光滑指标beta0,beta1,beta2
                    double beta0=13.0/12.0*(Fneg_ax[i+5][j][k]-2.0*Fneg_ax[i+4][j][k]+Fneg_ax[i+3][j][k])*
                                 (Fneg_ax[i+5][j][k]-2.0*Fneg_ax[i+4][j][k]+Fneg_ax[i+3][j][k])+
                                 1.0/4.0*(Fneg_ax[i+5][j][k]-4.0*Fneg_ax[i+4][j][k]+3.0*Fneg_ax[i+3][j][k])*
                                 (Fneg_ax[i+5][j][k]-4.0*Fneg_ax[i+4][j][k]+3.0*Fneg_ax[i+3][j][k]);

                    double beta1=13.0/12.0*(Fneg_ax[i+4][j][k]-2.0*Fneg_ax[i+3][j][k]+Fneg_ax[i+2][j][k])*
                                 (Fneg_ax[i+4][j][k]-2.0*Fneg_ax[i+3][j][k]+Fneg_ax[i+2][j][k])+
                                 1.0/4.0*(Fneg_ax[i+4][j][k]-Fneg_ax[i+2][j][k])*
                                 (Fneg_ax[i+4][j][k]-Fneg_ax[i+2][j][k]);

                    double beta2=13.0/12.0*(Fneg_ax[i+3][j][k]-2.0*Fneg_ax[i+2][j][k]+Fneg_ax[i+1][j][k])*
                                 (Fneg_ax[i+3][j][k]-2.0*Fneg_ax[i+2][j][k]+Fneg_ax[i+1][j][k])+
                                 1.0/4.0*(3.0*Fneg_ax[i+3][j][k]-4.0*Fneg_ax[i+2][j][k]+Fneg_ax[i+1][j][k])*
                                 (3.0*Fneg_ax[i+3][j][k]-4.0*Fneg_ax[i+2][j][k]+Fneg_ax[i+1][j][k]);

                    //计算非线性权重
                    double alpha0=rec_d0/( (epsilon+beta0)*(epsilon+beta0) );
                    double alpha1=rec_d1/( (epsilon+beta1)*(epsilon+beta1) );
                    double alpha2=rec_d2/( (epsilon+beta2)*(epsilon+beta2) );
                    double alpha_sum=alpha0+alpha1+alpha2;
                    double w0=alpha0/alpha_sum;
                    double w1=alpha1/alpha_sum;
                    double w2=alpha2/alpha_sum;
                    //计算重构值
                    Fhat_neg[i][j][k]=w0*( (2.0*Fneg_ax[i+5][j][k]-7.0*Fneg_ax[i+4][j][k]+11.0*Fneg_ax[i+3][j][k])/6.0 )+
                                 w1*( (-Fneg_ax[i+4][j][k]+5.0*Fneg_ax[i+3][j][k]+2.0*Fneg_ax[i+2][j][k])/6.0 )+
                                 w2*( (2.0*Fneg_ax[i+3][j][k]+5.0*Fneg_ax[i+2][j][k]-Fneg_ax[i+1][j][k])/6.0 );
                    }
                }
            }   
        }
        else if(recon_direction=="y"){
            for(int i=0;i<Nx;i++){
                for(int j=0;j<Ny-5;j++){
                    for(int k=0;k<4;k++){
                    ///计算三个子模板的光滑指标beta0,beta1,beta2
                    double beta0=13.0/12.0*(Fneg_ax[i][j+5][k]-2.0*Fneg_ax[i][j+4][k]+Fneg_ax[i][j+3][k])*
                                 (Fneg_ax[i][j+5][k]-2.0*Fneg_ax[i][j+4][k]+Fneg_ax[i][j+3][k])+
                                 1.0/4.0*(Fneg_ax[i][j+5][k]-4.0*Fneg_ax[i][j+4][k]+3.0*Fneg_ax[i][j+3][k])*
                                 (Fneg_ax[i][j+5][k]-4.0*Fneg_ax[i][j+4][k]+3.0*Fneg_ax[i][j+3][k]);

                    double beta1=13.0/12.0*(Fneg_ax[i][j+4][k]-2.0*Fneg_ax[i][j+3][k]+Fneg_ax[i][j+2][k])*
                                 (Fneg_ax[i][j+4][k]-2.0*Fneg_ax[i][j+3][k]+Fneg_ax[i][j+2][k])+
                                    1.0/4.0*(Fneg_ax[i][j+4][k]-Fneg_ax[i][j+2][k])*
                                    (Fneg_ax[i][j+4][k]-Fneg_ax[i][j+2][k]);

                    double beta2=13.0/12.0*(Fneg_ax[i][j+3][k]-2.0*Fneg_ax[i][j+2][k]+Fneg_ax[i][j+1][k])*
                                 (Fneg_ax[i][j+3][k]-2.0*Fneg_ax[i][j+2][k]+Fneg_ax[i][j+1][k])+
                                 1.0/4.0*(3.0*Fneg_ax[i][j+3][k]-4.0*Fneg_ax[i][j+2][k]+Fneg_ax[i][j+1][k])*
                                 (3.0*Fneg_ax[i][j+3][k]-4.0*Fneg_ax[i][j+2][k]+Fneg_ax[i][j+1][k]);

                    //计算非线性权重
                    double alpha0=rec_d0/( (epsilon+beta0)*(epsilon+beta0) );
                    double alpha1=rec_d1/( (epsilon+beta1)*(epsilon+beta1) );
                    double alpha2=rec_d2/( (epsilon+beta2)*(epsilon+beta2) );
                    double alpha_sum=alpha0+alpha1+alpha2;
                    double w0=alpha0/alpha_sum;
                    double w1=alpha1/alpha_sum;
                    double w2=alpha2/alpha_sum;
                    //计算重构值
                    Fhat_neg[i][j][k]=w0*( (2.0*Fneg_ax[i][j+5][k]-7.0*Fneg_ax[i][j+4][k]+11.0*Fneg_ax[i][j+3][k])/6.0 )+
                                 w1*( (-Fneg_ax[i][j+4][k]+5.0*Fneg_ax[i][j+3][k]+2.0*Fneg_ax[i][j+2][k])/6.0 )+
                                 w2*( (2.0*Fneg_ax[i][j+3][k]+5.0*Fneg_ax[i][j+2][k]-Fneg_ax[i][j+1][k])/6.0 );

                    }
                }
            }   
        }
        else{
            std::cout<<"Invalid recon_direction in WENO_FluxPos1d. Use 'x' or 'y'."<<std::endl;
        }
    }
}