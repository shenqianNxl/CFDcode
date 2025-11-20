//考虑单精度和双精度对计算效率的影响
//考虑在求计算网格到物理网格时，右上方ghost cell的处理
#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <string>
#include <fstream>
#include <iomanip>

#include "src/math_utils.h" 
#include "src/RK3_WENO_Euler2d.h"

int main(){

    //输入x、y方向的网格划分数目以及WENO格式和通量分裂方法
    std::cout << "Please input the number of subdivisions for the grid in the x-direction and y-direction." << std::endl;
    int Nx, Ny;
    std::cin >> Nx >> Ny;
    std::cout << "Please specify the desired WENO scheme. (Options: WENO-ZQ, WENO5-JS)" << std::endl;
    std::string method_WENO;
    std::cin >> method_WENO;
    std::cout<<"Please specify the desired flux splitting method.(Options:LF,LFlocal,SW)"<<std::endl;
    std::string method_splitflux;
    std::cin>>method_splitflux;
    //初始参数，另外注意考虑单精度与双精度的问题
    double gam=1.4;
    //以下两个参数可调，可以多试试
    double T=0.05;
    double CFL=0.6;
    //几何形状
    double l=1.0; //折角处x坐标
    double xmin=0.0, xmax=4;
    double ymin=0.0, ymax=2.4;
    const double pi = 3.14159265358979323846;
    double theta= pi/12; //15度折角

    double dx=(xmax-xmin)/Nx;
    double dy=(ymax-ymin)/Ny;

    //生成计算网格和物理网格
    //计算计算网格节点和单元中心坐标
    //Xc,Yc为一维向量，存储计算单元中心坐标，在计算域使用
    std::vector<double>X,Y;
    for(int i=0;i<=Nx;i++){
        X.push_back(xmin+i*dx);
    }
    for(int j=0;j<=Ny;j++){
        Y.push_back(ymin+j*dy);
    }
    std::vector<double> Xc, Yc;
    for(int i=0;i<Nx;i++){
        Xc.push_back(0.5*(X[i]+X[i+1]));
    }
    for(int j=0;j<Ny;j++){
        Yc.push_back(0.5*(Y[j]+Y[j+1]));
    }
    //计算网格变形到物理网格，Xc_p,Yc_p为二维向量，存储物理单元中心坐标
    std::vector<std::vector<double>>Xc_p,Yc_p;
    physicalMesh(Xc,Yc,tan(theta),ymax,l,Xc_p,Yc_p);

    //初值条件
    double t=0;
    std::vector<std::vector<std::vector<double>>> U(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    uInitialFunc(Xc_p,Yc_p,U,gam);

    //计算边界扩展后的边界点坐标
    std::vector<double> Xc_b(Nx+6,0.0); 
    std::vector<double> Yc_b(Ny+6,0.0);
    for(int i=0;i<Nx+6;i++){
        if(i<3){
            Xc_b[i]=Xc[0]-(3-i)*dx;
        }else if(i>=3 && i<=Nx+2){
            Xc_b[i]=Xc[i-3];
        }else{
            Xc_b[i]=Xc[Nx-1]+(i-(Nx+2))*dx;
        }
    }
    for(int j=0;j<Ny+6;j++){
        if(j<3){
            Yc_b[j]=Yc[0]-(3-j)*dy;
        }else if(j>=3 && j<=Ny+2){
            Yc_b[j]=Yc[j-3];
        }else{
            Yc_b[j]=Yc[Ny-1]+(j-(Ny+2))*dy;
        }
    }

    //网格变换下，物理域与计算域各变量之间的偏导数关系，并计算jacobi矩阵的逆的行列式Jinv
    std::vector<std::vector<double>>px_x(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>px_y(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>py_x(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>py_y(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>x_px(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>x_py(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>y_px(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>y_py(Nx+6,std::vector<double>(Ny+6,0.0));
    std::vector<std::vector<double>>Jinv(Nx+6,std::vector<double>(Ny+6,0.0));
    partialdiff_jacobi(Xc_b,Yc_b,tan(theta),l,ymax,px_x,px_y,py_x,py_y,
                        x_px,x_py,y_px,y_py,Jinv);

    //时间推进主循环
    std::time_t start = std::time(nullptr); //记录起始时间
    double lamda1,lamda2;
    double dt;
    while(1){
        //计算当前时间步长
        double dt_min=1e-12;
        lamda1= MaxNormEigenvalues_Euler2d(U,Jinv,x_px,x_py,y_px,y_py,gam,'x');
        lamda2= MaxNormEigenvalues_Euler2d(U,Jinv,x_px,x_py,y_px,y_py,gam,'y');
        dt=CFL/(lamda1/dx+lamda2/dy);
        if(t+dt>T){
            dt=T-t;
        }
        if(dt<dt_min){
            break;
        }

        //使用三阶Runge-Kutta方法进行时间推进
        std::vector<std::vector<std::vector<double>>> U_new(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
        RK3_WENO_Euler2d(U_new,t,dx,dy,dt,U,x_px,x_py,y_px,y_py,Jinv,gam,method_splitflux,method_WENO);
        t+=dt;
        std::cout<<"当前时间t="<<t<<std::endl;
        U=U_new;
    }

    //计算最终的物理量
    std::vector<std::vector<std::vector<double>>> U_p(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(4,0.0)));
    std::vector<std::vector<double>> rho(Nx, std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> u(Nx, std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> v(Nx, std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> E(Nx, std::vector<double>(Ny,0.0));
    std::vector<std::vector<double>> P(Nx, std::vector<double>(Ny,0.0));
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int m=0;m<4;m++){
                U_p[i][j][m]=U[i][j][m]/Jinv[i+3][j+3];
            }
            rho[i][j]=U_p[i][j][0];
            u[i][j]=U_p[i][j][1]/U_p[i][j][0];
            v[i][j]=U_p[i][j][2]/U_p[i][j][0];
            E[i][j]=U_p[i][j][3];
            P[i][j]=(gam-1.0)*(E[i][j]-0.5*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
        }
    }

    //输出参数到CSV文件
    // 设置输出精度
    std::ofstream file("fluid_data.csv");
    file << std::fixed << std::setprecision(6);
    // 写入表头
    file << "X,Y,rho,u,v,P\n";
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            file << Xc_p[i][j] << ","
                 << Yc_p[i][j] << ","
                 << rho[i][j] << ","
                 << u[i][j] << ","
                 << v[i][j] << ","
                 << P[i][j] << "\n";
        }
    }
    file.close();
    std::cout << "Data saved to fluid_data.csv" << std::endl;
    

    std::time_t end = std::time(nullptr); // 记录结束时间
    double duration = std::difftime(end, start); // 计算时间差，单位为秒
    std::cout << "操作耗时: " << duration << " 秒。" << std::endl;

    return 0;
}