#ifndef WENO_FLUX_H
#define WENO_FLUX_H
#include <iostream>
#include <vector>
#include <string>

void WENO_FluxPos1d(//输入
                    const std::vector<std::vector<std::vector<double>>> &Fpos_ax,
                    const std::string& method_WENO,
                    const std::string& recon_direction,
                    //输出
                    std::vector<std::vector<std::vector<double>>> &Fhat_pos);

void WENO_FluxNeg1d(//输入
                    const std::vector<std::vector<std::vector<double>>> &Fneg_ax,
                    const std::string& method_WENO,
                    const std::string& recon_direction,
                    //输出
                    std::vector<std::vector<std::vector<double>>> &Fhat_neg);
#endif // WENO_FLUX_H