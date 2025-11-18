#ifndef FLUXSPILT_H
#define FLUXSPILT_H

#include <vector>
#include <string>
#include <cmath>

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
                 const std::string& method_splitflux);

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
                 const std::string& method_splitflux);

#endif // FLUXSPILT_H