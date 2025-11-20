// src/RK3_WENO_Euler2d.h
#ifndef RK3_WENO_EULER2D_H
#define RK3_WENO_EULER2D_H

#include <vector>
#include <string>

void RK3_WENO_Euler2d(std::vector<std::vector<std::vector<double>>>& Unew,
                      double dx, double dy, double dt,
                      const std::vector<std::vector<std::vector<double>>>& U,
                      const std::vector<std::vector<double>>& x_px,
                      const std::vector<std::vector<double>>& x_py,
                      const std::vector<std::vector<double>>& y_px,
                      const std::vector<std::vector<double>>& y_py,
                      const std::vector<std::vector<double>>& Jinv,
                      double gam,
                      const std::string& method_splitflux,
                      const std::string& method_WENO);

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
                             const std::string& method_WENO);
#endif // RK3_WENO_EULER2D_H