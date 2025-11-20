// src/math_utils.h
//Including four functions:   physicalMesh,   uInitialFunc, 
//  partialdiff_jacobi,   MaxNormEigenvalues_Euler2d
#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <vector>


// 函数声明
void physicalMesh(const std::vector<double>& Xc, 
                  const std::vector<double>& Yc,
                  double k, double ymax,double l,
                  std::vector<std::vector<double>>& Xc_p,
                  std::vector<std::vector<double>>& Yc_p);

void uInitialFunc(const std::vector<std::vector<double>>& Xc_p,
                  const std::vector<std::vector<double>>& Yc_p,
                  std::vector<std::vector<std::vector<double>>>& U,
                  double gam);

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
                        std::vector<std::vector<double>>&Jinv);

double MaxNormEigenvalues_Euler2d(const std::vector<std::vector<std::vector<double>>>& U,
                                   const std::vector<std::vector<double>>& Jinv,
                                   const std::vector<std::vector<double>>& x_px,
                                   const std::vector<std::vector<double>>& x_py,
                                   const std::vector<std::vector<double>>& y_px,
                                   const std::vector<std::vector<double>>& y_py,
                                   double gam, char direction);
#endif // MATH_UTILS_H