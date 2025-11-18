#ifndef FLUXFUNCTIONS_EULER2D_H
#define FLUXFUNCTIONS_EULER2D_H
#include <vector>

void FluxFunctions_Euler2d(const std::vector<std::vector<std::vector<double>>>& U,
                    std::vector<std::vector<std::vector<double>>>& F,
                    std::vector<std::vector<std::vector<double>>>& G,
                    double gam);


#endif // FLUXFUNCTIONS_EULER2D_H