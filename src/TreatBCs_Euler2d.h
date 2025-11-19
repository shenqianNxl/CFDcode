#ifndef TREATBCS_EULER2D_H
#define TREATBCS_EULER2D_H

#include <vector>
#include <string>
#include <cmath>

void TreatBCs_Euler2d(std::vector<std::vector<std::vector<double>>>& U_p,
                        std::vector<std::vector<std::vector<double>>>& Ub_ax,
                        std::vector<std::vector<std::vector<double>>>& Ub_ay);

void TreatsBCsLeftRight_Euler(std::vector<std::vector<std::vector<double>>>& Ub_ax,
                              const std::vector<std::vector<std::vector<double>>>& U,
                              const std::string& lefttype,
                              const std::string& righttype);
                              
void TreatsBCsBottomUpper_Euler(std::vector<std::vector<std::vector<double>>>& Ub_ay,
                               const std::vector<std::vector<std::vector<double>>>& U,
                               const std::string& bottomtype,
                               const std::string& uppertype);

#endif // TREATBCS_EULER2D_H




