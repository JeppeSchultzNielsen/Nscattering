//
// Created by jeppe on 3/10/22.
//
#include <string>
#include <vector>
using namespace std;

#ifndef ALUSCATTERING_RUNNER_H
#define ALUSCATTERING_RUNNER_H

void createFile(string in);
std::vector<double> thickness(string in);
double gaussSum(double *x, double *par);
std::vector<double> findCurrent(string in);

#endif //ALUSCATTERING_RUNNER_H
