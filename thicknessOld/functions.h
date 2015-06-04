#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

double GetRandomNumber();

void GetScatteringAngles(double mu_x, double mu_y, double mu_z, double &mu_x_new, double &mu_y_new, double &mu_z_new, double g);

int GetCurrentPixel(double x, double y, double z, double x_res, double y_res, double z_res, int sizex, int sizey, int sizez);

double GetReflectionCoefficient(double alpha_i, double n_i, double n_t);

int GetCurrentElement(double x, double y, double z, double x_res, double y_res, double z_res, int sizex, int sizey, int sizez, double** &coords, int** &elems, int* &elementPixels);

double DeterminantCalc(double m00, double m01, double m02, double m03, double m10, double m11, double m12, double m13, double m20, double m21, double m22, double m23, double m30, double m31, double m32, double m33);

double* calculateInwardNormal(int element, int node0, int node1, int node2, double** &coords, double** &cents);

double distanceOfPointToPlane(double* &N, int node, double** &coords, double x, double y, double z, double mu_x, double mu_y, double mu_z);

double normalDistanceOfPointToPlane(double* &N, int node, double** &coords, double x, double y, double z);

int cornerEdgeChecker(double x, double y, double z, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);

int faceChecker(double x, double y, double z, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);

#endif
