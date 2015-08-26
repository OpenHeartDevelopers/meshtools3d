#ifndef _THICKEVAL_HPP
#define _THICKEVAL_HPP


#include <fstream>
#include <string>
#include <vector>
#include<map>
#include "Mesh.hpp"
#include "GetPot.hpp"
#include "LaplaceSolver.hpp"



struct GrahmOperatorOutput
{
  std::vector<double> N, projPt;
  double dist;
  GrahmOperatorOutput()
  :N(3,0),
  projPt(3,0),
  dist(0)
  {}
  inline double & Nx(){return N[0];};
  inline double & Ny(){return N[1];};
  inline double & Nz(){return N[3];};
};

class ThicknessEvaluation:
public  LaplaceSolver
{
  public:
    ThicknessEvaluation();
    ThicknessEvaluation(const Mesh * _mesh);
    ThicknessEvaluation(const GetPot & dfile, const Mesh * _mesh);
    ~ThicknessEvaluation();


  
  private:
    void evalThickness();
    bool isInElement(std::vector<double> xc ,size_t iElem);
    std::vector<double> waxpy(std::vector<double> & x, std::vector<double> & y, double alpha=1.0);
    GrahmOperatorOutput GrahmOperations( Point & p0, Point & p1, Point & p2, std::vector<double> & x_c);
    double distanceOfPointToPlane(std::vector<double> & N, Point & p0,  std::vector<double> & x_c, std::vector<double> & mu);
    double normalDistanceOfPointToPlane(std::vector<double> & N, Point & p0,  std::vector<double> & x_c);


};









#endif
