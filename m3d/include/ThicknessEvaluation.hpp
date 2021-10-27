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
    void setMesh(const Mesh * _mesh);
    void evalThickness();
    void writeThickness(std::string filename);
    inline const std::vector<double> & thickness() const {return _thickness;};
    inline const unsigned char & algorithm() const {return _algorithm;};
    


  
  private:
    void evalThicknessMethodMartin();
    void evalThicknessAlternativeMethod();
    bool isInElement(const std::vector<double> & xc , const size_t & iElem);
    GrahmOperatorOutput GrahmOperations( const Point & p0, const Point & p1, const Point & p2, const std::vector<double> & x_c);
    double distanceOfPointToPlane(const std::vector<double> & N, const Point & p0,  const std::vector<double> & x_c, const std::vector<double> & mu);
    double normalDistanceOfPointToPlane(const std::vector<double> & N, const Point & p0,  const std::vector<double> & x_c);
    double pointDistances( const Point & p0, const Point & p1);
    std::vector<double> waxpy(const std::vector<double> & x, const std::vector<double> & y, double alpha=1.0);
    std::vector<double> _thickness;
    unsigned char _algorithm;
    bool _swapregions;


};









#endif
