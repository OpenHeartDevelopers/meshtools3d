#ifndef _THICKEVAL_HPP
#define _THICKEVAL_HPP


#include <fstream>
#include <string>
#include <vector>
#include<map>
#include "Mesh.hpp"
#include "GetPot.hpp"
#include "LaplaceSolver.hpp"

class ThicknessEvaluation:
public  LaplaceSolver
{
  public:
    ThicknessEvaluation();
    ThicknessEvaluation(const Mesh * _mesh);
    ThicknessEvaluation(const GetPot & dfile, const Mesh * _mesh);
    ~ThicknessEvaluation();


  
  //private:


};









#endif
