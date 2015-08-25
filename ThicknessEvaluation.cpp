#include "ThicknessEvaluation.hpp"
#include<iostream>
#include "Chrono.hpp"
#include "mgmres.hpp"
#include<cmath>


ThicknessEvaluation::ThicknessEvaluation()
:LaplaceSolver()
{}

ThicknessEvaluation::ThicknessEvaluation(const Mesh *  _mesh)
:LaplaceSolver(_mesh)
{
}


ThicknessEvaluation::ThicknessEvaluation(const GetPot & dfile, const Mesh * _mesh)
:LaplaceSolver(dfile, _mesh)
{
}

ThicknessEvaluation::~ThicknessEvaluation()
{}
