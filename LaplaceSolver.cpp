#include "LaplaceSolver.hpp"
#include<iostream>
#include "Chrono.hpp"
#include "mgmres.hpp"

LaplaceSolver::LaplaceSolver()
:_consistentState(false),
_ptrmesh(NULL),
abs_toll(1e-6),
rel_toll(1e-6),
itr_max(100),
dimKrilovSp(10),
verbose(0)
{}

LaplaceSolver::LaplaceSolver(const Mesh *  _mesh)
:_consistentState(false),
_ptrmesh(NULL),
abs_toll(1e-6),
rel_toll(1e-6),
itr_max(100),
dimKrilovSp(10),
verbose(0)
{
  this->setMesh(_mesh);
}

void LaplaceSolver::setMesh(const Mesh *  _mesh)
{
  if(_mesh)
  {
    _ptrmesh=_mesh;
    _consistentState=true;
    this->evaldphi0();
    _sol.resize(_ptrmesh->nPt(),0);
    _RHS.resize(_ptrmesh->nPt(),0);
  }
}

LaplaceSolver::LaplaceSolver(const GetPot & dfile, const Mesh * _mesh)
:_consistentState(false),
_ptrmesh(NULL),
abs_toll(1e-6),
rel_toll(1e-6),
itr_max(100),
dimKrilovSp(10),
verbose(0)
{
  this->setMesh(_mesh);
  abs_toll    = dfile("laplacesolver/abs_toll",1.e-6);
  rel_toll    = dfile("laplacesolver/rel_toll",1.e-6);
  itr_max     = dfile("laplacesolver/itr_max",100);    
  dimKrilovSp = dfile("laplacesolver/dimKrilovSp",100);
  verbose     = dfile("laplacesolver/verbose",0);
  if(dimKrilovSp>static_cast<long int>(_ptrmesh->nPt()))
  {
    dimKrilovSp=_ptrmesh->nPt();
    std::cout<<"Warning: Krilov space dimension was bigger than matrix size; changed to "<<_ptrmesh->nPt()<<std::endl;
  }
}


void LaplaceSolver::setBCValue(std::set<size_t> region, double value)
{
  std::set<size_t>::iterator nodeiter;
  for(nodeiter=region.begin(); nodeiter!=region.end(); ++nodeiter)
  {
    DirichletBC.insert(std::pair<long int,double>(*nodeiter,value));
  }
}
void LaplaceSolver::setBCValue(BCContainerType BCS)
{
  BCContainerTypeIterator mapiter;
  for(mapiter=BCS.begin(); mapiter!=BCS.end(); ++mapiter)
  {
    DirichletBC.insert(*mapiter);
  }
}



void LaplaceSolver::eval_pattern()
{
  _Matrix._pattern.I.resize(1+_ptrmesh->nPt());
  std::vector<std::set<long int> > labelPerNode;
  labelPerNode.resize(_ptrmesh->nPt());
  // I create a vector of set
  // in each entry the set gives the connectivity
  // i.e. the label of neighb. nodes

  for(size_t iTet=0; iTet<_ptrmesh->nTet(); iTet++)  
  {
    for(short int iPt=0; iPt<4; iPt++)
    {
        for(short int jPt=iPt; jPt<4; jPt++)
        {
          size_t iNode=_ptrmesh->Tet(iTet).vertex[iPt];
          size_t jNode=_ptrmesh->Tet(iTet).vertex[jPt];
          labelPerNode[iNode].insert(jNode);
          labelPerNode[jNode].insert(iNode);
        }
    }
  }
  // I evaluate the number of non_zero elements into the pattern
  _Matrix._pattern.n_zero=0;
  for(size_t iPt=0; iPt<_ptrmesh->nPt(); iPt++)
  {
    _Matrix._pattern.n_zero=_Matrix._pattern.n_zero+labelPerNode[iPt].size();
  }
  _Matrix._pattern.J.resize(_Matrix._pattern.n_zero);
  _Matrix.K.resize(_Matrix._pattern.n_zero,0);
  // Here I  fill the pattern of the matrix
  _Matrix._pattern.I[0]=0;
  for(size_t iPt=0; iPt<_ptrmesh->nPt(); iPt++)
  {
    _Matrix._pattern.I[iPt+1]=_Matrix._pattern.I[iPt]+labelPerNode[iPt].size();
    
    unsigned counter=0;
    for(std::set<long int>::iterator it=labelPerNode[iPt].begin(); it!=labelPerNode[iPt].end(); ++it)
    {
      _Matrix._pattern.J[_Matrix._pattern.I[iPt]+counter]=*it;
      counter=counter+1;
    }
  }
}

std::vector<double> LaplaceSolver::localStiff(size_t iTet, double k)
{
  std::vector<double> lstif(16,0); //4X4 matrix
  double Volume=(_ptrmesh->VolTet(iTet));
  std::vector<double> invJt=(_ptrmesh->TetInvJacobianTransponse(iTet)); // 9X1
  std::vector<std::vector<double> > dphi;
  dphi.resize(4);
  double coef = k*Volume/6.0;
  for(short int jf=0; jf<4; jf++)
  {
    dphi[jf].resize(3,0);
    for(short int ic=0; ic<3; ic++)
    {
      for(short int jc=0; jc<3; jc++)
      {
        (dphi[jf])[ic]=(dphi[jf])[ic]+invJt[RMIndex(ic,jc,3)]  *(dphi0[jf])[jc];
      }
    }
  }
  // Now eval the matrix
  for(short int ir=0; ir<4; ir++)
  {
    for(short int jc=ir; jc<4; jc++)
    {
      double entry=0;
      for(short int icoor=0; icoor<3; icoor++)
      {
        entry=entry+((dphi[ir])[icoor])*((dphi[jc])[icoor]);
      }
      lstif[RMIndex(ir,jc,4)]=coef*entry;
      lstif[RMIndex(jc,ir,4)]=coef*entry;
    }
  }
  return(lstif);
}


LaplaceSolver::~LaplaceSolver()
{
  if(_ptrmesh)
  {
    _ptrmesh=NULL;
  }
  _Matrix.clear(true);
  dphi0.clear();
  _sol.clear();
  _RHS.clear();
}


short int LaplaceSolver::RMIndex(short int irow, short int jcol, short int rank)
{
  short int index=irow*rank+jcol;
  return(index);
}

void LaplaceSolver::evaldphi0()
{
  dphi0.resize(4);
  for(short int jf=0; jf<4; jf++)
  {
      dphi0[jf].resize(3);
  }
  for(short int jc=0; jc<3; jc++)
  {
    dphi0[0][jc]=-1;
    dphi0[jc+1][jc]=1;
  }
}




void LaplaceSolver::matrixAssembly(bool build_pattern)
{

  Chrono chrono;
  if(verbose)
  {
    std::cout<<"Assembling the matrix"<<std::endl;
  }
  chrono.start();
    if(build_pattern)
  {
    eval_pattern();
  }
  
  short int* reordering=new short int [4];
  for(size_t iTet=0; iTet<_ptrmesh->nTet(); iTet++)
  {
    std::vector<double>local_stiffness(16,0);
    reordering[0]=0;
    reordering[1]=1;
    reordering[2]=2;
    reordering[3]=3;
    std::vector<size_t> TetVertex(4,0);
    std::map<size_t,short int> vord;
    for(short int iv=0; iv<4; iv++)
    {
      TetVertex[iv]=_ptrmesh->Tet(iTet).vertex[iv];
      vord.insert(std::pair<size_t,short int>(TetVertex[iv],iv));
    }
    //Node reordering
    std::map<size_t,short int>::iterator it;
    short int countv=0;
    for(it=vord.begin(); it!=vord.end(); ++it)
    {
      reordering[countv]=it->second;
      countv++;
    }
    countv=0;
    vord.clear();
    local_stiffness=localStiff(iTet);
    std::vector<double> local_RHS(3,0);
    // local matrix and RHS
    for(short int iPt=0; iPt<4; iPt++)
    {
        long int global_vertex=TetVertex[iPt];
        std::map<long int, double>::iterator itv;
        itv=DirichletBC.find(global_vertex);
        //there is a boundary condition
        if(!(itv==DirichletBC.end()))
        {
          double aii=local_stiffness[RMIndex(iPt,iPt,4)];
          local_stiffness[RMIndex(iPt,iPt,4)]=0.0;
          //simmetric diagonalization bcs
          for(short int jPt=0; jPt<4; jPt++)
          {
            local_stiffness[RMIndex(iPt,jPt,4)]=0.0;  //rows to zero
            double aji=local_stiffness[RMIndex(jPt,iPt,4)];
            // simm. diag step
            local_RHS[jPt]=local_RHS[jPt]-aji*(itv->second); 
            local_stiffness[RMIndex(jPt,iPt,4)]=0.0;
          }
          local_stiffness[RMIndex(iPt,iPt,4)]=aii; 
          local_RHS[iPt]=aii*(itv->second);
        }
    }//end for on local points

    // matrix and RHS connectivity
    for(short int iPt=0; iPt<4; iPt++)
    {
      long int I_index=TetVertex[reordering[iPt]];
      _RHS[I_index]=_RHS[I_index]+local_RHS[reordering[iPt]];
      for(short int jPt=0; jPt<4; jPt++)
      {
        long int start=_Matrix._pattern.I[I_index];
        long int end=_Matrix._pattern.I[I_index+1]-1;
        long int J_index=TetVertex[reordering[jPt]];
        double matrix_entry=local_stiffness[RMIndex(reordering[iPt],reordering[jPt],4)];
        bool found=false;
        long int tmpIndex=0.5*(start+end);
        size_t numel=end-start;
        if(_Matrix._pattern.J[start]==J_index)
        {
          found=true;
          tmpIndex=start;
        }
        if(!found)
        {
          if(_Matrix._pattern.J[end]==J_index)
          found=true;
          tmpIndex=end;
        }
        if(!found)
        {
          if(_Matrix._pattern.J[tmpIndex]==J_index)
          found=true;
        }
        size_t counter=0;
        while(!found && counter<=numel)
        {
          if(J_index<_Matrix._pattern.J[tmpIndex])  
          {
            end=tmpIndex;
          }
          else
          {
            start=tmpIndex;
          }
          tmpIndex=0.5*(end+start);
          if(_Matrix._pattern.J[tmpIndex]==J_index)
          {
            found=true;
          }
          counter=counter+1;
        }
        if(found)
        {
          _Matrix.K[tmpIndex]=_Matrix.K[tmpIndex]+matrix_entry;
        }
        else
        {
          std::cerr<<"SOMETHING GOES WRONG"<<std::endl;
          exit(1);
        }
      }//end loop on jPt
    }//end loop on iPt
    
    

  }//end loop on tetra
  
  delete [] reordering;
  reordering=NULL;
  chrono.stop();
  std::cout<<" done in "<<chrono<<std::endl;
  chrono.reset();
}


void LaplaceSolver::solve()
{
  matrixAssembly();
  _sol=_RHS;
  Chrono chrono;
  long int nPt=_ptrmesh->nPt();
  std::cout<<"Solving the linear system (GMRES)"<<std::endl;
  chrono.start();
  pmgmres_ilu_cr( nPt, _Matrix._pattern.n_zero, _Matrix._pattern.I.data(), _Matrix._pattern.J.data(), _Matrix.K.data(), 
                    _sol.data(), _RHS.data(), itr_max, dimKrilovSp,  abs_toll, 
                    rel_toll, verbose );
  chrono.stop();
  std::cout<<" done in "<<chrono<<std::endl;
  chrono.reset();
}

