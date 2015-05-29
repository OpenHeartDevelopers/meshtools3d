#include "LaplaceSolver.hpp"
#include<iostream>
#include "Chrono.hpp"
#include "mgmres.hpp"
#include<cmath>

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


void LaplaceSolver::setBCValue(std::set<long int> region, double value)
{
  std::set<long int>::iterator nodeiter;
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


short int LaplaceSolver::RMIndex(short int irow, short int jcol, short int rank) const
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
    if(verbose)
    {
      std::cout<<"pattern evaluation...";
    }
    
    eval_pattern();
    if(verbose)
    {
      std::cout<<"done"<<std::endl;
    }
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
    std::map<size_t,short int>::iterator it;
    for(short int iv=0; iv<4; iv++)
    {
      TetVertex[iv]=_ptrmesh->Tet(iTet).vertex[iv];
      vord.insert(std::pair<size_t, short int >(TetVertex[iv],iv) );
    }
    //Node reordering
    
    short int countv=0;
    for(it=vord.begin(); it!=vord.end(); ++it)
    {
      reordering[countv]=it->second;
      countv=countv+1;
    }
    countv=0;
    vord.clear();
    local_stiffness=localStiff(iTet);
    std::vector<double> local_RHS(4,0);
    // local matrix and RHS (imposing bcs)
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



void LaplaceSolver::writeSolution(std::string filename)
{
  std::string fname=filename+"_potential.dat";
  std::ofstream fsol(fname.c_str());
  if(!fsol)
  {
    std::cerr<<"ERROR: FILE "<<fname<<" NOT OPENED"<<std::endl;
    exit(1);
  }
  size_t nPt=_ptrmesh->nPt();
  for(size_t iPt=0; iPt<nPt; iPt++)
  {
    fsol<<_sol[iPt]<<std::endl;
  }
  fsol.close();
}

void LaplaceSolver::writeVTKSolution(std::string filename, bool binary, double rescaling)
{
  // These types are defined for the binary output
  // I suppose that paraview look for float floating point
  // and int integers; with this type definition is easier to
  // change the output type
   
  typedef int vtkIntType;
  typedef float vtkFloatType;
  
  // For an unknown idiots reason (or perhaps because network uses always big endian), 
  // paraview needs output in big endian. So, first 
  // I determine the endianness of the current machine
  

  bool littleEndianMachine=isLittleEndian();
  
  std::string outputFileName=filename+"_potential.vtk";

  short int precision=12;

  size_t nPt=_ptrmesh->nPt();
  size_t nElem=0;
  bool is_3D=false;
  vtkIntType * vertex=NULL;
  short int nbV=0;
  vtkIntType typeCell=0;
  vtkIntType nbVoutput=0;
  if(_ptrmesh->nTet()>0)
  {
    nElem=_ptrmesh->nTet();
    is_3D=true;
    nbV=4;
    typeCell=10;
  }
  else
  {
    nElem=_ptrmesh->nTri();
    is_3D=false;
    nbV=3;
    typeCell=5;
  }
  nbVoutput=static_cast<vtkIntType>(nbV);

  if(binary && littleEndianMachine )
  {
    SwapBytes(&nbVoutput, sizeof(nbVoutput));
    SwapBytes(&typeCell, sizeof(typeCell));
  }
  vertex= new vtkIntType[nbV];
  std::ofstream VTKFile;
  
  if(binary)
  {
    VTKFile.open(outputFileName.c_str(),std::ios::out | std::ios::binary);
  }
  else
  {
    VTKFile.open(outputFileName.c_str(),std::ios::out );
  }
  
  VTKFile<<"# vtk DataFile Version 4.2"<<std::endl;
  VTKFile<<"Visualization of specified geometry"<<std::endl;
  if(binary)
  {
    VTKFile<<"BINARY"<<std::endl;
  }
  else
  {
    VTKFile<<"ASCII"<<std::endl;
  }
  VTKFile<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  
  //Points
  VTKFile<<"POINTS "<<nPt<<" float"<<std::endl;
  for(size_t iPt=0; iPt<nPt; iPt++)
  {
    const Point  & Pt=_ptrmesh->Pt(iPt);
    for(short int jc=0; jc<3; jc++)
    {
        vtkFloatType pcoord=rescaling*static_cast<vtkFloatType>(Pt.coord[jc]);
        if(binary)
        {
          if(littleEndianMachine)
          {
            SwapBytes(&pcoord, sizeof(pcoord));
          }
          VTKFile.write((char*) &pcoord,sizeof(vtkFloatType));
        }
        else
        {
          VTKFile<<std::setprecision(precision)<<pcoord;
          if(jc==2)
          {
            VTKFile<<std::endl;
          }
          else
          {
            VTKFile<<" ";
          }
        }
    }
  }//end loop on points
    
  if(binary)
  {
    VTKFile<<std::endl;
  }
  
  // Elements 
  VTKFile<<"CELLS "<<nElem<<" "<<(nbV+1)*nElem<<std::endl;
  for(size_t iElem=0; iElem<nElem; iElem++)
  {
    for(vtkIntType iV=0; iV<nbV; iV++ )
    {
      if(is_3D)
      {
        vertex[iV]=static_cast<vtkIntType>((_ptrmesh->Tet(iElem)).vertex[iV]);
      }
      else
      {
        vertex[iV]=static_cast<vtkIntType>((_ptrmesh->Tri(iElem)).vertex[iV]);
      }
    }
    if(binary)
    {
      VTKFile.write((char*) &nbVoutput,sizeof(vtkIntType)); 
      for(short int iV=0; iV<nbV; iV++ )
      {
        vtkIntType VertexCP=vertex[iV];
        if(littleEndianMachine)
        {
          SwapBytes(&VertexCP, sizeof(VertexCP));
        }
        VTKFile.write((char*) &VertexCP,sizeof(vtkIntType)); 
      }
    }
    else
    {
      VTKFile<<nbV;
      for(vtkIntType iV=0; iV<nbV; iV++ )
      {
        VTKFile<<" "<<vertex[iV];
      }
      VTKFile<<std::endl;
    }
  }//end loop on elements

  delete [] vertex;
  vertex=NULL;


  if(binary)
  {
    VTKFile<<std::endl;
  }

  //Cell types
  VTKFile<<"CELL_TYPES "<<nElem<<std::endl;
  for(size_t iElem=0; iElem<nElem; iElem++)
  {
    if(binary)
    {
      VTKFile.write((char*) &typeCell,sizeof(vtkIntType)); 
    }
    else
    {
      VTKFile<<typeCell;
      if(((1+iElem)%6) && (iElem<(nElem-1)))
      {
        VTKFile<<" ";
      }
      else
      {
        VTKFile<<std::endl;
      }

    }
  }
  if(binary)
  {
    VTKFile<<std::endl;
  }
  
  //now i detrmine point label
  VTKFile<<"POINT_DATA "<<nPt<<std::endl;
  VTKFile<<"SCALARS Potential float"<<std::endl;
  VTKFile<<"LOOKUP_TABLE default"<<std::endl;
  for(size_t iPt=0; iPt<nPt; iPt++)
  {
    vtkFloatType valueAtPt=static_cast<vtkFloatType>(_sol[iPt]);
    if(binary)
    {
      if(littleEndianMachine)
      {
        SwapBytes(&valueAtPt, sizeof(valueAtPt));
      }
      VTKFile.write((char*) &valueAtPt,sizeof(vtkFloatType));
    }
    else
    {
      VTKFile<<std::setprecision(precision)<<valueAtPt;
      if(((1+iPt)%6) && (iPt<nPt-1))
      {
        VTKFile<<" ";
      }
      else
      {
        VTKFile<<std::endl;
      }
    }
  }

  if(binary)
  {
    VTKFile<<std::endl;
  }
  /*if(is_3D)
  {
    VTKFile<<"CELL_DATA "<<nElem<<std::endl;
    VTKFile<<"VECTORS Gradient float"<<std::endl;
    VTKFile<<"LOOKUP_TABLE default"<<std::endl;
    for(size_t iElem=0; iElem<nElem; iElem++)
    {
      std::vector<double> gradient= ElementTetraGradient(iElem);
      vtkFloatType vx=static_cast<vtkFloatType>(gradient[0]);
      vtkFloatType vy=static_cast<vtkFloatType>(gradient[1]);
      vtkFloatType vz=static_cast<vtkFloatType>(gradient[2]);
      gradient.clear();
      if(binary)
      {
        if(littleEndianMachine)
        {
          SwapBytes(&vx, sizeof(vx));
          SwapBytes(&vy, sizeof(vy));
          SwapBytes(&vz, sizeof(vz));
        }
        VTKFile.write((char*) &vx,sizeof(vtkFloatType));
        VTKFile.write((char*) &vy,sizeof(vtkFloatType));
        VTKFile.write((char*) &vz,sizeof(vtkFloatType));
      
      }
      else
      {
        VTKFile<<std::setprecision(precision)<<vx<<" "
               <<std::setprecision(precision)<<vy<<" "
               <<std::setprecision(precision)<<vz;
        
        if(((1+iElem)%2) && (iElem<nElem-1))
        {
          VTKFile<<" ";
        }
        else
        {
          VTKFile<<std::endl;
        }
      }
    }
    if(binary)
    {
      VTKFile<<std::endl;
    }
  }*/
  
  VTKFile.close();
}


std::vector<double> LaplaceSolver::ElementTetraGradient(size_t iTet,bool normalize) const
{
  std::vector<double> gradient(3,0);
  double  gradient0[3];
  gradient0[0]=0.0;
  gradient0[1]=0.0;
  gradient0[2]=0.0;
  const Tetrahedron & Tet= _ptrmesh->Tet(iTet);
  std::vector<double> invJt=_ptrmesh->TetInvJacobianTransponse( iTet);
  //first: eval grad0
  for(short int iv=0; iv<4; iv++)
  {
    double sol_at_vertex=_sol[Tet.vertex[iv]];
    gradient0[0]=gradient0[0]+ sol_at_vertex*(dphi0[iv])[0];
    gradient0[1]=gradient0[1]+ sol_at_vertex*(dphi0[iv])[1];
    gradient0[2]=gradient0[2]+ sol_at_vertex*(dphi0[iv])[2];
  }
  //second: eval gradient by pre-multiplying bt J^-t
  double dx=0, dy=0, dz=0;
  for(short int jc=0; jc<3; jc++)
  {
    gradient[0]=gradient[0]+invJt[RMIndex(0,jc,3)]*gradient0[jc];
    gradient[1]=gradient[1]+invJt[RMIndex(1,jc,3)]*gradient0[jc];
    gradient[2]=gradient[2]+invJt[RMIndex(2,jc,3)]*gradient0[jc];
  }
  
  if(normalize)
  {
    double norm=std::sqrt(gradient[0]*gradient[0]+gradient[1]*gradient[1]+gradient[2]*gradient[2]);
    for(short int jc=0; jc<3; jc++)
    {
      gradient[jc]=gradient[jc]/norm;
    }
  
  }
  return(gradient);
}


//used to evaluate endianess
bool LaplaceSolver::isLittleEndian()
{
  bool littleEndianMachine=true;
  int x = 1;
	if(*(char *)&x == 1)
	{
	  littleEndianMachine=true;
	}
	else
	{
	  littleEndianMachine=false;
	}
  return(littleEndianMachine);
}


//used to convert to big endian
void LaplaceSolver::SwapBytes(void *pv, size_t n)
{
    char *p = static_cast<char*>(pv);
    size_t lo, hi;
    for(lo=0, hi=n-1; hi>lo; lo++, hi--)
    {
        char tmp=p[lo];
        p[lo] = p[hi];
        p[hi] = tmp;
    }
}



void LaplaceSolver::writeElementGradient(std::string filename)
{
  short int precision=12;
  short int width=11;
  std::string fname=filename+".grad";
  std::ofstream fgrad(fname.c_str());
  if(!fgrad)
  {
    std::cerr<<"ERROR: FILE "<<fname<<" NOT OPENED"<<std::endl;
    exit(1);
  }
  size_t nTet=_ptrmesh->nTet();
  for(size_t iTet=0; iTet<nTet; iTet++)
  {
    std::vector<double> gradient= ElementTetraGradient(iTet);
    fgrad<<std::setw(width)<<std::setprecision(precision)<<std::left << std::setfill('0')<<gradient[0]<<" "
         <<std::setw(width)<<std::setprecision(precision)<<std::left << std::setfill('0')<<gradient[1]<<" "
         <<std::setw(width)<<std::setprecision(precision)<<std::left << std::setfill('0')<<gradient[2]<<std::endl;
  }
  fgrad.close();




}




