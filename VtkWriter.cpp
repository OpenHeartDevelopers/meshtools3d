#include<stdio.h>
#include<cstdlib>
#include<string>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>
#include<cmath>

#include "VtkWriter.hpp"


VtkWriter::VtkWriter( bool binary)
:M_is_binary(binary),
precision(12),
M_littleEndianMachine(true),
_meshSetted(false),
_headvarWritten(false),
_mesh(NULL),
M_nbLocalDof(0),
typeCell(0),
nbVoutput(0)
{
   M_littleEndianMachine=isLittleEndian();
   M_dir="./";
   M_prefix="atrium";
}

VtkWriter::VtkWriter(Mesh * theMesh, bool binary)
:M_is_binary(binary),
precision(12),
M_littleEndianMachine(true),
_meshSetted(false),
_headvarWritten(false),
_mesh(NULL),
M_nbLocalDof(0),
typeCell(0),
nbVoutput(0)
{
  M_littleEndianMachine=isLittleEndian();
   M_dir="./";
   M_prefix="atrium";
  setMesh(theMesh);
}


VtkWriter::~VtkWriter()
{
  _mesh=NULL;
  _meshSetted=false;
  CloseFile();
}


void VtkWriter::openFileForOutput()
{
  
  if(_meshSetted)
  {
    std::string FileName= M_dir + "/" + M_prefix+ ".vtk";
    if(outputVTKFile.is_open())
    {
      std::cerr<<"Warning: a file is already opened; closing it"<<std::endl;
      outputVTKFile.close();
    }
    std::string command="mkdir -p " + M_dir;
    int ierr=system(command.c_str());
    if(ierr)
    {
        std::cerr<<"ERROR: directory "<< M_dir<<" not created"<<std::endl;
        exit(1);
    }
    std::cout<<"Opening File "<<FileName<<std::endl;
    if (M_is_binary)
    {
      outputVTKFile.open(FileName.c_str(),std::ios::out | std::ios::binary);
    }
    else
    {
      outputVTKFile.open(FileName.c_str(),std::ios::out );
    }
    std::cout<<"Writing header"<<std::endl;
    M_write_header(outputVTKFile, M_prefix );
    std::cout<<"Writing geometry"<<std::endl;
    M_write_geo(outputVTKFile);
  }
  else
  {
    std::cerr<<"ERROR: No mesh was assigned"<<std::endl;
    exit(1);
  }
}


void VtkWriter::CloseFile()
{
  if(outputVTKFile.is_open())
  {
    _headvarWritten = false;
    outputVTKFile.close();
  }
}


void VtkWriter::setMesh(Mesh * theMesh)
{
  _meshSetted=true;
  _mesh=theMesh;
  if(_mesh->nTet()>0)
  {
    typeCell=10;
    M_nbLocalDof=4;
  }
  else
  {
    if(_mesh->nTri()>0)
    {
        typeCell=5;
        M_nbLocalDof=3;
    }
  }
  nbVoutput=static_cast<vtkIntType>(M_nbLocalDof);
  if(M_is_binary && M_littleEndianMachine)
  {
    SwapBytes(&nbVoutput, sizeof(nbVoutput));
    SwapBytes(&typeCell, sizeof(typeCell));
  }
}



void VtkWriter::setOutput(std::string theDir,std::string thePrefix)
{
  setOutputDir(theDir);
  setPrefixName(thePrefix);
}

void VtkWriter::writeVariable(const unkVecType &  dvar, std::string nameVar, Type varType)
{
  if(!_headvarWritten)
  {
    M_write_var_head( outputVTKFile );
    _headvarWritten = true;  
  }
  if(varType==Scalar)
  {
    M_write_scalar( dvar, outputVTKFile, nameVar);
  }
  if(varType==Vector)
  {
    M_write_vector( dvar, outputVTKFile, nameVar);
  }
}


//Write the header of the file
void VtkWriter::M_write_header( std::ofstream & VTKFile,std::string infoStr )
{
  VTKFile<<"# vtk DataFile Version 4.2"<<std::endl;
  VTKFile<<infoStr<<std::endl;
  if(M_is_binary)
  {
    VTKFile<<"BINARY"<<std::endl;
  }
  else
  {
    VTKFile<<"ASCII"<<std::endl;
  }
}


//Write the Geometry
void VtkWriter::M_write_geo( std::ofstream & VTKFile )
{
  size_t np   = _mesh->nPt();
  size_t nTet = _mesh->nTet();
  size_t nTri = _mesh->nTri();
  size_t nElem=0;
  VTKFile<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  //Points
  VTKFile<<"POINTS "<<np<<" float"<<std::endl;
  for(size_t iPt=0; iPt<np; iPt++)
  {
    for(int jc=0; jc<3; jc++)
    {
      vtkFloatType pcoord=static_cast<vtkFloatType>(_mesh->Pt(iPt).coord[jc]);
      if(M_is_binary)
      {
        if(M_littleEndianMachine)
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
  
  if(M_is_binary)
  {
    VTKFile<<std::endl;
  }
  const Element * ElemVec;
  if(nTet>0)
  {
    nElem=nTet;
    ElemVec=(_mesh->Tet().data() );
  }
  else
  {
    nElem=nTri;
    ElemVec=(_mesh->Tri().data());
  }
  // Elements
  VTKFile<<"CELLS "<<nElem<<" "<<(M_nbLocalDof+1)*nElem<<std::endl;
  for(size_t iElem=0; iElem<nElem; iElem++)
  {
    if(M_is_binary)
    {
      VTKFile.write((char*) &nbVoutput,sizeof(vtkIntType));
    }
    else
    {
      VTKFile<<nbVoutput;
    }
    for(short int iV=0; iV<ElemVec[iElem].nbV(); iV++ )
    {
      vtkIntType Vertex=static_cast<vtkIntType>(ElemVec[iElem].vertex[iV]);
      if(M_is_binary)
      {
        if(M_littleEndianMachine)
        {
          SwapBytes(&Vertex, sizeof(Vertex));
        }
        VTKFile.write((char*) &Vertex,sizeof(vtkIntType));
      }
      else
      {
        VTKFile<<" "<<Vertex;
        if(iV==(ElemVec[iElem].nbV()-1))
        {
          VTKFile<<std::endl;
        }
      }
    }
  }//end loop on elements
  if(M_is_binary)
  {
    VTKFile<<std::endl;
  }
  //Cell types
  VTKFile<<"CELL_TYPES "<<nElem<<std::endl;
  for(size_t iElem=0; iElem<nElem; iElem++)
  {
    if(M_is_binary)
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
}

void VtkWriter::M_write_var_head( std::ofstream & VTKFile )
{
  VTKFile<<"POINT_DATA "<<_mesh->nPt()<<std::endl;
}


//used to evaluate endianess
bool VtkWriter::isLittleEndian()
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
void VtkWriter::SwapBytes(void *pv, size_t n)
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



void VtkWriter::M_write_scalar(const unkVecType &  dvar, std::ofstream & VTKFile, std::string nameVar)
{
  VTKFile<<"SCALARS "<<nameVar<<" FLOAT"<<std::endl;
  VTKFile<<"LOOKUP_TABLE default"<<std::endl;
  size_t dim = _mesh->nPt();
  for (size_t ip=0;ip<dim;++ip)
  {
    vtkFloatType valueAtPt=static_cast<vtkFloatType>( dvar[ip]);
    if(M_is_binary)
    {
      if(M_littleEndianMachine)
      {
        SwapBytes(&valueAtPt, sizeof(valueAtPt));
      }
      VTKFile.write((char*) &valueAtPt,sizeof(vtkFloatType));
    }
    else
    {
      VTKFile<<std::setprecision(precision)<<valueAtPt;
      if(((1+ip)%6) && (ip<dim-1))
      {
        VTKFile<<" ";
      }
      else
      {
        VTKFile<<std::endl;
      }
    }
  }
  VTKFile<<std::endl;
}


void VtkWriter::M_write_vector(const unkVecType &  dvar, std::ofstream & VTKFile, std::string nameVar)
{
  // consider vector variables as var_x,var_y,var_z (3 chuncks)
  VTKFile<<"VECTORS "<<nameVar<<" FLOAT"<<std::endl;
  VTKFile<<"LOOKUP_TABLE default"<<std::endl;
  size_t dim = _mesh->nPt();
  for (size_t ip=0;ip<dim;++ip)
  {
   for (short int jdim=0; jdim< 3;++jdim)
   {
     vtkFloatType valueAtPt=static_cast<vtkFloatType>( dvar[jdim*dim+ip]);
     if(M_is_binary)
     {
       if(M_littleEndianMachine)
       {
         SwapBytes(&valueAtPt, sizeof(valueAtPt));
       }
       VTKFile.write((char*) &valueAtPt,sizeof(vtkFloatType));
     }
     else
     {
       VTKFile<<std::setprecision(precision)<<valueAtPt;
       if(((1+ip)%6) && (ip<dim-1))
       {
         VTKFile<<" ";
       }
       else
       {
         VTKFile<<std::endl;
       }
     }
   }
  }
  VTKFile<<std::endl;
}


/*


void VtkWriter::M_write_output(std::string FileName, float ctime)
{
  std::ofstream VTKFile;
  std::ostringstream timestr;
  timestr<<ctime;
  if (M_is_binary)
  {
    VTKFile.open(FileName.c_str(),std::ios::out | std::ios::binary);
  }
  else
  {
    VTKFile.open(FileName.c_str(),std::ios::out );
  }
  M_write_header(VTKFile,timestr.str());
  M_write_geo(VTKFile);
  M_write_variables(VTKFile);
  VTKFile.close();
}

void VtkWriter::M_write_variables(std::ofstream & VTKFile)
{
  VTKFile<<"POINT_DATA "<<_mesh->num_Points()<<std::endl;
}










*/

