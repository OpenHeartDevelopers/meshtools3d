#ifndef _VTKWRITER_H_
#define _VTKWRITER_H_

#include "Mesh.hpp"
#include<stdio.h>
#include <fstream>
#include <sstream>
#include<string>
#include<vector>


class VtkWriter {
    typedef int vtkIntType;
    typedef float vtkFloatType;
    typedef std::vector<double>  unkVecType;
    
    public:
      enum Type{Scalar,Vector};
      VtkWriter( bool binary=false);
      VtkWriter(Mesh * theMesh, bool binary=false);
      ~VtkWriter();
      void setMesh(Mesh * theMesh);
      inline void setOutputDir(std::string theDir){M_dir=theDir;};
      inline void setPrefixName(std::string thePrefix){M_prefix=thePrefix;};
      void setOutput(std::string theDir,std::string thePrefix="atrium");
      void openFileForOutput();
      void writeVariable(const unkVecType &  dvar, std::string nameVar, Type varType );
      void CloseFile();
/*      
      
      
*/
      
    /*  
      
      
      void postProcess(const double & time);
  */

    
   private:
    bool isLittleEndian();
    void SwapBytes(void *pv, size_t n);
    
    // aux functions for writing parts of vtk file
    void M_write_header( std::ofstream & VTKFile,std::string infoStr );
    void M_write_geo( std::ofstream & VTKFile );
    void M_write_var_head( std::ofstream & VTKFile );
    void M_write_scalar(const unkVecType &  dvar, std::ofstream & VTKFile, std::string nameVar);
    void M_write_vector(const unkVecType &  dvar, std::ofstream & VTKFile, std::string nameVar);
    bool M_is_binary;
    short int precision;
    bool M_littleEndianMachine;
    bool _meshSetted;
    bool _headvarWritten;
    const Mesh *  _mesh;
    
    short int M_nbLocalDof;
    vtkIntType typeCell;
    vtkIntType nbVoutput;
    std::string M_dir;
    std::string M_prefix;
    std::ofstream outputVTKFile;

};
  
#endif
