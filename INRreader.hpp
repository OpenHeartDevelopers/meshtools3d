//  -*- c++ -*-
//  INRreader version 1.0                                     Oct/21/2016
//
//
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//
//  This library is used to extract data from .inr segmentations, the 
//  INRIA segmentation format
//
//  Developer: Cesare Corrado cesare.corrado@kcl.ac.uk
//==========================================================================

#ifndef _INRREADER_H_
#define _INRREADER_H_

#include<stdio.h>
#include <fstream>
#include <sstream>
#include<string>
#include<vector>
#include<string>
#include<iostream>
#include<iomanip>
#include<set>


typedef enum{
    VX_FLOAT, // floating point
    VX_FIXED, // integers
    VX_UFIXED, //unsigned integers
    VX_UNKNOWN //unknown type
}VOX_KIND;

struct IndexCoord
{
  size_t ix;
  size_t iy;
  size_t iz;
  size_t iv;
};

class INRInfo
{

  public:
    std::vector<size_t> SHAPE;       // dimensions (XDIM,YDIM,ZDIM)
    std::vector<double> RESOLUTION;  // voxel size (VX,VY,VZ)
    std::vector<double> TOFFSET;     // Translation offset (TX,TY,TZ))
    std::vector<double> ROTATION;    // rotation (RX,RY,RZ)
    std::vector<size_t> CENTER;      // Image Center (no idea what it is!!)
    size_t VDIM;                     // vectorial dimension (VDIM)
    VOX_KIND TYPE;
    size_t PIXSIZE;
    bool littleEndian;               // Endianness
    INRInfo()
    :SHAPE(3,0),
    RESOLUTION(3,0.0),
    TOFFSET(3,0.0),
    ROTATION(3,0.0),
    CENTER(3,0),
    VDIM(0),
    TYPE(VX_UNKNOWN),
    PIXSIZE(0),
    littleEndian(false)    
    {}
    inline size_t nbOfData(){return(SHAPE[0]*SHAPE[1]*SHAPE[2]*VDIM);}
    inline size_t nbOfPixels(){return(SHAPE[0]*SHAPE[1]*SHAPE[2]);}
    inline size_t PxPerPlane(){return(SHAPE[0]*SHAPE[1]);}
    void printInfo()
    {
      std::cout<<"SHAPE:";
      for(short int j=0; j<3; j++)
      {
        std::cout<<" "<<SHAPE[j];
      }
      std::cout<<std::endl;
      
      std::cout<<"\nRESOLUTION:";
      for(short int j=0; j<3; j++)
      {
        std::cout<<" "<<std::setprecision(6)<<RESOLUTION[j];
      }
      std::cout<<std::endl;

      std::cout<<"\nOFFSET:";
      for(short int j=0; j<3; j++)
      {
        std::cout<<" "<<std::setprecision(6)<<TOFFSET[j];
      }
      std::cout<<std::endl;

      std::cout<<"\nROTATION:";
      for(short int j=0; j<3; j++)
      {
        std::cout<<" "<<std::setprecision(6)<<ROTATION[j];
      }
      std::cout<<std::endl;
      std::cout<<"\nVDIM: "<<VDIM<<" PIXSIZE: "<<PIXSIZE<<" TYPE: ";
      switch(TYPE)
      {
        case VX_FLOAT:
          std::cout<<"float";
          break;
        case VX_FIXED:
          std::cout<<"signed fixed";
          break;
        case VX_UFIXED:
          std::cout<<"unsigned fixed";
          break;
        case VX_UNKNOWN:
          std::cout<<"unknown";
          break;
      }
      
      if(littleEndian)
      {
          std::cout<<" little endian";
      }
      else
      {
          std::cout<<" big endian";
      }
      std::cout<<std::endl;
    }

};



class INRreader
{
  public:
    INRreader(std::string filename);
    INRreader();
    ~INRreader();
    void readSegmentation(const std::string & filename);
    
    bool isPointInsideSegmentation(const double & x, const double & y, const double & z) const;
    IndexCoord voxelCoordInterp(const double & x, const double & y, const double & z) const; 
    IndexCoord voxelCoordInterpNonZero(const double & x, const double & y, const double & z);
    std::vector<double> interpolatedVoxelValue (const double & x, const double & y, const double & z);
    std::vector<double> interpolatedNonZeroVoxelValue (const double & x, const double & y, const double & z);
    double pickVoxelValue(const size_t & ix,const size_t & iy,const size_t & iz,const size_t iv=0) const;
    
    inline const INRInfo & info() const {return (_info);};
    inline const double & vx() const {return (_info.RESOLUTION[0]); };
    inline const double & vy() const {return (_info.RESOLUTION[1]); };
    inline const double & vz() const {return (_info.RESOLUTION[2]); };
    inline const size_t & xdim() const {return (_info.SHAPE[0]); };
    inline const size_t & ydim() const {return (_info.SHAPE[1]); };
    inline const size_t & zdim() const {return (_info.SHAPE[2]); };
    void printHeader();
    std::vector<double> extractPointCloud();
    std::vector<double> extractVoxelValues();
    
  private:
    bool readHeader(std::ifstream & ImageFile);
    bool readValues(std::ifstream & ImageFile);
    size_t index(const size_t & ix,const size_t & iy,const size_t & iz,const size_t & iv) const;
    IndexCoord reverseIndex(const size_t & index ) const;
    std::vector<double>evalBarycenter(const size_t & index);
    double EuclideanDist(const double * P1, const double * P2);
    bool isLittleEndian();
    void SwapBytes(void *pv, size_t n);
    INRInfo _info;
    size_t nb_Of_Pixels;
    size_t px_per_Plane;
    size_t byteLen;
    char * data;
    std::set<size_t> nzeroEntryIndexes;

};
#endif



