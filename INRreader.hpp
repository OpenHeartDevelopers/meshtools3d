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
#include<map>
#include<limits>
#include<cfloat>
#include<cmath>

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
  IndexCoord()
  :ix(0),
  iy(0),
  iz(0),
  iv(0)
  {}
};

struct LinearRegression2DData
{
  //ax+by+cz+d=0
  std::vector<double> pnorm;
  std::vector<double> xc;
  std::vector<double> zminmax;
  std::vector<double> rhominmax;
  bool isAhole;
  LinearRegression2DData()
  :pnorm(3,0.0),
  xc(3,0.0),
  zminmax(2,0.0),
  rhominmax(2,0.0),
  isAhole(false)
  {}

};

struct eigenData
{
  std::vector<double> eigenValues;
  std::vector<double> eigenVectors;
  eigenData()
  :eigenValues(3,0.0),
  eigenVectors(9,0.0)
  {}
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
    
    ~INRInfo()
    {
      freeMemory();
    }
    void freeMemory()
    {
      SHAPE.clear();
      RESOLUTION.clear();
      TOFFSET.clear();
      ROTATION.clear();
      CENTER.clear();
      VDIM=0;
      TYPE=VX_UNKNOWN;
      PIXSIZE=0;
      littleEndian=false;
    }    

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

class BoundingBox
{
  public:
    BoundingBox()
    :_bbox(3)
    {
      _bbox.resize(3);
      for(unsigned char j=0; j<3; j++)
      {
        _bbox[j].clear();
        _bbox[j].resize(2);
        (_bbox[j])[0]=1.e32;
        (_bbox[j])[1]=0.0;
      }
    }
    ~BoundingBox()
    {
      emptyMemory();
    }
    void emptyMemory()
    {
      _bbox.clear();
    }
    inline std::vector<std::vector<double> > & bbox() {return(_bbox);};
    inline const std::vector<std::vector<double> > & bbox() const {return(_bbox);};
    inline const std::vector<double> & x() const {return(_bbox[0]);};
    inline const std::vector<double> & y() const {return(_bbox[1]);};
    inline const std::vector<double> & z() const {return(_bbox[2]);};
  private:
    std::vector<std::vector<double> >_bbox;
};


class INRreader
{
  public:
    typedef std::map<double,BoundingBox> BboxMapType;
    typedef BboxMapType::iterator BboxMapIterType;
    typedef BboxMapType::const_iterator BboxMapCIterType;
    
    INRreader(std::string filename);
    INRreader();
    ~INRreader();
    void readSegmentation(const std::string & filename);
    void freeMemory();
    void createBoundingBoxes();
    bool isPointInsideSegmentation(const double & x, const double & y, const double & z) const;
    IndexCoord voxelCoordInterp(const double & x, const double & y, const double & z) const; 
    IndexCoord voxelCoordInterpNonZero(const double & x, const double & y, const double & z);
    std::vector<double> interpolatedVoxelValue (const double & x, const double & y, const double & z);
    std::vector<double> interpolatedNonZeroVoxelValue (const double & x, const double & y, const double & z);
    inline std::vector<double> interpolatedVoxelValue(std::vector<double>& Pt) {return(interpolatedVoxelValue(Pt[0],Pt[1],Pt[2]) ); };
    inline std::vector<double> interpolatedNonZeroVoxelValue(std::vector<double>& Pt) {return(interpolatedNonZeroVoxelValue(Pt[0],Pt[1],Pt[2]) ); };
    double pickVoxelValue(const size_t & ix,const size_t & iy,const size_t & iz,const size_t iv=0) const;
    inline const INRInfo & info() const {return (_info);};
    inline const double & vx() const {return (_info.RESOLUTION[0]); };
    inline const double & vy() const {return (_info.RESOLUTION[1]); };
    inline const double & vz() const {return (_info.RESOLUTION[2]); };
    inline const size_t & xdim() const {return (_info.SHAPE[0]); };
    inline const size_t & ydim() const {return (_info.SHAPE[1]); };
    inline const size_t & zdim() const {return (_info.SHAPE[2]); };
    inline const BboxMapType & bboxlabels() const {return(_bboxlabels); };
    void printHeader();
    std::vector<double> extractPointCloud();
    std::vector<double> extractVoxelValues();
    
  private:
    bool readHeader(std::ifstream & ImageFile);
    bool readValues(std::ifstream & ImageFile);
    void evalLabeledRegionsBounds();
    double pickValue(const size_t & _index) const;
    void setValue(const size_t & _index, const double & _value);
    size_t index(const size_t & ix,const size_t & iy,const size_t & iz,const size_t & iv) const;
    IndexCoord reverseIndex(const size_t & index ) const;
    std::vector<double>evalBarycenter(const size_t & index);
    double EuclideanDist(const double * P1, const double * P2);
    LinearRegression2DData evalRegressionPlane(std::set<size_t> & pointlist);
    bool isLittleEndian();
    void SwapBytes(void *pv, size_t n);
    eigenData eigen_decomposition3X3(std::vector<double> & Amatr);
    std::vector<double> InvertA3X3(const std::vector<double> & Mat0) const;
    unsigned char RM3X3Ind(const unsigned char & irow, const unsigned char & jcol) const;
    inline double hypot2(const double & x, const double & y) const {return sqrt(x*x+y*y);};
    bool _isAllocated;
    INRInfo _info;
    size_t nb_Of_Pixels;
    size_t px_per_Plane;
    size_t byteLen;
    char * data;
    std::set<size_t> nzeroEntryIndexes;
    BboxMapType _bboxlabels;

};
#endif



