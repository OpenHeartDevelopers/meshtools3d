#include<stdio.h>
#include <fstream>
#include <sstream>
#include<string.h>
#include<vector>
#include<string>
#include<iostream>
#include <limits>
#include <algorithm>
#include <iomanip>
#include<cmath>

#include "INRreader.hpp"
#ifndef DELTAMAX
#define DELTAMAX 5
#endif

INRreader::INRreader(std::string filename)
:nb_Of_Pixels(0),
px_per_Plane(0),
byteLen(0),
data(NULL)
{
    nzeroEntryIndexes.clear();
    std::ifstream INRfile(filename.c_str(),std::ios::in |std::ios::binary);
    if(!INRfile)
    {
        std::cerr<<"ERROR: FILE "<<filename<<" NOT OPENED! "<<std::endl;
        exit(1);
    }
    
    bool headerRead=readHeader(INRfile);
    if(!headerRead)
    {
            std::cerr<<"ERROR: PROBLEMS WITH THE IMAGE HEADER"<<std::endl;
            exit(1);
    }
    
    nb_Of_Pixels=_info.nbOfPixels();
    px_per_Plane=_info.PxPerPlane();
    byteLen=(_info.PIXSIZE)/8;
    data = new char [byteLen*(_info.nbOfData())];
    bool valuesRead=readValues(INRfile);
    if(!valuesRead)
    {
            std::cerr<<"ERROR: PROBLEMS WITH THE IMAGE VALUES"<<std::endl;
            exit(1);
    }
    

    INRfile.close();
    
}

INRreader::~INRreader()
{
  delete [] data;
  data=NULL;
  byteLen=0;
  nzeroEntryIndexes.clear();
  nb_Of_Pixels=0;
  px_per_Plane=0;
}

IndexCoord INRreader::voxelCoordInterp(const double & x, const double & y, const double & z)
{
  IndexCoord Ixyz;
  Ixyz.iv=0;
  Ixyz.ix=static_cast<size_t>(x/_info.RESOLUTION[0]);
  Ixyz.iy=static_cast<size_t>(y/_info.RESOLUTION[1]);
  Ixyz.iz=static_cast<size_t>(z/_info.RESOLUTION[2]);
  if(Ixyz.ix>=_info.SHAPE[0])
  {
    Ixyz.ix=_info.SHAPE[0]-1;
  }

  if(Ixyz.iy>=_info.SHAPE[1])
  {
    Ixyz.iy=_info.SHAPE[1]-1;
  }

  if(Ixyz.iz>=_info.SHAPE[2])
  {
    Ixyz.iz=_info.SHAPE[2]-1;
  }
  return(Ixyz);
}

IndexCoord INRreader::voxelCoordInterpNonZero(const double & x, const double & y, const double & z)
{
    IndexCoord Ixyz=voxelCoordInterp(x,y,z);
    bool iszero=true;
    for(size_t iv=0; iv<_info.VDIM; iv++)
    {
      long long int voxval=static_cast<long long int>(pickVoxelValue(Ixyz.ix,Ixyz.iy,Ixyz.iz,iv));
      if(voxval>0)
      {
        iszero=false;
        break;
      }
    }
    if(iszero)
    {
      std::vector<size_t> ixrange(2,0),iyrange(2,0),izrange(2,0);
      if(Ixyz.ix<DELTAMAX)
      {
        ixrange[0]=0;
      }
      else
      {
        ixrange[0]=Ixyz.ix-DELTAMAX;
      }
      
      if(Ixyz.iy<DELTAMAX)
      {
        iyrange[0]=0;
      }
      else
      {
        iyrange[0]=Ixyz.iy-DELTAMAX;
      }


      if(Ixyz.iz<DELTAMAX)
      {
        izrange[0]=0;
      }
      else
      {
        izrange[0]=Ixyz.iz-DELTAMAX;
      }
      ixrange[1]=std::min(ixrange[0]+2*DELTAMAX,_info.SHAPE[0]);
      iyrange[1]=std::min(iyrange[0]+2*DELTAMAX,_info.SHAPE[1]);
      izrange[1]=std::min(izrange[0]+2*DELTAMAX,_info.SHAPE[2]);
      std::set<size_t> researchInterval;
      for(size_t indx=ixrange[0]; indx<ixrange[1]; indx++ )
      {
          for(size_t indy=iyrange[0]; indy<iyrange[1]; indy++ )
          {
              for(size_t indz=izrange[0]; indz<izrange[1]; indz++ )
              {
                for(size_t iv=0; iv<_info.VDIM; iv++)
                {
                    size_t indexCandidate=index(indx,indy,indz, iv);
                    if((nzeroEntryIndexes.find(indexCandidate))!=nzeroEntryIndexes.end())
                    {
                        researchInterval.insert(indexCandidate);
                    }
                }
              }
          }
      }
      
      if(researchInterval.empty())
      {
          std::cerr<<" SET INTERSECTION EMPTY!!"<<std::endl;
          exit(1);
      }
      double * P1 = new double [3];
      P1[0]=x;
      P1[1]=y;
      P1[2]=z;
      std::set<size_t>::iterator it=researchInterval.begin();
      std::vector<double> bar=evalBarycenter(*it);
      size_t ind=*it;
      double dist0=EuclideanDist(P1, bar.data());
      for(it=researchInterval.begin();it!=researchInterval.end(); ++it)
      {
          bar=evalBarycenter(*it);
          double dist1=EuclideanDist(P1, bar.data());
          if(dist1<dist0 )
          {
              dist0=dist1;
              ind=*it;
          }
      }
      delete [] P1;
      P1=NULL;
      Ixyz=reverseIndex(ind);
    }
    return(Ixyz);
}

std::vector<double> INRreader::interpolatedVoxelValue (const double & x, const double & y, const double & z)
{
  std::vector<double> voxvalue(_info.VDIM,0.0);
  IndexCoord voxcoord=voxelCoordInterp(x,y,z);
  for(size_t iv=0; iv<_info.VDIM; iv++)
  {
      voxvalue[iv]=pickVoxelValue(voxcoord.ix,voxcoord.iy,voxcoord.iz, iv);
  }
  return(voxvalue);
}

std::vector<double> INRreader::interpolatedNonZeroVoxelValue (const double & x, const double & y, const double & z)
{
  std::vector<double> voxvalue(_info.VDIM,0.0);
  IndexCoord nzerocoord=voxelCoordInterpNonZero(x,y,z);
  for(size_t iv=0; iv<_info.VDIM; iv++)
  {
      voxvalue[iv]=pickVoxelValue(nzerocoord.ix,nzerocoord.iy,nzerocoord.iz, iv);
  }
  return(voxvalue);
}

double INRreader::pickVoxelValue(const size_t & ix,const size_t & iy,const size_t & iz,size_t iv)
{
    char * value= new char[byteLen];
    size_t byte_offset=byteLen*index(ix,iy,iz, iv);
    for(size_t ibyte=0; ibyte<byteLen; ibyte++ )
    {
        value[ibyte]=data[byte_offset+ibyte];
    }
    double res=0.0;    

    switch(_info.TYPE)
    {
        case VX_FLOAT:
        {
          memcpy (&res, value, byteLen);      
          break;
        }
        case VX_FIXED:
        {
          long long int lintres=0;
          memcpy (&lintres, value, byteLen);      
          res=static_cast<double>(lintres);
          break;
        }
        case VX_UFIXED:
        {
          long long unsigned luires=0;
          memcpy (&luires, value, byteLen);
          res=static_cast<double>(luires);      
          break;
        }
        default:
        {
          std::cerr<<"unknown type"<<std::endl;
          exit(1);
          break;
        }
    }
    delete [] value;
    value = NULL;
    return(res);
}

void INRreader::printHeader()
{
    _info.printInfo();
}

std::vector<double> INRreader::extractPointCloud()
{
  std::vector<double> barycentres(3*nzeroEntryIndexes.size(),0.0);
  size_t counter=0;
  for(std::set<size_t>::iterator it=nzeroEntryIndexes.begin(); it!=nzeroEntryIndexes.end(); ++it)
  {
      std::vector<double> localBar=evalBarycenter(*it);
      barycentres[3*counter]=localBar[0];
      barycentres[3*counter+1]=localBar[1];
      barycentres[3*counter+2]=localBar[2];
      counter=counter+1;
  }
  return(barycentres);
}

std::vector<double> INRreader::extractVoxelValues()
{

  std::vector<double> voxValues(nzeroEntryIndexes.size(),0.0);
  
  size_t counter=0;
  for(std::set<size_t>::iterator it=nzeroEntryIndexes.begin(); it!=nzeroEntryIndexes.end(); ++it)
  {
      IndexCoord intcoord=reverseIndex(*it );
      double Value=pickVoxelValue(intcoord.ix,intcoord.iy,intcoord.iz,intcoord.iv);
      voxValues[counter]=Value;
      counter=counter+1;
  }
  return(voxValues);
}


// Private member functions

bool INRreader::readHeader(std::ifstream & ImageFile)
{
    bool readerisok=false;
    bool headerbeginok=false;
    bool headerendok=false;
    std::string lineHeader="";
    const std::string endOfHeader="##}";
    do{
        std::getline(ImageFile,lineHeader,'\n');
        
        if(!headerbeginok)
        {
            if(lineHeader.compare("#INRIMAGE-4#{")==0)
            {
              headerbeginok=true;
            }
        }
        if(ImageFile.eof())
        {
            std::cerr<<"ERROR: PROBLEMS WITH THE IMAGE FILE"<<std::endl;
            exit(1);
        }
        
        std::string delimiter="=";
        size_t pos=lineHeader.find(delimiter);
        if(pos != std::string::npos)
        {
            //substr(size_t pos, size_t n ) is a substring starting at pos an with length n
            std::string token = lineHeader.substr(0, pos);
            std::string value = lineHeader.substr(pos+ delimiter.length(),std::string::npos);  
            //move token to uppercase
            std::transform(token.begin(), token.end(),token.begin(), ::toupper);
            
            // image size (SHAPE)
            if (token.compare("XDIM") == 0)
            {
              std::stringstream tonumber(value);
              tonumber>>(_info.SHAPE[0]);
            }
            else if (token.compare("YDIM") == 0)
            { 
              std::stringstream tonumber(value);
              tonumber>>(_info.SHAPE[1]);
            }
            else if (token.compare("ZDIM") == 0)
            { 
              std::stringstream tonumber(value);
              tonumber>>(_info.SHAPE[2]);
            }
            
            //voxel size (RESOLUTION)
            else if (token.compare("VX") == 0)
            {
              std::stringstream tonumber(value);
              tonumber>>(_info.RESOLUTION[0]);
            }
            else if (token.compare("VY") == 0)
            { 
              std::stringstream tonumber(value);
              tonumber>>(_info.RESOLUTION[1]);
            }
            else if (token.compare("VZ") == 0)
            { 
               std::stringstream tonumber(value);
               tonumber>>(_info.RESOLUTION[2]);
            }
            
            //translation offset (TOFFSET)
            else if (token.compare("TX") == 0)
            {
              std::stringstream tonumber(value);
              tonumber>>(_info.TOFFSET[0]);
            }
            else if (token.compare("TY") == 0)
            { 
              std::stringstream tonumber(value);
              tonumber>>(_info.TOFFSET[1]);
            }
            else if (token.compare("TZ") == 0)
            { 
              std::stringstream tonumber(value);
              tonumber>>(_info.TOFFSET[2]);
            }
            
            //image rotation (ROTATION)
            else if (token.compare("RX") == 0)
            {
              std::stringstream tonumber(value);
              tonumber>>(_info.ROTATION[0]);
            }
            else if (token.compare("RY") == 0)
            { 
              std::stringstream tonumber(value);
              tonumber>>(_info.ROTATION[1]);
            }
            else if (token.compare("RZ") == 0)
            { 
              std::stringstream tonumber(value);
              tonumber>>(_info.ROTATION[2]);
            }
            
            // vectorial dimension (VDIM)
            else if (token.compare("VDIM") == 0)
            {
              std::stringstream tonumber(value);
              tonumber>>(_info.VDIM);
            
            }
            //scale; don't know what is it
            else if (token.compare("SCALE") == 0)
            {
                //do nothing
            }
            //type of the pixel (TYPE)
            else if (token.compare("TYPE") == 0)
            {
                std::transform(value.begin(), value.end(),value.begin(), ::toupper);
                size_t pos2=value.find(" ");
                std::string datatype=value.substr(0,pos2);
                if(datatype.compare("FLOAT") == 0) 
                {
                  _info.TYPE=VX_FLOAT;
                }
                else if(datatype.compare("SIGNED") == 0)
                {
                  _info.TYPE=VX_FIXED;
                } 
                else if(datatype.compare("UNSIGNED") == 0) 
                {
                  _info.TYPE=VX_UFIXED;
                }
                else
                {
                  _info.TYPE=VX_UNKNOWN;
                }
            }
            //size of the pixel (PIXSIZE)
            else if (token.compare("PIXSIZE") == 0)
            {
              size_t pos2=value.find(" ");
              std::string bitsize=value.substr(0,pos2);
              std::stringstream tonumber(bitsize);
              tonumber>>(_info.PIXSIZE);
            }
            
            //endianess (CPU)
            else if (token.compare("CPU") == 0)
            {
                std::transform(value.begin(), value.end(),value.begin(), ::toupper);
                if( (value.compare("DECM") == 0) ||  (value.compare("ALPHA") == 0) ||  (value.compare("PC") == 0))
                {
                    _info.littleEndian=true;
                }
                else
                {
                    _info.littleEndian=false;
                }
            }
        }
    }while(endOfHeader.compare(lineHeader)!=0);
    
    if(endOfHeader.compare(lineHeader)==0) 
    {
        headerendok=true;
    }
    readerisok=(headerendok && headerbeginok);
    return(readerisok);
}

bool INRreader::readValues(std::ifstream & ImageFile)
{
  bool readerisok=true;
  bool _amIlittleEndian=isLittleEndian();
  bool _swapBites=true;
  
  if((_amIlittleEndian && _info.littleEndian) ||(!(_amIlittleEndian || _info.littleEndian)) )
  {
      _swapBites=false;
  }
  size_t nb_of_elements=_info.nbOfData();
  char * value= new char[byteLen];
  std::cout<<"Reading "<<nb_of_elements<<" Values..."<<std::endl;
  for(size_t iEntry=0; iEntry<nb_of_elements; iEntry++)
  {
      if(ImageFile.eof())
      {
            readerisok=false;
            break;
      }
      
      ImageFile.read(value,byteLen);
      if(_swapBites)
      {
          SwapBytes(value, byteLen);
      }
      size_t byte_offset=iEntry*byteLen;
      for(size_t ibyte=0; ibyte<byteLen; ibyte++ )
      {
          data[byte_offset+ibyte]=value[ibyte];
      }
      
      switch(_info.TYPE)
      {
        case VX_FLOAT:
        {
          double restmp=0;
          memcpy (&restmp, value, byteLen);      
          long long int res=static_cast<long long int>(restmp);
          if(res>0)
          {
            nzeroEntryIndexes.insert(iEntry);
          }
          break;
        }
        case VX_FIXED:
        {
          long long int res=0;
          memcpy (&res, value, byteLen);      
          if(res>0)
          {
              nzeroEntryIndexes.insert(iEntry);
          }
          break;
        }
        case VX_UFIXED:
        {
          long long unsigned res=0;
          memcpy (&res, value, byteLen);      
          if(res>0)
          {
              nzeroEntryIndexes.insert(iEntry);
          }
          break;
        }
        default:
        {
          std::cerr<<"unknown type"<<std::endl;
          exit(1);
          break;
        }
      }
  }
  if(readerisok)
  {
      std::cout<<"...done; there are "<<nzeroEntryIndexes.size()<<" non-zero pixels"<<std::endl;
  }
  
  delete [] value;
  value=NULL;
  return(readerisok);
}


size_t INRreader::index(const size_t & ix,const size_t & iy,const size_t & iz,const size_t & iv)
{
  //Returns the index within the data array. The figure is usually written as follows:
  // Each component is composed by nb_of_pixels data; so (nb_of_pixels)_comp1,..,(nb_of_pixels)compNV
  // Pixels are divided by plane (z direction); so (px_per_Plane)_plane1,...,(px_per_Plane)_planeNZ
  // For each plane, the first ny component i read will constitute the first colum and so on: it is written column_wise
  size_t _index=0;
  size_t vdimcompStart=(iv*nb_Of_Pixels);
  size_t zPlaneStart= px_per_Plane*iz;
  _index=vdimcompStart+zPlaneStart+(iy*_info.SHAPE[0])+ix  ;
  return(_index);
}


IndexCoord INRreader::reverseIndex(const size_t & index )
{
  IndexCoord reverseIndexing;
  // first the vector component;
  reverseIndexing.iv=std::floor(index/nb_Of_Pixels);
  
  //now the plane
  size_t remainder=index-(reverseIndexing.iv*nb_Of_Pixels);
  reverseIndexing.iz=std::floor(remainder/px_per_Plane);
  
  //now the rows
  size_t remainder2=remainder-(reverseIndexing.iz*px_per_Plane);
  reverseIndexing.iy=std::floor(remainder2/(_info.SHAPE[0]));
  //now the colums
  reverseIndexing.ix=remainder2-(reverseIndexing.iy*_info.SHAPE[0]);
  
  return(reverseIndexing);
}

std::vector<double> INRreader::evalBarycenter(const size_t & index)
{
    std::vector<double> barycenter(3,0.0);
    IndexCoord intcoord=reverseIndex(index );
    barycenter[0]=(static_cast<double>(intcoord.ix)+0.5)*_info.RESOLUTION[0];
    barycenter[1]=(static_cast<double>(intcoord.iy)+0.5)*_info.RESOLUTION[1];
    barycenter[2]=(static_cast<double>(intcoord.iz)+0.5)*_info.RESOLUTION[2];
    return(barycenter);
}

double INRreader::EuclideanDist(const double * P1, const double * P2)
{
    double dist=std::sqrt((P1[0]-P2[0])*(P1[0]-P2[0])+(P1[1]-P2[1])*(P1[1]-P2[1])+(P1[2]-P2[2])*(P1[2]-P2[2]));
    return(dist);               
}

//used to evaluate endianess
bool INRreader::isLittleEndian()
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
void INRreader::SwapBytes(void *pv, size_t n)
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



