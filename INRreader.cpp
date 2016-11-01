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
#include<map>

#include "INRreader.hpp"
#ifndef DELTAMAX
#define DELTAMAX 5
#endif

INRreader::INRreader(std::string filename)
:_isAllocated(false),
nb_Of_Pixels(0),
px_per_Plane(0),
byteLen(0),
data(NULL)
{
    nzeroEntryIndexes.clear();
    readSegmentation(filename);
    _bboxlabels.clear();
}

INRreader::INRreader()
:_isAllocated(false),
nb_Of_Pixels(0),
px_per_Plane(0),
byteLen(0),
data(NULL)
{
    nzeroEntryIndexes.clear();
    _bboxlabels.clear();
}

INRreader::~INRreader()
{
  freeMemory();
}

void INRreader::freeMemory()
{
  if(_isAllocated)
  {
    _isAllocated=false;
    delete [] data;
    data=NULL;
    byteLen=0;
    nzeroEntryIndexes.clear();
    _bboxlabels.clear();
    nb_Of_Pixels=0;
    px_per_Plane=0;
    _info.freeMemory();
  }

}

void INRreader::readSegmentation(const std::string & filename)
{
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
    _isAllocated=true;
}

void INRreader::createBoundingBoxes()
{
  if(_isAllocated)
  {
      evalLabeledRegionsBounds();
  }
}

bool INRreader::isPointInsideSegmentation(const double & x, const double & y, const double & z) const
{
  bool is_inside=false;
  if(_isAllocated)
  {
    IndexCoord Ixyz=voxelCoordInterp(x,y,z);
    for(size_t iv=0; iv<_info.VDIM; iv++)
    {
        double vox=pickVoxelValue(Ixyz.ix,Ixyz.iy,Ixyz.iz,iv);
        long long int voxval=static_cast<long long int>(vox);
        if(voxval>0)
        {
          is_inside=true;
          break;
        }
    }
  }
  return(is_inside);
}

IndexCoord INRreader::voxelCoordInterp(const double & x, const double & y, const double & z) const
{
  IndexCoord Ixyz;
  Ixyz.iv=0;
  
  if(_isAllocated)
  {
    Ixyz.ix=static_cast<size_t>(std::floor(x/_info.RESOLUTION[0]));
    Ixyz.iy=static_cast<size_t>(std::floor(y/_info.RESOLUTION[1]));
    Ixyz.iz=static_cast<size_t>(std::floor(z/_info.RESOLUTION[2]));
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
  }
  return(Ixyz);
}

IndexCoord INRreader::voxelCoordInterpNonZero(const double & x, const double & y, const double & z)
{
    IndexCoord Ixyz=voxelCoordInterp(x,y,z);
    if(_isAllocated)
    {
      bool iszero=!(isPointInsideSegmentation(x,y,z));
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
    }    
    return(Ixyz);
}

std::vector<double> INRreader::interpolatedVoxelValue (const double & x, const double & y, const double & z)
{
  std::vector<double> voxvalue(_info.VDIM,0.0);
  if(_isAllocated)
  {
    IndexCoord voxcoord=voxelCoordInterp(x,y,z);
    for(size_t iv=0; iv<_info.VDIM; iv++)
    {
        voxvalue[iv]=pickVoxelValue(voxcoord.ix,voxcoord.iy,voxcoord.iz, iv);
    }
  }
  return(voxvalue);
}

std::vector<double> INRreader::interpolatedNonZeroVoxelValue (const double & x, const double & y, const double & z)
{
  std::vector<double> voxvalue(_info.VDIM,0.0);
  if(_isAllocated)
  {
    IndexCoord nzerocoord=voxelCoordInterpNonZero(x,y,z);
    for(size_t iv=0; iv<_info.VDIM; iv++)
    {
        voxvalue[iv]=pickVoxelValue(nzerocoord.ix,nzerocoord.iy,nzerocoord.iz, iv);
    }
  }
  return(voxvalue);
}

double INRreader::pickVoxelValue(const size_t & ix,const size_t & iy,const size_t & iz,const size_t iv) const
{
    return(pickValue(index(ix,iy,iz, iv)));
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

void INRreader::evalLabeledRegionsBounds()
{
  // Evaluates bounding box of regions with label different from 0 and 1;
  // Implemented for VDIM=1 only
  // _bboxlabels is a map that, for each region label different from 1 assign a bounding box
  _bboxlabels.clear();
  if(_info.VDIM == 1)
  {
    typedef long long int regionLabeltype;
    
    typedef std::set<size_t> setPointType;
    typedef setPointType::iterator setPointTypeIterator;
    typedef std::set<regionLabeltype> regionSetType;
    
    typedef std::map<size_t, regionLabeltype > voxelToRegionMapType;
    typedef std::map<regionLabeltype, setPointType> regionSubdivisionType;
    typedef regionSubdivisionType::iterator regionSubdivisionTypeIterator;
    typedef std::map<size_t,setPointType> voxConnectType;

    regionSubdivisionType voxelRegions;
    regionSetType regionLabels;
    
    //voxelToregionMap is a map between a  voxel and its regions
    voxelToRegionMapType voxelToregionMap;
    
    // Initialize the map setpointtype that describes the set of indices with the same label
    for (setPointTypeIterator it=nzeroEntryIndexes.begin(); it!=nzeroEntryIndexes.end(); ++it)
    {
      regionLabeltype voxLabel=static_cast<regionLabeltype >(pickValue(*it));
      voxelRegions[voxLabel].insert(*it);
    }
    // Delete elements with index label =1
    voxelRegions.erase(voxelRegions.find(1));
    // create a set of only non-zero and non-one elements
    setPointType boundVoxels;
    boundVoxels.clear();
    regionLabels.clear();
    for(regionSubdivisionTypeIterator mapit=voxelRegions.begin(); mapit!=voxelRegions.end();++mapit)
    {
      regionLabels.insert(mapit->first);
      for(setPointTypeIterator itset=((mapit->second).begin());itset!=((mapit->second).end());++itset)
      {
        boundVoxels.insert(*itset);
        voxelToregionMap.insert(std::pair<size_t,regionLabeltype>(*itset, mapit->first) );
      }
    }
    // Now I create a voxel connectivity
    voxConnectType connectivity;
    for(setPointTypeIterator it=boundVoxels.begin(); it!=boundVoxels.end(); ++it)
    {
      IndexCoord Ixyz=reverseIndex(*it);
      Ixyz.iv=0;
      if(Ixyz.ix>0)
      {
        size_t jind=index((Ixyz.ix-1),Ixyz.iy,Ixyz.iz,Ixyz.iv);
        if(boundVoxels.find(jind)!=boundVoxels.end())
        {
          connectivity[*it].insert(jind);
          connectivity[jind].insert(*it);
        }
        
        if(Ixyz.iy>0)
        {
          jind=index((Ixyz.ix-1),(Ixyz.iy-1),Ixyz.iz,Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          if(Ixyz.iz>0)
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy-1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iz<(_info.SHAPE[2]-1))
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy-1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
        }// and on xy, y>0
        
        
        if(Ixyz.iy<(_info.SHAPE[1]-1))
        {
          jind=index((Ixyz.ix-1),(Ixyz.iy+1),Ixyz.iz,Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          if(Ixyz.iz>0)
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy+1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iz<(_info.SHAPE[2]-1))
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy+1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
        }// and on xy, y<NYDIM
        
        if(Ixyz.iz>0)
        {
          jind=index((Ixyz.ix-1),Ixyz.iy,(Ixyz.iz-1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          
          /*if(Ixyz.iy>0)
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy-1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iy<(_info.SHAPE[1]-1))
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy+1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }*/
        }//and of xz, z>0
        
        
        if(Ixyz.iz<(_info.SHAPE[2]-1))
        {
          jind=index((Ixyz.ix-1),Ixyz.iy,(Ixyz.iz+1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          /*if(Ixyz.iy>0)
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy-1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iy<(_info.SHAPE[1]-1))
          {
            jind=index((Ixyz.ix-1),(Ixyz.iy+1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }*/
        }//and of xz, z<NZDIM
      } //(end of ix >0)


      if(Ixyz.ix<(_info.SHAPE[0]-1))
      {
        size_t jind=index((Ixyz.ix+1),Ixyz.iy,Ixyz.iz,Ixyz.iv);
        if(boundVoxels.find(jind)!=boundVoxels.end())
        {
          connectivity[*it].insert(jind);
          connectivity[jind].insert(*it);
        }
        
        if(Ixyz.iy>0)
        {
          jind=index((Ixyz.ix+1),(Ixyz.iy-1),Ixyz.iz,Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          if(Ixyz.iz>0)
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy-1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iz<(_info.SHAPE[2]-1))
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy-1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
        }// and on xy, y>0
        
        
        if(Ixyz.iy<(_info.SHAPE[1]-1))
        {
          jind=index((Ixyz.ix+1),(Ixyz.iy+1),Ixyz.iz,Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          if(Ixyz.iz>0)
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy+1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iz<(_info.SHAPE[2]-1))
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy+1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
        }// and on xy, y<NYDIM
        
        if(Ixyz.iz>0)
        {
          jind=index((Ixyz.ix+1),Ixyz.iy,(Ixyz.iz-1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          
          /*if(Ixyz.iy>0)
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy-1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iy<(_info.SHAPE[1]-1))
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy+1),(Ixyz.iz-1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }*/
        }//and of xz, z>0
        
        if(Ixyz.iz<(_info.SHAPE[2]-1))
        {
          jind=index((Ixyz.ix+1),Ixyz.iy,(Ixyz.iz+1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
          /*if(Ixyz.iy>0)
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy-1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }
          if(Ixyz.iy<(_info.SHAPE[1]-1))
          {
            jind=index((Ixyz.ix+1),(Ixyz.iy+1),(Ixyz.iz+1),Ixyz.iv);
            if(boundVoxels.find(jind)!=boundVoxels.end())
            {
              connectivity[*it].insert(jind);
              connectivity[jind].insert(*it);
            }
          }*/
        }//and of xz, z<NZDIM
      }// (end of ix<NXDIM)
      
      
      if(Ixyz.iy>0)
      {
        size_t jind=index(Ixyz.ix,(Ixyz.iy-1),Ixyz.iz,Ixyz.iv);
        if(boundVoxels.find(jind)!=boundVoxels.end())
        {
          connectivity[*it].insert(jind);
          connectivity[jind].insert(*it);
        }
        
        if(Ixyz.iz>0)
        {
          jind=index(Ixyz.ix,(Ixyz.iy-1),(Ixyz.iz-1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
        }//end on iy, iz>0
        
        if(Ixyz.iz<(_info.SHAPE[2]-1))
        {
          jind=index(Ixyz.ix,(Ixyz.iy-1),(Ixyz.iz+1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
        }//end on iy, iz<NZDIM
      }//end on iy>0
      

      if(Ixyz.iy<(_info.SHAPE[1]-1))
      {
        size_t jind=index(Ixyz.ix,(Ixyz.iy+1),Ixyz.iz,Ixyz.iv);
        if(boundVoxels.find(jind)!=boundVoxels.end())
        {
          connectivity[*it].insert(jind);
          connectivity[jind].insert(*it);
        }

        if(Ixyz.iz>0)
        {
          jind=index(Ixyz.ix,(Ixyz.iy+1),(Ixyz.iz-1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
        }//end on iy, iz>0
        
        if(Ixyz.iz<(_info.SHAPE[2]-1))
        {
          jind=index(Ixyz.ix,(Ixyz.iy+1),(Ixyz.iz+1),Ixyz.iv);
          if(boundVoxels.find(jind)!=boundVoxels.end())
          {
            connectivity[*it].insert(jind);
            connectivity[jind].insert(*it);
          }
        }//end on iy, iz<NZDIM
      }// end on iy<NYDIM
      
      if(Ixyz.iz>0)
      {
        size_t jind=index(Ixyz.ix,Ixyz.iy,(Ixyz.iz-1),Ixyz.iv);
        if(boundVoxels.find(jind)!=boundVoxels.end())
        {
          connectivity[*it].insert(jind);
          connectivity[jind].insert(*it);
        }
      }//end on iz>0
      if(Ixyz.iz<(_info.SHAPE[2]-1))
      {
        size_t jind=index(Ixyz.ix,Ixyz.iy,(Ixyz.iz+1),Ixyz.iv);
        if(boundVoxels.find(jind)!=boundVoxels.end())
        {
          connectivity[*it].insert(jind);
          connectivity[jind].insert(*it);
        }
      }// end on iz<NZDIM
    }
    
    // At this point, I have a voxel connectivity structure of the
    // Domain items that falls on the PVs and the MV
    // Now the growing algorithm on labeled regions only
    for(regionSubdivisionTypeIterator itRegToNodeMap=voxelRegions.begin(); itRegToNodeMap!=voxelRegions.end(); ++itRegToNodeMap)
    {
      regionLabeltype labelOfRegions=itRegToNodeMap->first;
      // iterate on voxels belonging to the region labelOfRegions
      // extract the local connectivity of the region regionconnect
      voxConnectType regionconnect;
      regionconnect.clear();
      for(setPointTypeIterator cRegIt=(itRegToNodeMap->second).begin(); cRegIt!=(itRegToNodeMap->second).end(); ++cRegIt)
      {
        // node belonging to region
        size_t voxel=*cRegIt;
        // explore connectivity of node; extract only connections belonging to the same region
        setPointType connections;
        for(setPointTypeIterator connectiter=connectivity[voxel].begin(); connectiter!=connectivity[voxel].end(); ++connectiter)
        {
          //voxelToregionMap[*connectiter]: set<long long int> of region labels
          // connectiter belongs to the region under study
          if((voxelToregionMap[*connectiter])==labelOfRegions)
          {
            connections.insert(*connectiter);
          }
        } //end loop with connectiter
        regionconnect.insert(std::pair<size_t, setPointType> (voxel,connections) );
      }// end loop with cRegIt iterator
      //first: determine a seed point
      size_t seed= *((itRegToNodeMap->second).begin());

      setPointType RegionVoxels, Queue;
      Queue.insert(seed);
      RegionVoxels.insert(seed);
      /* Here the algorithm grows: expands the region RegionVoxels using the regional connectivity
       it iterates until all the connected voxels are covered
       */
      while(!Queue.empty())
      {
        size_t candidate= *(Queue.begin());
        Queue.erase(candidate);
        //regionconnect : map<size_t, connectSetType> isw the local connectivity on the region
        for(setPointTypeIterator countVoxel=regionconnect.at(candidate).begin(); countVoxel!=regionconnect.at(candidate).end(); ++countVoxel)
        {
          //check if countNode is inside the region
          if(RegionVoxels.find(*countVoxel)==RegionVoxels.end())
          {
            //if not member of nodesRegion, add it
            RegionVoxels.insert(*countVoxel);
            Queue.insert(*countVoxel);
          }
        }
      }//end of while
      //RegionVoxels: voxels of the connected region with label itRegToNodeMap->first
      //check if ... itRegToNodeMap->second : list of points belonging to (set)
      //if is the case, there are some point to move
      if((itRegToNodeMap->second).size() !=RegionVoxels.size()  )
      {
        //first: create a new label region
        regionLabeltype newRegionLabel=1+(*(regionLabels.rbegin()));
        regionLabels.insert(newRegionLabel);

        //copy inside newSet the whole set of point belonging to the current region
        setPointType newSet=(itRegToNodeMap->second);
        
        //Delete points identified to the current region
        for(setPointTypeIterator reg_iter=RegionVoxels.begin();reg_iter!=RegionVoxels.end(); ++reg_iter)
        {
          newSet.erase(*reg_iter);
        }
        //insert the new region inside pointRegions
        voxelRegions.insert(std::pair<regionLabeltype, setPointType>(newRegionLabel,newSet) );
        //nodeToRegionMapType = <size_t, set<long long int> >
        //remove point not belonging to region
        for(setPointTypeIterator reg_iter=newSet.begin();reg_iter!=newSet.end(); ++reg_iter)
        {
          (itRegToNodeMap->second).erase(*reg_iter);
          voxelToregionMap[*reg_iter]=newRegionLabel;
        }
      }// end if on size
    }//algo division ends
    //now re-labeling of voxels
    for(regionSubdivisionTypeIterator itRegToNodeMap=voxelRegions.begin(); itRegToNodeMap!=voxelRegions.end(); ++itRegToNodeMap)
    {
       double voxValue = static_cast<double>(itRegToNodeMap->first);
       BoundingBox localbbox;
       for(unsigned char jcoord=0; jcoord<3; jcoord++)
       {
        ( localbbox.bbox()[jcoord] )[0]=2.0*_info.RESOLUTION[jcoord]*_info.SHAPE[jcoord];
        ( localbbox.bbox()[jcoord] )[1]=0.0;
       }
       for(setPointTypeIterator vox_iter=(itRegToNodeMap->second).begin();vox_iter!=(itRegToNodeMap->second).end(); ++vox_iter)
       {
         setValue(*vox_iter, voxValue);
         std::vector<double>  barycenter=evalBarycenter(*vox_iter);
         for(unsigned char jcoord=0; jcoord<3; jcoord++)
         {
            if(barycenter[jcoord]< ( (localbbox.bbox()[jcoord])[0] )  )
            {
                (localbbox.bbox()[jcoord])[0]=barycenter[jcoord];
            }
            
            if(barycenter[jcoord]> ( (localbbox.bbox()[jcoord])[1] ) )
            {
                (localbbox.bbox()[jcoord])[1]=barycenter[jcoord];
            }
         }
       }
       for(unsigned char jcoord=0; jcoord<3; jcoord++)
       {
         (localbbox.bbox()[jcoord])[0]=(localbbox.bbox()[jcoord])[0]-0.5*_info.RESOLUTION[jcoord];
         (localbbox.bbox()[jcoord])[1]=(localbbox.bbox()[jcoord])[1]+0.5*_info.RESOLUTION[jcoord];
       }

       _bboxlabels.insert(std::pair<double,BoundingBox >(voxValue,localbbox) );
    }
    

    for(regionSubdivisionTypeIterator mapit=voxelRegions.begin();mapit!=voxelRegions.end();++mapit)
    {
        LinearRegression2DData regdata=evalRegressionPlane(mapit->second);
        if(regdata.isAhole)
        {
          std::vector<double> pnorm=regdata.pnorm;
          double norm=sqrt(pnorm[0]*pnorm[0]+pnorm[1]*pnorm[1]+pnorm[2]*pnorm[2]);
          for(unsigned char cc=0; cc<3; cc++)
          {
              pnorm[cc]=pnorm[cc]/norm;
          }
          std::vector<double> xc=regdata.xc;
          for(std::set<size_t>::iterator it=nzeroEntryIndexes.begin();it!=nzeroEntryIndexes.end(); ++it)
          {
              std::vector<double> barp=evalBarycenter(*it);
              for(unsigned char cc=0; cc<3; cc++)
              {
                  barp[cc]=barp[cc]-xc[cc];
              }
              double dproj=dotprod(barp,pnorm);
              
              if((dproj>=regdata.zminmax[0]) &&  (dproj<=regdata.zminmax[1]) )
              {
                //projected on plane
                for(unsigned char cc=0; cc<3; cc++)
                {
                    barp[cc]=barp[cc]-dproj*pnorm[cc];
                }
                double P_xc=EuclideanDist(barp.data(), xc.data());
                //if point voxel is between the two planes zmin and zmax, mark the region
                if(P_xc<=regdata.rhominmax[1])
                {
                  double voxValue = static_cast<double>(mapit->first);
                  setValue(*it, voxValue);
                }
              } 
          }
        }
    }
  }
}

double INRreader::pickValue(const size_t & _index) const
{
  
  double res=0.0;
  if(_isAllocated)
  {
    char * value= new char[byteLen];
    size_t byte_offset=byteLen*_index;
    for(size_t ibyte=0; ibyte<byteLen; ibyte++ )
    {
      value[ibyte]=data[byte_offset+ibyte];
    }
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
  }
  return(res);
}

void INRreader::setValue(const size_t & _index, const double & _value)
{
  char * buffer= new char[byteLen];
  switch(_info.TYPE)
  {
    case VX_FLOAT:
    {
      memcpy (buffer, &_value, byteLen);
      break;
    }
    case VX_FIXED:
    {
      long long int lintres=static_cast<long long int>(_value);
      memcpy ( buffer, &lintres, byteLen);
      break;
    }
    case VX_UFIXED:
    {
      long long unsigned luires=static_cast<long long unsigned>(_value);
      memcpy ( buffer, &luires, byteLen);
      break;
    }
    default:
    {
      std::cerr<<"unknown type"<<std::endl;
      exit(1);
      break;
    }
  }
  size_t byte_offset=byteLen*_index;
  for(size_t ibyte=0; ibyte<byteLen; ibyte++ )
  {
    data[byte_offset+ibyte]=buffer[ibyte];
  }
  delete [] buffer;
  buffer=NULL;
}

size_t INRreader::index(const size_t & ix,const size_t & iy,const size_t & iz,const size_t & iv) const
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

IndexCoord INRreader::reverseIndex(const size_t & index ) const
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

double INRreader::dotprod(const std::vector<double> & v1, const std::vector<double> & v2)
{
  double prod=0.0;
  for(size_t ii=0; ii<v1.size(); ii++)
  {
      prod=prod+v1[ii]*v2[ii];
  }
  return(prod);
}

LinearRegression2DData INRreader::evalRegressionPlane(std::set<size_t> & pointlist)
{
  LinearRegression2DData datareg;
  size_t n=pointlist.size();
  if(n>=3)
  {
    size_t counter=0;
    std::vector<std::vector<double> > pointcoord(n);
    std::vector<double> r_G(3,0);
    for(std::set<size_t>::iterator it=pointlist.begin(); it!=pointlist.end();++it)
    {
      std::vector<double> bar=evalBarycenter(*it);
      pointcoord[counter]=bar;
      for(unsigned char cc=0; cc<3; cc++)
      {
        r_G[cc]=r_G[cc]+bar[cc];
      }
      counter++;
    }
    for(unsigned char cc=0; cc<3; cc++)
    {
      r_G[cc]=r_G[cc]/static_cast<double>(n);
      for(size_t pt=0; pt<n; pt++)
      {
        (pointcoord[pt])[cc]=(pointcoord[pt])[cc]-r_G[cc];
      }
    }
    datareg.xc=r_G;
    std::vector<double> LsqMat(9,0);
    for(unsigned char ic=0; ic<3; ic++)
    {
      double matrixEntry=0.0;
      for(unsigned char jc=ic; jc<3; jc++)
      {
        for(size_t pt=0; pt<n; pt++)
        {
          matrixEntry=matrixEntry+(pointcoord[pt])[ic]*(pointcoord[pt])[jc];
        }
        LsqMat[RM3X3Ind(ic,jc)]=matrixEntry;
        LsqMat[RM3X3Ind(jc,ic)]=matrixEntry;
      }
    }
    
    eigenData eigens=eigen_decomposition3X3(LsqMat);
    for(unsigned char ic=0; ic<3; ic++)
    {
      datareg.pnorm[ic]=eigens.eigenVectors[RM3X3Ind(ic,0)]; //extract the 1st col
    }

    
    std::vector<double> normdir=datareg.pnorm;
    double normnorm=sqrt(normdir[0]*normdir[0]+normdir[1]*normdir[1]+normdir[2]*normdir[2]);
    for(unsigned char ic=0; ic<3; ic++)
    {
      normdir[ic]=normdir[ic]/normnorm;
    }
    
    std::vector<double> signed_dist_from_plane(n,0.0);
    std::vector<double> proj_dist_from_center(n,0.0);
    for(size_t ipt=0; ipt<n; ipt++)
    {
      // dot product brtween normal versor and coordinates
      double tmp=0.0;
      std::vector<double> localcoord=(pointcoord[ipt]);
      for(unsigned char ic=0; ic<3; ic++)
      {
        tmp=tmp+normdir[ic]*localcoord[ic];
      }
      signed_dist_from_plane[ipt]=tmp;
      // Plane projection
      tmp=0.0;
      for(unsigned char ic=0; ic<3; ic++)
      {
        localcoord[ic]=localcoord[ic]-signed_dist_from_plane[ipt]*normdir[ic];
        tmp=tmp+localcoord[ic]*localcoord[ic];
      }
      proj_dist_from_center[ipt]=sqrt(tmp);
    }
    // proj_dist_from_center: in-plane circular distribution of the points
    // signed_dist_from_plane point locations above/below the plane
    double zmin=1e32,zmax=-1e3,rhomin=zmin,rhomax=zmax;
    for(size_t ipt=0; ipt<n; ipt++)
    {
      if(signed_dist_from_plane[ipt]>zmax)
      {
        zmax=signed_dist_from_plane[ipt];
      }
      if(signed_dist_from_plane[ipt]<zmin)
      {
        zmin=signed_dist_from_plane[ipt];
      }

      if(proj_dist_from_center[ipt]>rhomax)
      {
        rhomax=proj_dist_from_center[ipt];
      }
      if(proj_dist_from_center[ipt]<rhomin)
      {
        rhomin=proj_dist_from_center[ipt];
      }
    }
    datareg.zminmax[0]=zmin;
    datareg.zminmax[1]=zmax;
    datareg.rhominmax[0]=rhomin;
    datareg.rhominmax[1]=rhomax;
    datareg.isAhole=false;
    if(!isPointInsideSegmentation(datareg.xc[0],datareg.xc[1],datareg.xc[2]))
    {
      datareg.isAhole=true;
    }

  }
  return datareg;
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

eigenData INRreader::eigen_decomposition3X3(std::vector<double> &Amatr)
{
  eigenData eigens;
  std::vector<double> V=Amatr;
  std::vector<double> e(3,0.0), d(3,0.0);
  for (unsigned char  cc = 0; cc < 3; cc++)
  {
    d[cc] = V[RM3X3Ind(2,cc)];
  }
  // Householder reduction to tridiagonal form.
  for (unsigned char cc = 1; cc<2; cc++)
  {
    unsigned char iindex=3-cc;
    // Scale to avoid under/overflow.
    double scale = 0.0;
    double h = 0.0;
    for (unsigned char k = 0; k < iindex; k++)
    {
      scale = scale + std::sqrt(d[k]*d[k]);
    }
    if (scale == 0.0)
    {
      e[iindex] = d[iindex-1];
      for (unsigned char j = 0; j < iindex; j++)
      {
        d[j] = V[RM3X3Ind(iindex-1,j)];
        V[RM3X3Ind(iindex,j)]  = 0.0;
        V[RM3X3Ind(j,iindex)]  = 0.0;
      }
    }
    else
    {
      // Generate Householder vector.
      for (unsigned char k = 0; k < iindex; k++)
      {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[iindex-1];
      double g = sqrt(h);
      if (f > 0)
      {
        g = -g;
      }
      e[iindex] = scale * g;
      h = h - f * g;
      d[iindex-1] = f - g;
      for (unsigned char j = 0; j < iindex; j++)
      {
        e[j] = 0.0;
      }
      // Apply similarity transformation to remaining columns.
      for (unsigned char j = 0; j < iindex; j++)
      {
        f = d[j];
        V[RM3X3Ind(j,iindex)] = f;
        g = e[j] + V[RM3X3Ind(j,j)] * f;
        for (unsigned char k = j+1; k <= iindex-1; k++)
        {
          g += V[RM3X3Ind(k,j)] * d[k];
          e[k] += V[RM3X3Ind(k,j)] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (unsigned char j = 0; j < iindex; j++)
      {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for (unsigned char j = 0; j < iindex; j++)
      {
        e[j] -= hh * d[j];
      }
      for (unsigned char j = 0; j < iindex; j++)
      {
        f = d[j];
        g = e[j];
        for (unsigned char k = j; k <= iindex-1; k++)
        {
          V[RM3X3Ind(k,j)] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[RM3X3Ind(iindex-1,j)];
        V[RM3X3Ind(iindex,j)] = 0.0;
      }
    }
    d[iindex] = h;
  }
  // Accumulate transformations.
  for (unsigned char i = 0; i < 2; i++)
  {
    V[RM3X3Ind(2,i)]=V[RM3X3Ind(i,i)];
    V[RM3X3Ind(i,i)] = 1.0;
    double h = d[i+1];
    if (h != 0.0)
    {
      for (unsigned char k = 0; k <= i; k++)
      {
        d[k] = V[RM3X3Ind(k,i+1)]/ h;
      }
      for (unsigned char j = 0; j <= i; j++)
      {
        double g = 0.0;
        for (unsigned char k = 0; k <= i; k++)
        {
          g += V[RM3X3Ind(k,i+1)]*V[RM3X3Ind(k,j)];
        }
        for (unsigned char k = 0; k <= i; k++)
        {
          V[RM3X3Ind(k,j)] -= g * d[k];
        }
      }
    }
    for (unsigned char k = 0; k <= i; k++)
    {
      V[RM3X3Ind(k,i+1)] = 0.0;
    }
  }
  for (unsigned char j = 0; j < 3; j++)
  {
    d[j] = V[RM3X3Ind(2,j)];
    V[RM3X3Ind(2,j)] = 0.0;
  }
  V[RM3X3Ind(2,2)] = 1.0;
  e[0] = 0.0;
  for (unsigned char i = 1; i < 3; i++)
  {
    e[i-1] = e[i];
  }
  e[2] = 0.0;
  double f = 0.0;
  double tst1 = 0.0;
  double eps = 1.0e-52;
  for (unsigned char l = 0; l < 3; l++)
  {
    // Find small subdiagonal element
    double tmp=std::sqrt(d[l]*d[l])+std::sqrt(e[l]*e[l]);
    tst1 = std::max(tst1,tmp);
    unsigned char m = l;
    while (m < 2)
    {
      if (std::sqrt(e[m]*e[m]) <= eps*tst1)
      {
        break;
      }
      m++;
    }
    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.
    if (m > l)
    {
      int iter = 0;
      do
      {
        iter = iter + 1;  // (Could check iteration count here.)
        if(iter>500)
        {
          break;
        }
        // Compute implicit shift
        double g = d[l];
        double p = (d[l+1] - g) / (2.0 * e[l]);
        double r = hypot2(p,1.0);
        if (p < 0)
        {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        double dl1 = d[l+1];
        double h = g - d[l];
        for (unsigned char i = l+2; i < 3; i++)
        {
          d[i] -= h;
        }
        f = f + h;
        // Implicit QL transformation.
        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l+1];
        double s = 0.0;
        double s2 = 0.0;
        for (signed char i = m-1; i >= l; i--)
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);
          
          // Accumulate transformation.
          
          for (unsigned char k = 0; k < 3; k++)
          {
            h = V[RM3X3Ind(k,i+1)];
            V[RM3X3Ind(k,i+1)] = s * V[RM3X3Ind(k,i)] + c * h;
            V[RM3X3Ind(k,i)] = c * V[RM3X3Ind(k,i)] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        // Check for convergence.
      } while (std::sqrt(e[l]*e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  // Sort eigenvalues and corresponding vectors.
  for (unsigned char i = 0; i < 2; i++)
  {
    unsigned char k = i;
    double p = d[i];
    for (unsigned char j = i+1; j < 3; j++)
    {
      if (d[j] < p)
      {
        k = j;
        p = d[j];
      }
    }
    if (k != i)
    {
      d[k] = d[i];
      d[i] = p;
      for (unsigned char j = 0; j < 3; j++)
      {
        p = V[RM3X3Ind(j,i)];
        V[RM3X3Ind(j,i)] = V[RM3X3Ind(j,k)];
        V[RM3X3Ind(j,k)] = p;
      }
    }
  }
  eigens.eigenValues=d;
  eigens.eigenVectors=V;
  return(eigens);
  
}

std::vector<double> INRreader::InvertA3X3(const std::vector<double> & Mat0) const
{
  //matrix is considered as row major
  std::vector<double> inverted(9,0);
  double det =    Mat0[RM3X3Ind(0,0)]*( Mat0[RM3X3Ind(1,1)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(1,2)]*Mat0[RM3X3Ind(2,1)] ) +
  -1.0*Mat0[RM3X3Ind(0,1)]*( Mat0[RM3X3Ind(1,0)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(1,2)]*Mat0[RM3X3Ind(2,0)] ) +
  Mat0[RM3X3Ind(0,2)]*( Mat0[RM3X3Ind(1,0)]*Mat0[RM3X3Ind(2,1)] - Mat0[RM3X3Ind(2,0)]*Mat0[RM3X3Ind(1,1)] ) ;
  
  if(sqrt(det*det)>0.0)
  {
    inverted[RM3X3Ind(0,0)] =        Mat0[RM3X3Ind(1,1)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(1,2)]*Mat0[RM3X3Ind(2,1)];
    inverted[RM3X3Ind(0,1)] = -1.0*( Mat0[RM3X3Ind(0,1)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(2,1)]*Mat0[RM3X3Ind(0,2)]);
    inverted[RM3X3Ind(0,2)] =        Mat0[RM3X3Ind(0,1)]*Mat0[RM3X3Ind(1,2)] - Mat0[RM3X3Ind(0,2)]*Mat0[RM3X3Ind(1,1)];
    
    inverted[RM3X3Ind(1,0)] = -1.0*( Mat0[RM3X3Ind(1,0)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(1,2)]*Mat0[RM3X3Ind(2,0)] );
    inverted[RM3X3Ind(1,1)] =        Mat0[RM3X3Ind(0,0)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(0,2)]*Mat0[RM3X3Ind(2,0)];
    inverted[RM3X3Ind(1,2)] = -1.0*( Mat0[RM3X3Ind(0,0)]*Mat0[RM3X3Ind(1,2)] - Mat0[RM3X3Ind(0,2)]*Mat0[RM3X3Ind(1,0)]);
    
    inverted[RM3X3Ind(2,0)] =        Mat0[RM3X3Ind(1,0)]*Mat0[RM3X3Ind(2,1)] - Mat0[RM3X3Ind(2,0)]*Mat0[RM3X3Ind(1,1)];
    inverted[RM3X3Ind(2,1)] = -1.0*( Mat0[RM3X3Ind(0,0)]*Mat0[RM3X3Ind(2,1)] - Mat0[RM3X3Ind(2,0)]*Mat0[RM3X3Ind(0,1)]);
    inverted[RM3X3Ind(2,2)] =        Mat0[RM3X3Ind(1,1)]*Mat0[RM3X3Ind(0,0)] - Mat0[RM3X3Ind(1,0)]*Mat0[RM3X3Ind(0,1)];
    
    std::vector<double>::iterator it;
    for(it=inverted.begin(); it!=inverted.end(); ++it)
    {
      *it=*it/det;
    }
  }
  else
  {
    inverted.clear();
    inverted.resize(9,0);
  }
  return(inverted);
}

unsigned char INRreader::RM3X3Ind(const unsigned char & irow, const unsigned char & jcol) const
{
  unsigned char index=3*irow + jcol;
  return(index);
}

