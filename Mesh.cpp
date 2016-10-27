//  -*- c++ -*-
//  Mesh.cpp version 1.0                                     Oct/27/2016
//
//
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//
//  Implementation of the Mesh class
//
//  Developer: Cesare Corrado cesare.corrado@kcl.ac.uk
//==========================================================================

#include "Mesh.hpp"
#include "VtkWriter.hpp"
#include<string>
#include<iostream>
#include<iomanip>
#include<limits>
#include<fstream>
#include<sstream>
#include<stdio.h>
#include<cstdlib>
#include<cmath>
#include<map>
#include<set>

#ifndef HEADER_SIZE
#define HEADER_SIZE 1024
#endif

#ifndef TREECSIZE
#define TREECSIZE 201
#endif

Mesh::Mesh()
:consistentState(false),
points(0),
triangles(0),
tetrahedra(0),
triaToTet(0),
_endoLabel(-1),
_epiLabel(-1)
{
  _faceToFace.clear();
  _elemBoundaryTris.clear();
  _endoTris.clear();
  _epiTris.clear();
  _searchTree.clear();
}

Mesh::Mesh(const  std::string & inputFileame)
:consistentState(false),
points(0),
triangles(0),
tetrahedra(0),
triaToTet(0),
_endoLabel(-1),
_epiLabel(-1)
{
  _faceToFace.clear();
  _elemBoundaryTris.clear();
  _endoTris.clear();
  _epiTris.clear();
  _searchTree.clear();
  readFromFile(inputFileame);
}

Mesh::Mesh(const  Mesh & srcMesh)
:consistentState(false)
{
  _faceToFace.clear(); 
  _elemBoundaryTris.clear();
  _endoTris.clear();
  _epiTris.clear();
  _searchTree.clear();
  points.resize(srcMesh.nPt());
  triangles.resize(srcMesh.nTri());
  triaToTet.resize(srcMesh.nTri());
  tetrahedra.resize(srcMesh.nTet());
  for(short int jdim=0; jdim<3; jdim++)
  {
    info.baricenter[jdim]=srcMesh.Info().baricenter[jdim];
    info.checksum[jdim]=srcMesh.Info().checksum[jdim];
    for(short int ib=0; ib<2; ib++)
    {
      (info.bbox[jdim])[ib]=(srcMesh.Info().bbox[jdim])[ib];
    }
  }
  
  for(size_t iPt=0; iPt<this->nPt(); iPt++)
  {
    for(short int ic=0; ic<3; ic++)
    {
      points[iPt].coord[ic]=srcMesh.Pt(iPt).coord[ic];
    }
  }

  for(size_t iTri=0; iTri<this->nTri(); iTri++)
  {
    triangles[iTri].vertex[0]=srcMesh.Tri(iTri).vertex[0];
    triangles[iTri].vertex[1]=srcMesh.Tri(iTri).vertex[1];
    triangles[iTri].vertex[2]=srcMesh.Tri(iTri).vertex[2];
    triangles[iTri].regionLabel=srcMesh.Tri(iTri).regionLabel;
    triaToTet[iTri]=srcMesh.triaToTetMap(iTri);
  }

  for(size_t iTet=0; iTet<this->nTet(); iTet++)
  {
    tetrahedra[iTet].vertex[0]=srcMesh.Tet(iTet).vertex[0];
    tetrahedra[iTet].vertex[1]=srcMesh.Tet(iTet).vertex[1];
    tetrahedra[iTet].vertex[2]=srcMesh.Tet(iTet).vertex[2];
    tetrahedra[iTet].vertex[3]=srcMesh.Tet(iTet).vertex[3];
    tetrahedra[iTet].regionLabel=srcMesh.Tet(iTet).regionLabel;
  }
  consistentState=true;
}


void Mesh::readFromFile(const  std::string & inputFileame)
{
  std::string nodeFileName = inputFileame + ".pts";
  std::string elemFileName = inputFileame + ".elem";
  std::ifstream fnode(nodeFileName.c_str());
  if(!fnode)
  {
    std::cerr<<"ERROR: FILE "<<nodeFileName<<" DOES NOT EXIST OR IS CORRUPTED"<<std::endl;
    exit(1);
  }
  else
  {
    std::cout<<"Reading nodes"<<std::endl;
    size_t nbPt=0;
    fnode>>nbPt;
    points.resize(nbPt);
    for(size_t iPt=0; iPt<nbPt; iPt++)
    {
      for(short int ic=0; ic<3; ic++)
      {
        fnode>>points[iPt].coord[ic];
        info.baricenter[ic]=info.baricenter[ic]+(points[iPt].coord[ic]);
        info.checksum[ic]=info.checksum[ic]+static_cast<float>(points[iPt].coord[ic]);
        if(info.bbox[ic][0]>points[iPt].coord[ic])
        {
          info.bbox[ic][0]=points[iPt].coord[ic];
        }
        if(info.bbox[ic][1]<points[iPt].coord[ic])
        {
          info.bbox[ic][1]=points[iPt].coord[ic];
        }
      }
    }
    for(short int ic=0; ic<3; ic++)
    {
      info.baricenter[ic]=info.baricenter[ic]/nbPt;
    }
  }
  fnode.close();
  
  std::ifstream felem(elemFileName.c_str());
  if(felem)
  {
    std::cout<<"Reading elements"<<std::endl;
    typedef std::pair <mapfacetype::iterator, mapfacetype::iterator> range_iter_type;
    size_t nbTri=0;
    size_t nbTet=0;
    size_t nbElem=0;
    felem>>nbElem;
    triangles.resize(nbElem);
    tetrahedra.resize(nbElem);
    
    std::vector<short int * > permutation;
    permutation.resize(4);
    for(short int jface=0; jface<4; jface++)
    {
      permutation[jface] = new short int [3];
      for(short int jf=0; jf<3; jf++)
      {
        short int asum=jface+jf;
        short int rem= asum%4;
        permutation[jface][jf]=rem;
      }
    }
    
    mapfacetype bound_faces;
    for(size_t iElem=0; iElem<nbElem; iElem++)
    {
      std::string TypeOfElem;
      felem>>TypeOfElem;
      if(TypeOfElem=="Tt")
      {
        nbTet=nbTet+1;
        
        std::vector<facetype> ordered_faces;
        ordered_faces.resize(4);
        for(short int iv=0; iv<4; iv++)
        {
#ifndef NDEBUG
          felem>>tetrahedra.at(nbTet-1).vertex[iv];
#else
          felem>>tetrahedra[nbTet-1].vertex[iv];
#endif          
        }
#ifndef NDEBUG        
        felem>>tetrahedra.at(nbTet-1).regionLabel;
#else
        felem>>tetrahedra[nbTet-1].regionLabel;
#endif        

        //ordered_faces: the 4 faces with nodes ordered
        for(short int jface=0; jface<4; jface++)
        {
          for(short int jf=0; jf<3; jf++)
          {
            ordered_faces[jface].insert(tetrahedra[nbTet-1].vertex[permutation[jface][jf]]);
          }
        }
        // now I use as key the value of the first node
        // if the face is inside bound_faces, I delete it
        // if the face is not inside bound_faces I add it
        //bound_faces
        for(short int jface=0; jface<4; jface++)
        {
          size_t key=*(ordered_faces[jface].begin());
          bool value_to_insert=true;
          if(bound_faces.count(key))          
          {
              // range of elements having the same key
              range_iter_type itmapRange=bound_faces.equal_range(key);
              // scroll on faces having the same key 
              for(mapfacetype::iterator itmap=itmapRange.first; itmap!=itmapRange.second; ++itmap)
              {
                if(ordered_faces[jface]==(itmap->second).second)
                {
                  value_to_insert=false;
                  bound_faces.erase(itmap);
                  break;                
                }
              }
          }
          if(value_to_insert)
          {
            facetype face=ordered_faces[jface];
            faceLabtype labeledFace=std::make_pair(nbTet-1,face);
            bound_faces.insert(std::pair<size_t, faceLabtype > (key,labeledFace) );
          }
        }//end loop on jface
      }
      else
      {
        if(TypeOfElem=="Tr")
        {
          nbTri=nbTri+1;
          for(short int iv=0; iv<3; iv++)
          {
#ifndef NDEBUG
            felem>>triangles.at(nbTri-1).vertex[iv];
#else            
            felem>>triangles[nbTri-1].vertex[iv];
#endif
          }
#ifndef NDEBUG
          felem>>triangles.at(nbTri-1).regionLabel;
#else            
          felem>>triangles[nbTri-1].regionLabel;
#endif          
        }
        else
        {
          std::cout<<"element type "<<TypeOfElem<<"ignored"<<std::endl;
          felem.ignore(std::numeric_limits<std::streamsize>::max(), '\n' );
        }
      }
    }
    felem.close();
    if(nbTet<nbElem)
    {
      size_t delta=nbElem-nbTet;
      tetrahedra.erase(tetrahedra.begin()+nbTet,tetrahedra.begin()+nbTet+delta);
    }
    
    if(nbTri==0)
    {
      evalTriangles(bound_faces,  nbTri,false,true);
    }
    else
    {
      if(nbTri<nbElem)
      {
        size_t delta=nbElem-nbTri;
        triangles.erase(triangles.begin()+nbTri,triangles.begin()+nbTri+delta);
      }
    }
  
    for(short int jface=0; jface<4; jface++)
    {
      delete[] permutation[jface];
      permutation[jface] = NULL;
    }
    permutation.clear();
    
  }
  else
  {
    std::cerr<<"ERROR: FILE "<<elemFileName<<" DOES NOT EXIST OR IS CORRUPTED"<<std::endl;
    exit(1);
  }
  consistentState=true;
}



void Mesh::evalTriangles(mapfacetype bound_faces, size_t & nbTri,bool build_search_tree, bool outwardNormOnBoundary)
{
  /* This routine extracts the triangles (faces) on the boundary and fills the triangle 
     attribute and the triaToTet attribute, that maps the triangle to the tetrahedra the face 
     belongs to; the region label of the triangle is the same of the tetra it belongs to
     if Flags outwardNormOnBoundary is settet to True, then the triangle nodes are re-oriented
     in such a way the normal vector points outside the tetraedron
     in case build_search_tree is true, it also produces the search tree: a structure that,
     given a point with coordinate x,y,z gives a set of nearest triangle candidates.
  */
  if(build_search_tree)
  {
    _searchTree.resize(TREECSIZE*TREECSIZE*TREECSIZE);
  }
  nbTri=bound_faces.size();
  triangles.resize(nbTri);
  triaToTet.resize(nbTri);
  mapfacetype::iterator itmap;
  size_t jcount=0;  
  for(itmap=bound_faces.begin();itmap!=bound_faces.end();++itmap)
  {
    //extract the tetrahedra owners of the triangle
    const Tetrahedron & Tet=tetrahedra[(itmap->second).first];
    facetype_iter fiter;
    // extract triangle node labels
    Triangle Tria;
    unsigned char ivertex=0;
    for(fiter=(itmap->second).second.begin();fiter!=(itmap->second).second.end();++fiter)
    {
      Tria.vertex[ivertex]=*fiter;
      ivertex=ivertex+1;
    }
    Tria.regionLabel=Tet.regionLabel;
    
    if(outwardNormOnBoundary)
    {
      short int notInTriaIndex=-1;
      for(short int ivTe=0; ivTe<4; ivTe++)
      {
            notInTriaIndex=ivTe;
            for(short int ivTr=0; ivTr<3; ivTr++)
            {
              if(Tria.vertex[ivTr]==Tet.vertex[ivTe])
              {
                notInTriaIndex=-1;
                break;
              }
            }
            if(!(notInTriaIndex<0))
            {
              break;
            }
      }
      //tetVertex: [0,1,2] on triangle; 3 external to triangle
      std::vector<Point> teVertex(4);
      for(short int ivTr=0; ivTr<3; ivTr++)
      {
        for(short int jcoord=0;jcoord<3;jcoord++)
        {
          (teVertex[ivTr]).coord[jcoord]=(points[Tria.vertex[ivTr]]).coord[jcoord];
        }
      }
      std::vector<double> v1(3,0),v2(3,0),v3(3,0);
      for(short int jcoord=0;jcoord<3;jcoord++)
      {
        (teVertex[3]).coord[jcoord]=(points[Tet.vertex[notInTriaIndex]]).coord[jcoord];
        v1[jcoord] = teVertex[1].coord[jcoord]-teVertex[0].coord[jcoord];
        v2[jcoord] = teVertex[2].coord[jcoord]-teVertex[0].coord[jcoord];
        v3[jcoord] = teVertex[3].coord[jcoord]-teVertex[0].coord[jcoord];
      }
      std::vector<double> normal(3,0);
      normal[0]=v1[1]*v2[2]-v1[2]*v2[1];
      normal[1]=-1.0*(v1[0]*v2[2]-v1[2]*v2[0]);
      normal[2]=v1[0]*v2[1]-v1[1]*v2[0];
      double vectormod=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      for(short int jcoord=0;jcoord<3;jcoord++)
      {
        normal[jcoord]=normal[jcoord]/vectormod;
      }
      //now evaluate the projection onto normal.
      // if projection>0, inward normal -> change point order
      double proj=normal[0]*v3[0]+normal[1]*v3[1]+normal[2]*v3[2];
      if(proj>0)
      {
        size_t tmp= Tria.vertex[1];
        Tria.vertex[1]=Tria.vertex[2];
        Tria.vertex[2]=tmp;
      }
    }//end if on outwardNormOnBoundary
    
    if(build_search_tree)
    {
        std::vector<double > centroid=ElemCentroid(Tria);
        size_t hkey=evalHashKey(centroid);
        _searchTree[hkey].insert(jcount);
    }
    triaToTet[jcount]=(itmap->second).first;
    //copy into triangle vector
    for(ivertex=0;ivertex<3;ivertex++)
    {
      triangles[jcount].vertex[ivertex]=Tria.vertex[ivertex];
    }
    triangles[jcount].regionLabel=Tria.regionLabel;
    jcount=jcount+1;        
  }//end loop on itmap iterator

}

double Mesh::hTri(size_t iTri )
{
  //! returns the diameter of the circumscribed circle
  double h=0.0;
  if(consistentState && nTri())
  {
    double area=AreaTri(iTri);
    Triangle & Tria=triangles[iTri];
    std::vector<Point> TriaPts(3);
    for(short int iPt=0; iPt<3; iPt++)
    {
      for(short int ic=0; ic<3; ic++)
      {
        TriaPts[iPt].coord[ic]=points[Tria.vertex[iPt]].coord[ic];
      }
      //TriaPts[iPt].x=points[Tria.vertex[iPt]].x;
      //TriaPts[iPt].y=points[Tria.vertex[iPt]].y;
      //TriaPts[iPt].z=points[Tria.vertex[iPt]].z;
    }
    short int * permutation = new short int[3];
    permutation[0]=1;
    permutation[1]=2;
    permutation[2]=0;
    double egdgeproduct=1.0;
    for(short int iEdge=0; iEdge<3; iEdge++)
    {
      double dx=(TriaPts[iEdge].x()-TriaPts[permutation[iEdge]].x());
      double dy=(TriaPts[iEdge].y()-TriaPts[permutation[iEdge]].y());
      double dz=(TriaPts[iEdge].z()-TriaPts[permutation[iEdge]].z());
      egdgeproduct=egdgeproduct*sqrt(dx*dx+dy*dy+dz*dz);
    }
    delete [] permutation;  
    h=0.5*egdgeproduct/area;
  }
  return(h);

}


double Mesh::rhoTri(size_t iTri)
{
  //! returns the radius of the inscribed circle
  double rho=0.0;
  if(consistentState && nTri())
  {
    double area=AreaTri(iTri);
    Triangle & Tria=triangles[iTri];
    std::vector<Point> TriaPts(3);
    for(short int iPt=0; iPt<3; iPt++)
    {
      for(short int ic=0; ic<3; ic++)
      {
        TriaPts[iPt].coord[ic]=points[Tria.vertex[iPt]].coord[ic];
      }
    }
  
    short int * permutation = new short int[3];
    permutation[0]=1;
    permutation[1]=2;
    permutation[2]=0;
    double perimeter=0.0;
    for(short int iEdge=0; iEdge<3; iEdge++)
    {
      double dx=(TriaPts[iEdge].x()-TriaPts[permutation[iEdge]].x());
      double dy=(TriaPts[iEdge].y()-TriaPts[permutation[iEdge]].y());
      double dz=(TriaPts[iEdge].z()-TriaPts[permutation[iEdge]].z());
      perimeter=perimeter+sqrt(dx*dx+dy*dy+dz*dz);
    }
    delete [] permutation;
    rho=area/perimeter;
  }
  return(rho);
}


double Mesh::qTri(size_t iTri)
{
  double quality=0.0;
  if(consistentState && nTri())
  {
    quality=0.5*hTri(iTri)/rhoTri(iTri);  
  }
  
  return(quality);
}

double Mesh::AreaTri(size_t iTri) const
{
  //! returns the area of a triangle
  double area=0.0;
  if(consistentState && nTri())
  {
    const Triangle & Tria=triangles[iTri];
    std::vector<Point> TriaPts(3);
    for(short int iPt=0; iPt<3; iPt++)
    {
      for(short int ic=0; ic<3; ic++)
      {
        TriaPts[iPt].coord[ic]=points[Tria.vertex[iPt]].coord[ic];
      }
    }
    std::vector<double> v1(3,0), v2(3,0),vprod(3,0);
    for(short int ic=0; ic<3; ic++)
    {
      v1[ic]=TriaPts[1].coord[ic]-TriaPts[0].coord[ic];
      v2[ic]=TriaPts[2].coord[ic]-TriaPts[0].coord[ic];        
    }
    vprod[0]=v1[1]*v2[2]-v2[1]*v1[2];
    vprod[1]=v1[2]*v2[0]-v2[2]*v1[0];
    vprod[2]=v1[0]*v2[1]-v2[0]*v1[1];
  
    area=0.5*sqrt(vprod[0]*vprod[0]+vprod[1]*vprod[1]+vprod[2]*vprod[2]);
  }
  return(area);
}

double Mesh::hTet(size_t iTet)
{

  double h=0.0;
  if(consistentState && nTet())
  {
    double volume=VolTet(iTet);
    Tetrahedron Tet=tetrahedra[iTet];
  }
  return(h);
}


double Mesh::rhoTet(size_t iTet)
{
  double rho=0.0;
  if(consistentState && nTet())
  {
    double volume=VolTet(iTet);
    double area=AreaTet(iTet);
    rho=3.0*volume/area;
  }
  return(rho);
}


double Mesh::AreaTet(size_t iTet)
{
  double area=0.0;
  if(consistentState && nTet())
  {
    Tetrahedron Tet=tetrahedra[iTet];
    std::vector<Point> TetPts(4);
    for(short int iPt=0; iPt<4; iPt++)
    {
      for(short int ic=0; ic<3; ic++)
      {
        TetPts[iPt].coord[ic]=points[Tet.vertex[iPt]].coord[ic];
      }
    }
    std::vector<double> v1(3,0), v2(3,0),v3(3,0);
    std::vector<double> v4(3,0), v5(3,0);
    std::vector<double> vprod1(3,0), vprod2(3,0),vprod3(3,0),vprod4(3,0);
    
    for(short int ic=0; ic<3; ic++)
    {
      v1[ic]=TetPts[1].coord[ic]-TetPts[0].coord[ic];
      v2[ic]=TetPts[2].coord[ic]-TetPts[0].coord[ic];
      v3[ic]=TetPts[3].coord[ic]-TetPts[0].coord[ic];

      v4[ic]=TetPts[2].coord[ic]-TetPts[1].coord[ic];
      v5[ic]=TetPts[3].coord[ic]-TetPts[1].coord[ic];
    }

    vprod1[0]=v1[1]*v2[2]-v2[1]*v1[2];
    vprod1[1]=v1[2]*v2[0]-v2[2]*v1[0];
    vprod1[2]=v1[0]*v2[1]-v2[0]*v1[1];
    double area1=0.5*sqrt(vprod1[0]*vprod1[0]+vprod1[1]*vprod1[1]+vprod1[2]*vprod1[2]);
  
    vprod2[0]=v1[1]*v3[2]-v3[1]*v1[2];
    vprod2[1]=v1[2]*v3[0]-v3[2]*v1[0];
    vprod2[2]=v1[0]*v3[1]-v3[0]*v1[1];
    double area2=0.5*sqrt(vprod2[0]*vprod2[0]+vprod2[1]*vprod2[1]+vprod2[2]*vprod2[2]);

    vprod3[0]=v2[1]*v3[2]-v3[1]*v2[2];
    vprod3[1]=v2[2]*v3[0]-v3[2]*v2[0];
    vprod3[2]=v2[0]*v3[1]-v3[0]*v2[1];
    double area3=0.5*sqrt(vprod3[0]*vprod3[0]+vprod3[1]*vprod3[1]+vprod3[2]*vprod3[2]);

    vprod4[0]=v4[1]*v5[2]-v5[1]*v4[2];
    vprod4[1]=v4[2]*v5[0]-v5[2]*v4[0];
    vprod4[2]=v4[0]*v5[1]-v5[0]*v4[1];
    double area4=0.5*sqrt(vprod4[0]*vprod4[0]+vprod4[1]*vprod4[1]+vprod4[2]*vprod4[2]);

    area=area1+area2+area3+area4;
  }
  return(area);
}


double Mesh::VolTet(size_t iTet) const
{
  //! returns the volume of a tetra
  double volume=0.0;
  if(consistentState && nTet())
  {
    const Tetrahedron & Tet=tetrahedra[iTet];
    std::vector<Point> TetPts(4);
    for(short int iPt=0; iPt<4; iPt++)
    {
      for(short int ic=0; ic<3; ic++)
      {
        TetPts[iPt].coord[ic]=points[Tet.vertex[iPt]].coord[ic];
      }
    }
    std::vector<double> v1(3,0), v2(3,0),v3(3,0),vprod(3,0);
    for(short int ic=0; ic<3; ic++)
    {
      v1[ic]=TetPts[1].coord[ic]-TetPts[0].coord[ic];
      v2[ic]=TetPts[2].coord[ic]-TetPts[0].coord[ic];
      v3[ic]=TetPts[3].coord[ic]-TetPts[0].coord[ic];                
    }
    //2*face area normal
    vprod[0] =       v1[1]*v2[2] - v1[2]*v2[1];
    vprod[1] = -1.0*(v1[0]*v2[2] - v1[2]*v2[0]);
    vprod[2] =       v1[0]*v2[1] - v1[1]*v2[0];
    
    volume=(vprod[0]*v3[0]+vprod[1]*v3[1]+vprod[2]*v3[2])/6.0;
    volume=sqrt(volume*volume);
  }
  return(volume);
}

void Mesh::evalBoundaryLabels(bool debug,std::string debugDir,size_t print_interval)
{
  /*
    This routine evaluates the boundary labels
    - regionLabels is the set with the list of all the boundary label values
    - pointRegions is a map (region subdivision type) that, for each region label, 
        maps the point IDs with that label. Two regions can share points a priori
  */
  
  if(consistentState && nTri())
  {
    nbElToRegionLab.clear();
    pointRegions.clear();
    regionLabels.clear();
    typedef std::map<size_t, regionSetType > nodeToRegionMapType;
    typedef std::map<size_t, connectSetType> connectivityType;
        
    connectivityType connectivity; //connectivity of triangles (surface connectivity)
    
    //NodeToregionMap is a map between a  point and its regions; a point at the beginning can have more than one region
    nodeToRegionMapType NodeToregionMap;
    
    // determine initial region labels; 
    // determine initial region to node mapping
    // fill connectivity of surface points (on triangles) and node to region mapping
    // INITIALIZATION
    for(size_t iTri=0; iTri<nTri(); iTri++)
    {
      int labOfRegion=triangles[iTri].regionLabel;
      regionLabels.insert(labOfRegion);
      for(short int iv=0; iv<3; iv++)
      {
          size_t p0=triangles[iTri].vertex[iv];
          (pointRegions[labOfRegion]).insert(p0);
          (NodeToregionMap[p0]).insert(labOfRegion);
          for(short int jv=iv+1; jv<3; jv++)
          {
            size_t p1=triangles[iTri].vertex[jv];
            connectivity[p0].insert(p1);
            connectivity[p1].insert(p0);
          }
      }
    }//end loop on triangles
    
    std::set<int> regionsDebug;
    if(debug)
    {
        regionSubdivisionTypeIterator itRegToNodeMap=pointRegions.begin();
        int ireg=itRegToNodeMap->first;
        size_t nMax=(itRegToNodeMap->second).size();
        
        
        for(itRegToNodeMap=pointRegions.begin(); itRegToNodeMap!=pointRegions.end(); ++itRegToNodeMap)
        {
            size_t nPtReg=(itRegToNodeMap->second).size();
            if(nPtReg>nMax)
            {
                ireg=itRegToNodeMap->first;    
            }
        }
        
        regionsDebug.insert(ireg);
    }
        
    //algo on divisions starts; iterate on the region labels
    for(regionSubdivisionTypeIterator itRegToNodeMap=pointRegions.begin(); itRegToNodeMap!=pointRegions.end(); ++itRegToNodeMap)
    {
      int labelOfRegions=itRegToNodeMap->first;
      // iterate on points belonging to the region labelOfRegions
      // extract the local connectivity of the region regionconnect
      connectivityType regionconnect;
      regionconnect.clear();
      for(connectSetTypeIterator cRegIt=(itRegToNodeMap->second).begin(); cRegIt!=(itRegToNodeMap->second).end(); ++cRegIt)
      {
        // node belonging to region
        size_t node=*cRegIt;
        // explore connectivity of node; extract only connections belonging to the same region
        connectSetType connections;

        for(connectSetTypeIterator connectiter=connectivity[node].begin(); connectiter!=connectivity[node].end(); ++connectiter)
        {
          //NodeToregionMap[*connectiter]: set<int> of region labels
          // connectiter belongs to the region under study
          if((NodeToregionMap[*connectiter].find(labelOfRegions))!=(NodeToregionMap[*connectiter].end()))
          {
            connections.insert(*connectiter);
          }
        } //end loop with connectiter
        regionconnect.insert(std::pair<size_t, connectSetType> (node,connections) );
      }// end loop with cRegIt iterator

      //first: determine a seed point
      size_t seed= *((itRegToNodeMap->second).begin());
      for(connectSetTypeIterator seedIter=(itRegToNodeMap->second).begin(); seedIter!=(itRegToNodeMap->second).end(); ++seedIter)
      {
        seed = *seedIter;
        // full inside the region (not on a boundary of two)
        if(NodeToregionMap[seed].size()==1)
        {
          break;
        }
      }
      std::set<size_t> RegionNodes, Queue;
      Queue.insert(seed);
      RegionNodes.insert(seed);
      /* Here the algorithm grows: expands the region RegionNodes using the regional connectivity
         it iterates until all the nodes connected are covered
      */
      
      size_t localCounter=0;
      while(!Queue.empty())
      {
        size_t candidate= *(Queue.begin());
        Queue.erase(candidate);
        //regionconnect : map<size_t, connectSetType> isw the local connectivity on the region
        for(connectSetTypeIterator countNode=regionconnect.at(candidate).begin(); countNode!=regionconnect.at(candidate).end(); ++countNode)
        {
          //check if countNode is inside the region
          if(RegionNodes.find(*countNode)==RegionNodes.end())
          {
            //if not member of nodesRegion, add it
            RegionNodes.insert(*countNode);
            Queue.insert(*countNode);
          }
        }
        //debug (printing of the output)
        if(debug)
        {
          if (regionsDebug.find(labelOfRegions) != regionsDebug.end()) //region inside the list of regions to post-process
          {
              bool printSol=false;
              size_t nPtLabeled=RegionNodes.size();
              size_t nQuotient= std::floor(nPtLabeled/print_interval);
              if(nQuotient>localCounter)
              {
                  localCounter=nQuotient;
                  printSol=true;
              }
              if(Queue.empty())
              {
                localCounter=localCounter+1;
                printSol=true;
              }
              if(printSol )
              { 
                  // label the points
                  std::vector<double> meshLabels(this->nPt(),NAN);
                  for(std::set<size_t>::iterator debugit = RegionNodes.begin(); debugit!=RegionNodes.end(); ++debugit)
                  {
                      meshLabels[*debugit]=labelOfRegions;
                  }
                  VtkWriter writerVTK(this,true);
                  writerVTK.setOutputDir(debugDir);
                  std::ostringstream indexIte, indexReg;
                  indexReg<<labelOfRegions;
                  indexIte<<nQuotient;
                  std::string fName="regions_debug_"+indexReg.str()+"_"+indexIte.str();
                  writerVTK.setPrefixName(fName);
                  // now open the files and write the mesh geometry
                  writerVTK.openFileForOutput();
                  std::string varName="debugVar_"+indexReg.str();
                  writerVTK.writeVariable(meshLabels, varName,VtkWriter::Scalar);
                  writerVTK.CloseFile();
              }
          }
        }
      }//end of while
      //RegionNodes: nodes of the connected region with label itRegToNodeMap->first
      //check if ... itRegToNodeMap->second : list of points belonging to (set)
      //if is the case, there are some point to move
      if((itRegToNodeMap->second).size() !=RegionNodes.size()  )
      {
        //first: create a new label region
        int newRegionLabel=1+(*(regionLabels.rbegin()));
        regionLabels.insert(newRegionLabel);
        if(debug)
        {
          if (regionsDebug.find(labelOfRegions) != regionsDebug.end()) //"mother region" is inside the list of regions to post-process
          {
              regionsDebug.insert(newRegionLabel);
          }
        }
        
        //copy inside newSet the whole set of point belonging to the current region 
        connectSetType newSet=(itRegToNodeMap->second);

        //Delete points identified to the current region
        for(connectSetTypeIterator reg_iter=RegionNodes.begin();reg_iter!=RegionNodes.end(); ++reg_iter)
        {
          newSet.erase(*reg_iter);
        }
        //insert the new region inside pointRegions
        pointRegions.insert(std::pair<int, connectSetType>(newRegionLabel,newSet) );
        
        //nodeToRegionMapType = <size_t, set<int> >
        //remove point not belonging to region
         for(connectSetTypeIterator reg_iter=newSet.begin();reg_iter!=newSet.end(); ++reg_iter)
         {
          (itRegToNodeMap->second).erase(*reg_iter);
          NodeToregionMap[*reg_iter].erase(labelOfRegions);
          NodeToregionMap[*reg_iter].insert(newRegionLabel);
         }
      }// end if on size
    }//algo on divisions ends
    //now re-labeling of triangles
    for(size_t iTri=0; iTri<nTri(); iTri++)
    {
      for(short int iv=0; iv<3; iv++)
      {
        size_t p0=triangles[iTri].vertex[iv];
        if(NodeToregionMap[p0].size()==1)
        {
          triangles[iTri].regionLabel=*(NodeToregionMap[p0].begin());
          break;
        }
      }
    }
    

    for(regionSubdivisionTypeIterator iter=pointRegions.begin(); iter!=pointRegions.end(); ++iter)
    {
      int labreg=iter->first;
      size_t sizereg=(iter->second).size();
      nbElToRegionLab.insert(std::pair<size_t,int>(sizereg,labreg));
    }
    //now endo and epi of type long int for Laplace solver
    std::multimap<size_t,int>::const_reverse_iterator it=++(nbElToRegionLab.rbegin());
    _endoLabel =it->second;
    for(connectSetTypeIterator ite=pointRegions.at(_endoLabel).begin();ite!=pointRegions.at(_endoLabel).end(); ++ite)
    {
      _Endo.insert(static_cast<long int>(*ite));
    }

    it=(nbElToRegionLab.rbegin());
    _epiLabel =it->second;
    for(connectSetTypeIterator ite=pointRegions.at(_epiLabel).begin();ite!=pointRegions.at(_epiLabel).end(); ++ite)
    {
      _Epi.insert(static_cast<long int>(*ite));
    }
    
    //now check that endo and epi are correct
    
    Point endoBar, epiBar;
    for(std::set<long int>::iterator endo_iter=_Endo.begin(); endo_iter !=_Endo.end(); ++endo_iter)
    {
      for(short int  jcoord=0; jcoord<3; jcoord++)
      {
        endoBar.coord[jcoord] = endoBar.coord[jcoord] + (points[*endo_iter]).coord[jcoord];
      }
    }
    
    for(std::set<long int>::iterator epi_iter=_Epi.begin(); epi_iter !=_Epi.end(); ++epi_iter)
    {
      for(short int  jcoord=0; jcoord<3; jcoord++)
      {
        epiBar.coord[jcoord] = epiBar.coord[jcoord] + (points[*epi_iter]).coord[jcoord];
      }
    }
    
    for(short int  jcoord=0; jcoord<3; jcoord++)
    {
      endoBar.coord[jcoord] = endoBar.coord[jcoord]/_Endo.size();
      epiBar.coord[jcoord] =  epiBar.coord[jcoord]/_Epi.size();
    }

    //now the distances
    std::vector<double>endodist(2,1.e32);
    std::vector<double>epidist(2,1.e32);
    
    for(std::set<long int>::iterator endo_iter=_Endo.begin(); endo_iter !=_Endo.end(); ++endo_iter)
    {
      std::vector<double> dist(2,0.0);
      for(short int  jcoord=0; jcoord<3; jcoord++)
      {
        dist[0] = dist[0] + (points[*endo_iter].coord[jcoord]-endoBar.coord[jcoord])*(points[*endo_iter].coord[jcoord]-endoBar.coord[jcoord]);
        dist[1] = dist[1] + (points[*endo_iter].coord[jcoord]-epiBar.coord[jcoord])*(points[*endo_iter].coord[jcoord]-epiBar.coord[jcoord]);
      }
      dist[0] = sqrt(dist[0]);
      dist[1] = sqrt(dist[1]);
      for(short int jdist=0; jdist<2; jdist++)
      {
        if(dist[jdist]<endodist[jdist])
        {
          endodist[jdist] = dist[jdist];
        }
      }
    }
    
    for(std::set<long int>::iterator epi_iter=_Epi.begin(); epi_iter !=_Epi.end(); ++epi_iter)
    {
      std::vector<double> dist(2,0.0);
      for(short int  jcoord=0; jcoord<3; jcoord++)
      {
        dist[0] = dist[0] + (points[*epi_iter].coord[jcoord]-endoBar.coord[jcoord])*(points[*epi_iter].coord[jcoord]-endoBar.coord[jcoord]);
        dist[1] = dist[1] + (points[*epi_iter].coord[jcoord]-epiBar.coord[jcoord])*(points[*epi_iter].coord[jcoord]-epiBar.coord[jcoord]);
      }
      dist[0] = sqrt(dist[0]);
      dist[1] = sqrt(dist[1]);
      for(short int jdist=0; jdist<2; jdist++)
      {
        if(dist[jdist]<endodist[jdist])
        {
          epidist[jdist] = dist[jdist];
        }
      }
    }

    //now: endodist: distance of endocardium point form endo center and epicenter
    //now: epidist: distance of endocardium point form endo center and epicenter
    std::vector<bool> cmpVec(2,false);
    cmpVec[0]=(endodist[0]<=epidist[0]);
    cmpVec[1]=(endodist[1]<=epidist[1]);
    if(!(cmpVec[0] | cmpVec[1]))
    {
      std::set<long int> _endotmp(_Epi);
      _Epi.clear();
      _Epi=_Endo;
      _Endo.clear();
      _Endo=_endotmp;
       _endotmp.clear();
    }

  }// end if on consistence of mesh
}

//clear all except endo and epi; used to free memory
void Mesh::unsetBoundaryLabels()
{
  pointRegions.clear();
  regionLabels.clear();
  nbElToRegionLab.clear();
  _searchTree.clear();
  info.emptyMem();
}

void Mesh::unsetEndoEpiSets()
{
  _Endo.clear();
  _Epi.clear();
}

void Mesh::writeBoundaryLabels(std::string & fileDir, std::string & FileName)
{
    for(regionSubdivisionTypeIterator regIter=pointRegions.begin();regIter!=pointRegions.end();++regIter)
    {
        std::ofstream surfLabFile;
        std::string SurfLabelFilename;
        if(((regIter->first)!=_epiLabel) &&  ((regIter->first)!=_endoLabel)   )
        {
          std::ostringstream regionLabel;
          regionLabel << regIter->first;
          SurfLabelFilename=fileDir+ "/"+FileName+"_"+regionLabel.str()+".vtx";
        }
        else
        {
          if((regIter->first)==_epiLabel)
          {
            SurfLabelFilename=fileDir+ "/"+FileName+"_"+"epi"+".vtx";
          }
          else
          {
            SurfLabelFilename=fileDir+ "/"+FileName+"_"+"endo"+".vtx";
          }
        
        }
        
        surfLabFile.open(SurfLabelFilename.c_str());
        surfLabFile<<(regIter->second).size()<<std::endl;
        for(connectSetTypeIterator pointSetIter=(regIter->second).begin(); pointSetIter !=(regIter->second).end(); ++pointSetIter )
        {
          surfLabFile<<*pointSetIter<<std::endl;
        }
        surfLabFile.close();
    }
}

void Mesh::initializePtVector(size_t dim)
{
  points.clear();
  points.resize(dim);
  if(consistentState==false && tetrahedra.size())
  {
    consistentState=true;
  }
}


void Mesh::initializeTetraVector(size_t dim)
{
  tetrahedra.clear();
  tetrahedra.resize(dim);
  if(consistentState==false && points.size())
  {
    consistentState=true;
  }
}

void Mesh::initializeConnectivities()
{
  if(consistentState==true)
  {

    _faceToFace.resize(this->nTet());
    _elemBoundaryTris.resize(this->nTet());
    pointToElemConnectionType _conn_nodes; // not used
    std::cout << " Producing node-element connectivity array..." <<std::endl;
    _conn_nodes.resize(this->nPt());
    for(size_t iTetra=0; iTetra<nTet(); iTetra++)
    {
      const Tetrahedron & Tetra = tetrahedra[iTetra];
      for(short int iVertex=0; iVertex<4; iVertex++)
      {
        _conn_nodes[Tetra.vertex[iVertex]].insert(iTetra);
      }
    }
    
    pointToElemConnectionType _conn_nodesTris; //not used
    std::cout << " Producing node-triangle connectivity array..." <<std::endl;
    _conn_nodesTris.resize(this->nPt());
    for(size_t iTria=0; iTria<nTri(); iTria++)
    {
      const Triangle & Tria = triangles[iTria];
      for(short int iVertex=0; iVertex<3; iVertex++)
      {
        _conn_nodesTris[Tria.vertex[iVertex]].insert(iTria);
      }
      const Tetrahedron & Tetra = tetrahedra[triaToTet[iTria]];
      std::vector<bool> checkSide(4,false);
      short int matchingCounter=0;
      for(short int iTetV=0; iTetV<4; iTetV++)
      {
        size_t iPt=Tetra.vertex[iTetV];
        checkSide[iTetV]=false;
        for(short int iTriV=0; iTriV < 3; iTriV++)
        {
          if(Tria.vertex[iTriV]==iPt)
          {
            checkSide[iTetV]=true;
            matchingCounter=matchingCounter+1;
            break;
          }
        }
      }
      //<short int,size_t>
      
      if(matchingCounter==3)
      {
        for(short int iTetV=0; iTetV<4; iTetV++)
        {
          if(!checkSide[iTetV])
          {
            _elemBoundaryTris[triaToTet[iTria]].insert(std::pair<short int, size_t >(iTetV,iTria) );    
            break;
          }
        }
      }
    }
    
    std::cout << "Constructing element face-to-face look-up table... "<<std::endl;
    for(size_t iTetra=0; iTetra<nTet(); iTetra++)
    {
      connectSetType elemConnElem;
      connectSetTypeIterator it_ce;
      elemConnElem.clear();
      const Tetrahedron & Tetra = tetrahedra[iTetra];
      // Iterates over each node in current element
		  for(short int iVertex=0; iVertex<4; iVertex++)
		  {
			  size_t node = Tetra.vertex[iVertex];
			  for(connectSetTypeIterator it =_conn_nodes[node].begin(); it != _conn_nodes[node].end(); ++it  )
			  {
          elemConnElem.insert(*it);
			  }
		  }
      elemConnElem.erase(iTetra);
      // Creates a set with all nodes defining element
		  connectSetType elemNodes;
		  elemNodes.clear();
		  for(short int iVertex=0; iVertex<4; iVertex++)
		  {
		    elemNodes.insert(Tetra.vertex[iVertex]);
		  }
			// Iterates over each triangle face
			for(short int iVertex=0; iVertex<4; iVertex++)
			{
			  // Erases one node (which then defines a triangle)
			  elemNodes.erase(Tetra.vertex[iVertex]);
			  
			  // Iterates-over all elements in list of connected elements
  			for(it_ce=elemConnElem.begin();it_ce!=elemConnElem.end();it_ce++)
	  		{
  				
  				size_t connElem = *it_ce;
  				
  				// Tries to find all 3 current triangle nodes (in turn) within this connected element
				  short int p = 0;
				  for(short int gVertex= 0 ;gVertex < 4; gVertex++)
				  { 
					    if(elemNodes.find(tetrahedra[connElem].vertex[gVertex]) != elemNodes.end())
					    {
					      p++;
					    }
				  }
          // If p = 3 the we've found this triangle and thus this element must connect via the current element face
				  if(p == 3)
				  {
					  _faceToFace[iTetra].insert(std::pair<short int, size_t > (iVertex, connElem) );
					  break;
				  }
	  		}
			  
			  // Puts the erased node back before we try to move-on to the next triangle
			  elemNodes.insert(Tetra.vertex[iVertex]);
      }//end iteration on each triangle face    
    }// end loop on itetra
  



    size_t badElement = 0;
	  for(size_t iTetra=0;iTetra<nTet();iTetra++)
	  {
	    size_t codim = (4-_faceToFace[iTetra].size());
	    _elemBoundaryTris[iTetra].size();
		  for(int h=0;h<4;h++)
		  { 
			  if( (codim>0)  && (codim!=_elemBoundaryTris[iTetra].size()))
			  	badElement++;
		  }
	  }
	  std::cout << "Number of bad elements = " << badElement <<std::endl;

    _endoTris.clear();
    _epiTris.clear();
    typedef std::set<long int>::iterator surfiterType;

    for(surfiterType endoIter=_Endo.begin(); endoIter!=_Endo.end(); ++endoIter)
    {
      for(connectSetTypeIterator itsubset=_conn_nodesTris[*endoIter].begin(); itsubset !=_conn_nodesTris[*endoIter].end(); ++itsubset)
      {
        _endoTris.insert(*itsubset);
      }
    }
    for(surfiterType epiIter=_Epi.begin(); epiIter!=_Epi.end(); ++epiIter)
    {
      for(connectSetTypeIterator itsubset=_conn_nodesTris[*epiIter].begin(); itsubset !=_conn_nodesTris[*epiIter].end(); ++itsubset)
      {
        _epiTris.insert(*itsubset);
      }
    }

  }//end if on consistent state
}

void Mesh::writeCarpMesh(std::string outputFileName, bool binary, double rescaling)
{
  writePoints(outputFileName,rescaling, binary);  
  writeElements(outputFileName, binary);
  writeFakeFiebers(outputFileName, binary);
}

void Mesh::writePoints(std::string outputFileName, double rescaling, bool binary)
{
  short int precision=12;
  short int width=11;
  typedef float binaryPtsReal;
  
  size_t nPt=points.size();
  std::ofstream ptFile;
  if(binary)
  {
    std::string pfileName=outputFileName+".bpts";
    ptFile.open(pfileName.c_str(),std::ios::out | std::ios::binary);
    int npInt= static_cast<int>(nPt);
    int is_big_endian=static_cast<int>(!isLittleEndian());
    float chksum_pts = static_cast<float>(rescaling)*(info.chksum_pts());
    char * header = new char[HEADER_SIZE];
    sprintf(header,"%d %d %g # npts is_big chksum\n",npInt,is_big_endian,chksum_pts);
    //ptFile.write((char*) &header,HEADER_SIZE*sizeof(char));
    ptFile<<header;
    delete [] header;
    header=NULL;
  }
  else
  {
    std::string pfileName=outputFileName+".pts";
    ptFile.open(pfileName.c_str(),std::ios::out );
    ptFile<<nPt<<std::endl;
  }
  
  for(size_t iPt=0; iPt<nPt; iPt++)
  {
    Point  & Pt=points[iPt];
    if(binary)
    {
      binaryPtsReal x=static_cast<binaryPtsReal>(Pt.x()*rescaling);
      binaryPtsReal y=static_cast<binaryPtsReal>(Pt.y()*rescaling);
      binaryPtsReal z=static_cast<binaryPtsReal>(Pt.z()*rescaling);

      
      ptFile.write((char*) &x,sizeof(binaryPtsReal));
      ptFile.write((char*) &y,sizeof(binaryPtsReal));
      ptFile.write((char*) &z,sizeof(binaryPtsReal));
    }
    else
    {
      ptFile<<" "<<std::setw(width)<<std::setprecision(precision)<<rescaling*Pt.x()<<" "
                 <<std::setw(width)<<std::setprecision(precision)<<rescaling*Pt.y()<<" "
                 <<std::setw(width)<<std::setprecision(precision)<<rescaling*Pt.z()<<std::endl;
    }
  }
  ptFile.close();
}


void Mesh::writeElements(std::string outputFileName, bool binary)
{
  typedef enum elem_t {Tetra, Hexa, Octa, Pyramid, Prism, Quad, Tri, Line} Elem_t;
  Elem_t elemType;
  size_t nElem=0;
  bool is_3D=false;
  std::string elType;
  size_t * vertex=NULL;
  unsigned char nbV=0;
  
  if(this->nTet()>0)
  {
    nElem=this->nTet();
    is_3D=true;
    elType="Tt";
    elemType=Tetra;
    nbV=4;
  }
  else
  {
    nElem=this->nTri();
    elType="Tr";
    elemType=Tri;
    is_3D=false;
    nbV=3;
  }
  vertex= new size_t[nbV];
  std::ofstream elFile;

  long int chksum_elems = 0;
  long int chksum_types = static_cast<long int>((elemType)*nElem);
  long int chksum_tags  = 0;
  long int chksum=chksum_elems+chksum_types+chksum_tags;

  if(binary)
  {
    std::string elfileName=outputFileName+".belem";
    elFile.open(elfileName.c_str(),std::ios::out | std::ios::binary);
    int nelInt= static_cast<int>(nElem);
    int is_big_endian=static_cast<int>(!isLittleEndian());
    char * header = new char[HEADER_SIZE];
    sprintf(header, "%d %d %ld # nelems is_big chksum\n",nelInt,is_big_endian,chksum);
    //elFile.write((char*) &header,HEADER_SIZE*sizeof(char));
    elFile<<header;
    delete [] header;
    header=NULL;
  }
  else
  {
    std::string elfileName=outputFileName+".elem";
    elFile.open(elfileName.c_str(),std::ios::out );
    elFile<<nElem<<std::endl;
  }
  
  
  
  for(size_t iElem=0; iElem<nElem; iElem++)
  {
    int regionLabel=0;
    for(unsigned char iV=0; iV<nbV; iV++ )
    {
      if(is_3D)
      {
        vertex[iV]=tetrahedra[iElem].vertex[iV];
        regionLabel=tetrahedra[iElem].regionLabel;
      }
      else
      {
        vertex[iV]=triangles[iElem].vertex[iV];
        regionLabel=triangles[iElem].regionLabel;
      }
    }
    
  
    if(binary) 
    {
      
      chksum_tags=chksum_tags+regionLabel;
      elFile.write((char*) &elemType, sizeof(int));
      for(unsigned char iV=0; iV<nbV; iV++ )
      {
        chksum_elems=chksum_elems+vertex[iV];
        int iVertex=static_cast<int>(vertex[iV]);
        elFile.write((char*) &iVertex,sizeof(int));
      }
      elFile.write((char*) &regionLabel,sizeof(int));
      
    }
    else
    {
      elFile<<elType;
      for(unsigned char iV=0; iV<nbV; iV++)
      {
        elFile<<" "<<vertex[iV];
      }
      elFile<<" "<<regionLabel<<std::endl;
    }
  }
  delete [] vertex;
  vertex=NULL;
  if(binary) 
  {
    chksum=chksum_elems+chksum_types+chksum_tags;
    int nelInt= static_cast<int>(nElem);
    int is_big_endian=static_cast<int>(!isLittleEndian());
    char * header = new char[HEADER_SIZE];
    sprintf(header, "%d %d %ld # nelems is_big chksum\n",nelInt,is_big_endian,chksum);
    elFile.seekp (0);
    //elFile.write((char*) &header,HEADER_SIZE*sizeof(char));
    elFile<<header;
    delete [] header;
    header=NULL;

  }
  elFile.close();
}



void Mesh::writeFakeFiebers(std::string outputFileName, bool binary)
{
  short int precision=12;
  short int width=11;
  size_t nElem=0;
  
  if(this->nTet()>0)
  {
    nElem=this->nTet();
  }
  else
  {
    nElem=this->nTri();
  }

  std::ofstream fFibFile;

  if(binary)
  {
    std::string elfileName=outputFileName+".blon";
    fFibFile.open(elfileName.c_str(),std::ios::out | std::ios::binary);
    fFibFile<<"1"<<std::endl;
  }
  else
  {
    std::string elfileName=outputFileName+".lon";
    fFibFile.open(elfileName.c_str(),std::ios::out );
    fFibFile<<"1"<<std::endl;
  }
  
  for(size_t iElem=0; iElem<nElem; iElem++)
  {
    
    std::vector<double> fakefib(3,0);
    fakefib[0]=1.0;
    for(short int icomp=0; icomp<2; icomp++)
    {
      double component=fakefib[icomp];
      if(binary) 
      {
        fFibFile.write((char*) &component, sizeof(double));
      }
      else
      {
        fFibFile<<std::setw(width)<<std::setprecision(precision)<<component<<" ";
      }
    }
    double component=fakefib[2];
      if(binary) 
      {
        fFibFile.write((char*) &component, sizeof(double));
      }
      else
      {
        fFibFile<<std::setw(width)<<std::setprecision(precision)<<component<<std::endl;
      }
  }
  fFibFile.close();
}

void Mesh::writeVTKMesh(std::string outputFileName, bool binary, double rescaling)
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
  
  std::string fileName=outputFileName+".vtk";

  short int precision=12;
  

  size_t nPt=this->nPt();
  size_t nElem=0;
  bool is_3D=false;
  vtkIntType * vertex=NULL;
  short int nbV=0;
  vtkIntType typeCell=0;
  vtkIntType nbVoutput=0;
  if(this->nTet()>0)
  {
    nElem=this->nTet();
    is_3D=true;
    nbV=4;
    typeCell=10;
  }
  else
  {
    nElem=this->nTri();
    is_3D=false;
    nbV=3;
    typeCell=5;
  }
  nbVoutput=static_cast<vtkIntType>(nbV);


  if(binary && littleEndianMachine )
  {
    SwapBytes(&typeCell, sizeof(typeCell));
    SwapBytes(&nbVoutput, sizeof(nbVoutput));
  }


  
  vertex= new vtkIntType[nbV];
  std::ofstream VTKFile;
  
  if(binary)
  {
    VTKFile.open(fileName.c_str(),std::ios::out | std::ios::binary);
  }
  else
  {
    VTKFile.open(fileName.c_str(),std::ios::out );
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
    Point  & Pt=points[iPt];
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
        vertex[iV]=static_cast<vtkIntType>(tetrahedra[iElem].vertex[iV]);
      }
      else
      {
        vertex[iV]=static_cast<vtkIntType>(triangles[iElem].vertex[iV]);
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
      VTKFile<<nbVoutput;
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
  std::vector<int> plab(nPt,0);
  for(regionSubdivisionTypeIterator it=pointRegions.begin(); it!= pointRegions.end(); ++it)
  {
    for(connectSetTypeIterator itP=(it->second).begin(); itP!=(it->second).end(); ++itP)
    {
      plab[*itP]=(it->first);
    }
    
  }
  
  VTKFile<<"POINT_DATA "<<nPt<<std::endl;
  VTKFile<<"SCALARS NodeLabels FLOAT"<<std::endl;
  VTKFile<<"LOOKUP_TABLE default"<<std::endl;
  for(size_t iPt=0; iPt<nPt; iPt++)
  {
    vtkFloatType labelOfPt=static_cast<vtkFloatType>(plab[iPt]);
    
    if(binary)
    {
      if(littleEndianMachine)
      {
        SwapBytes(&labelOfPt, sizeof(labelOfPt));
      }
      VTKFile.write((char*) &labelOfPt,sizeof(vtkFloatType));
    }
    else
    {
      VTKFile<<labelOfPt;
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

  
  VTKFile.close();
}


std::vector<double> Mesh::copyLabelVectorForVTKOutput()
{
  std::vector<double> plab(nPt(),0);
  for(regionSubdivisionTypeIterator it=pointRegions.begin(); it!= pointRegions.end(); ++it)
  {
    for(connectSetTypeIterator itP=(it->second).begin(); itP!=(it->second).end(); ++itP)
    {
      plab[*itP]=static_cast<double>(it->first);
    }
  }
  
  return(plab);

}


void Mesh::clear()
{
  info.clear();
  points.clear();
  triangles.clear();
  tetrahedra.clear();
  triaToTet.clear();
  pointRegions.clear();
  regionLabels.clear();
  nbElToRegionLab.clear();
  _faceToFace.clear();
  _elemBoundaryTris.clear();
  _endoTris.clear();
  _epiTris.clear();
  _searchTree.clear();

  consistentState=false;
}


Mesh::~Mesh()
{
  clear();
}


void Mesh::extractBoundary(bool build_search_tree)
{
  /* This routine extracts the boundary faces and fill a bound_faces
     type  that is a multi-map; for each face the entries are:
     key:    among the faces nodes, is the smaller node global label
     value:  a pair, composed by:
             - the label of the Tetrahedra element the face belongs to
             - the list of the Global label of the nodes of the face
    This routine iterates on all the tetra, then:
    - extract the faces of the tetra;
    -- if the face of the tetra is not in the boundary faces list, the face is added
    -- if the face of the tetra is in the boundary faces list, the face is removed
    At the end of the execution, the remaining faces are those on the boundary
  */
  
  typedef std::pair <mapfacetype::iterator, mapfacetype::iterator> range_iter_type;
  if(consistentState)
  { 
    size_t nbTri=0;
    size_t nTet=tetrahedra.size();
    std::vector<short int * > permutation;
    permutation.resize(4);
    for(short int jface=0; jface<4; jface++)
    {
      permutation[jface] = new short int [3];
      for(short int jf=0; jf<3; jf++)
      {
        short int asum=jface+jf;
        short int rem= asum%4;
        permutation[jface][jf]=rem;
      }
    }
    mapfacetype bound_faces;
    // a multi-list of pairs:
    // the first argument is a search key corresponding to the point with smaller global label
    // the second argument is compose by a pair: element lable and node global labels
    bound_faces.clear();
    for(size_t iTet=0; iTet<nTet; iTet++)
    {
      std::vector<facetype> ordered_faces;    
      ordered_faces.resize(4);
      //ordered_faces: the 4 faces with nodes ordered
      for(short int jface=0; jface<4; jface++)
      {
        for(short int jf=0; jf<3; jf++)
        {
          ordered_faces[jface].insert(tetrahedra[iTet].vertex[(permutation[jface][jf])] );
        }
      }
      // now I use as key the value of the first node
      // if the face is inside bound_faces, I delete it
      // if the face is not inside bound_faces I add it
      for(short int jface=0; jface<4; jface++)
      {
        size_t key=*(ordered_faces[jface].begin());
        bool value_to_insert=true;
        if((bound_faces.count(key))>0)
        {
            // range of elements having the same key
            range_iter_type itmapRange=bound_faces.equal_range(key);
            // scroll on faces having the same key 
            for(mapfacetype::iterator itmap=itmapRange.first; itmap!=itmapRange.second; ++itmap)
            {
              if(ordered_faces[jface]==(itmap->second).second)
              {
                value_to_insert=false;
                bound_faces.erase(itmap);
                break;                
              }
            }
        }
        
        if(value_to_insert)
        {
          facetype face=ordered_faces[jface];
          faceLabtype labeledFace(iTet,face);
          bound_faces.insert(std::pair<size_t, faceLabtype > (key,labeledFace) );
        }
      }//end loop on jface

    }//end loop on Tetra

    for(short int jface=0; jface<4; jface++)
    {
      delete[] permutation[jface];
      permutation[jface] = NULL;
    }
    permutation.clear();
    evalTriangles(bound_faces, nbTri,build_search_tree, true);
  }//end if on consistent state

}

Mesh::IDsetType Mesh::extractTrianglesFromBBOX(const std::vector<std::vector<double> > & bb)
{
  return(extractSIDFromBBOX(bb));
}

//used to evaluate endianess
bool Mesh::isLittleEndian()
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
void Mesh::SwapBytes(void *pv, size_t n)
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


void Mesh::checkConnectivity()
{
  if(consistentState)
  {
    //first extract all the vertex labels within tetrahedra and triangles
    std::set<size_t> connected_vertices;
    for(size_t iTet=0; iTet<nTet(); iTet++)
    {
      for(short int iv=0; iv<4; iv++)
      {
        connected_vertices.insert(tetrahedra[iTet].vertex[iv]);
      }
    }
  
    for(size_t iTri=0; iTri<nTri(); iTri++)
    {
      for(short int iv=0; iv<3; iv++)
      {
        connected_vertices.insert(triangles[iTri].vertex[iv]);
      }
    }
    //if there are not connected points
    if(nPt()!=connected_vertices.size() )
    {
      size_t count=0;
      std::set<size_t>::iterator it;
      std::map<size_t,size_t> renumbering;
      std::vector<bool> is_connected(nPt(),false);
      //I creeate a mapping between the connected nodes and the new reordering
      for(it=connected_vertices.begin();it!=connected_vertices.end(); ++it)
      {
        renumbering.insert(std::pair<size_t,size_t>(*it,count));
        count++;
      }
      //First I renumber the Tetra Vertices
      for(size_t iTet=0; iTet<nTet(); iTet++)
      {
        for(short int iv=0; iv<4; iv++)
        {
          size_t old_index=tetrahedra[iTet].vertex[iv];
          is_connected[old_index]=true;
#ifndef NDEBUG
          size_t new_index=(renumbering.at(old_index));          
#else
          size_t new_index=(renumbering[old_index]);          
#endif
          tetrahedra[iTet].vertex[iv]=new_index;
        }
      }
      //After I renumber Tria Vertex
      for(size_t iTri=0; iTri<nTri(); iTri++)
      {
        for(short int iv=0; iv<3; iv++)
        {
          size_t old_index=triangles[iTri].vertex[iv];
          is_connected[old_index]=true;
#ifndef NDEBUG
          size_t new_index=(renumbering.at(old_index));
#else
          size_t new_index=(renumbering[old_index]);
#endif
          triangles[iTri].vertex[iv]=new_index;
        }
      }
      //extract not connected indices
      std::set<size_t> not_connected;
      for(size_t iPt=0; iPt<is_connected.size(); iPt++)
      {
        if(!is_connected[iPt])
        {
          not_connected.insert(iPt);
        }
      }
      //now i can delete unconnected vertices, starting by the end
      std::set<size_t>::reverse_iterator rit;
      for(rit=not_connected.rbegin(); rit!=not_connected.rend(); ++rit)
      {
        size_t index_to_remove=*rit;
        points.erase(points.begin()+index_to_remove);
      }
    }//end if not connected
  }
}



void Mesh::preprocessingOperations()
{
  this->checkConnectivity();
}



std::vector<double> Mesh::TetJacobian(size_t iTet) const
{
  std::vector<double> Jacobian(9,0.0);
  if(consistentState)
  {
    const Tetrahedron & Tetra=tetrahedra[iTet];
    for(short int jc=0; jc<3; jc++)
    {
      for(short int jv=0; jv<3; jv++)
      {
        Jacobian[3*(jc)+jv]=points[Tetra.vertex[1+jv]].coord[jc]-points[Tetra.vertex[0]].coord[jc];
      }
    }

  }
  return(Jacobian);
}



std::vector<double> Mesh::TetJacobianTransponse(size_t iTet) const
{
  std::vector<double> JacobianT(9,0.0);
  if(consistentState)
  {
    const Tetrahedron & Tetra=tetrahedra[iTet];
    for(short int jc=0; jc<3; jc++)
    {
      for(short int jv=0; jv<3; jv++)
      {
        JacobianT[3*(jv)+jc]=points[Tetra.vertex[1+jv]].coord[jc]-points[Tetra.vertex[0]].coord[jc];
      }
    }
  }
 
  return(JacobianT);
}

std::vector<double>Mesh::InvertA3X3(const std::vector<double> & Mat0) const
{
  //matrix is considered as row major
  
  std::vector<double> inverted(9,0);
  
  //rowMajor3X3Index(short int irow, short int jcol)
  
  double det =    Mat0[RM3X3Ind(0,0)]*( Mat0[RM3X3Ind(1,1)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(1,2)]*Mat0[RM3X3Ind(2,1)] ) +
             -1.0*Mat0[RM3X3Ind(0,1)]*( Mat0[RM3X3Ind(1,0)]*Mat0[RM3X3Ind(2,2)] - Mat0[RM3X3Ind(1,2)]*Mat0[RM3X3Ind(2,0)] ) +
                  Mat0[RM3X3Ind(0,2)]*( Mat0[RM3X3Ind(1,0)]*Mat0[RM3X3Ind(2,1)] - Mat0[RM3X3Ind(2,0)]*Mat0[RM3X3Ind(1,1)] ) ;

  
  
  
  /*for(short int ir=0; ir<3; ir++)
  {
    for(short int jc=0; jc<3; jc++)
    {
      double entry = (Mat0[RM3X3Ind((ir+1)%3,(jc+1)%3)]*Mat0[RM3X3Ind((ir+2)%3,(jc+2)%3)]) - (Mat0[RM3X3Ind((ir+1)%3,(jc+2)%3)]*Mat0[RM3X3Ind((ir+2)%3,(jc+1)%3)] );
      inverted[RM3X3Ind(jc,ir)]=entry;
    }
  }*/

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

short int Mesh::RM3X3Ind(short int irow, short int jcol) const
{
  short int index=3*irow + jcol;
  return(index);
}

size_t Mesh::TensorIJtoIndex(const size_t & I, const size_t & J,  const size_t & JDIM ) const
{
  size_t index= J+(JDIM*I);
  return(index);
}

size_t Mesh::TensorIJKtoIndex(const size_t & I, const size_t & J,  const size_t & K,const size_t & IDIM, const size_t & JDIM ) const
{
  size_t index= JDIM*(K*IDIM+I)+ J;
  return(index);
}

std::vector<double> Mesh::dimensionlessCoord(const std::vector<double> & coordVec) const
{
  std::vector<double> aCoord(3,0);
  for(unsigned char jc=0; jc<3; jc++)
  {
    aCoord[jc]=(coordVec[jc]-(info.bbox[jc])[0])/( (info.bbox[jc])[1]-(info.bbox[jc])[0] );
  }
  return(aCoord);
}

size_t Mesh::evalHashKey(const std::vector<double> & coordVec) const 
{
    size_t hashKey=0;
    std::vector<double> aCoord=dimensionlessCoord(coordVec);
    std::vector<size_t> IJK(3,0);
    for(unsigned char jc=0; jc<3; jc++)
    {
        IJK[jc]=std::floor((TREECSIZE-1)*aCoord[jc]);
    }
    hashKey=TensorIJKtoIndex(IJK[0], IJK[1], IJK[2],TREECSIZE, TREECSIZE);
    return(hashKey);
}


Mesh::IDsetType Mesh::extractSIDFromBBOX( const std::vector<std::vector<double> > & bb ) const
{
  IDsetType setID;
  std::vector<std::vector<size_t> > IJKbbox(3);
  std::vector<double> c0(3,0),c1(3,0);
  for(unsigned char jc=0; jc<3; jc++)
  {
    c0[jc]=(bb[jc])[0];
    c1[jc]=(bb[jc])[1];
    IJKbbox[jc].resize(2,0);
  }
  std::vector<double> c0adim=dimensionlessCoord(c0);
  std::vector<double> c1adim=dimensionlessCoord(c1);
  for(unsigned char jc=0; jc<3; jc++)
  {
    (IJKbbox[jc])[0]=std::floor((TREECSIZE-1)*c0adim[jc]);
    (IJKbbox[jc])[1]=std::floor((TREECSIZE-1)*c1adim[jc]);
  }
  for(size_t I=(IJKbbox[0])[0]; I<=(IJKbbox[0])[1]; I++)
  {
      for(size_t J=(IJKbbox[1])[0]; J<=(IJKbbox[1])[1]; J++)
      {
        for(size_t K=(IJKbbox[2])[0]; K<=(IJKbbox[2])[1]; K++)
        {
            size_t index=TensorIJKtoIndex(I, J, K,TREECSIZE, TREECSIZE);
            for(IDiteratorType it=_searchTree[index].begin();it!=_searchTree[index].end(); ++it)
            {
                setID.insert(*it);
            }
        }
      }
  }
  return(setID);

}

std::vector<double> Mesh::TetInvJacobian(size_t iTet) const
{
  std::vector<double> invJacobian(9,0);
  if(consistentState)
  {
    const std::vector<double> Jacobian=TetJacobian(iTet);
    invJacobian=InvertA3X3(Jacobian);
  }
  return(invJacobian);
}


std::vector<double> Mesh::TetInvJacobianTransponse(size_t iTet) const
{
  std::vector<double> invJacobianT(9,0);
  if(consistentState)
  {
    const std::vector<double> JacobianT=TetJacobianTransponse(iTet);
    invJacobianT=InvertA3X3(JacobianT);
  }
  return(invJacobianT);
}

std::vector<double> Mesh::TetraCentroid(size_t iTet) const
{
  std::vector<double> centroid(3,0);
  if(consistentState && this->nTet())
  {
    const Tetrahedron & Tetra= tetrahedra[iTet];
    for(unsigned char jv=0; jv<4; jv++)
    {
      for(unsigned char ic=0; ic<3; ic++)
      {
        centroid[ic]=centroid[ic]+0.25*points[Tetra.vertex[jv]].coord[ic];
      }
    
    }
  }
  return(centroid);
}


std::vector<double> Mesh::TriaCentroid(size_t iTri) const
{
  std::vector<double> centroid(3,0);
  if(consistentState && this->nTri())
  {
    const Triangle & Tria= triangles[iTri];
    for(unsigned char jv=0; jv<3; jv++)
    {
      for(unsigned char ic=0; ic<3; ic++)
      {
        centroid[ic]=centroid[ic]+points[Tria.vertex[jv]].coord[ic];
      }
    }
    for(unsigned char ic=0; ic<3; ic++)
    {
      centroid[ic]=centroid[ic]/3.0;
    }
  }
  return(centroid);
}

std::vector<double> Mesh::ElemCentroid(const Element & Elem) const
{
  std::vector<double> centroid(3,0);
  if(consistentState)
  {
      for(unsigned char jv=0; jv<Elem.nbV(); jv++)
      {
        for(unsigned char ic=0; ic<3; ic++)
        {
          centroid[ic]=centroid[ic]+points[Elem.vertex[jv]].coord[ic];
        }
      }

      for(unsigned char ic=0; ic<3; ic++)
      {
        centroid[ic]=centroid[ic]/static_cast<double>(Elem.nbV());
      }
  }

  return(centroid);
}


void Mesh::writeTetraCentroids(std::string outputFileName, double rescaling)
{
  short int precision=12;
  size_t nTet=this->nTet();
  if(consistentState && nTet)
  {
    std::string fname=outputFileName+".cpts";
    std::ofstream TetCentroidFile(fname.c_str());
    TetCentroidFile<<nTet<<std::endl;
    for(size_t iTet=0; iTet<nTet; iTet++)
    {
      std::vector<double> centroid=TetraCentroid(iTet);
      centroid[0]=centroid[0]*rescaling;
      centroid[1]=centroid[1]*rescaling;
      centroid[2]=centroid[2]*rescaling;
      TetCentroidFile<<std::setprecision(precision)<<centroid[0]<<" "
                     <<std::setprecision(precision)<<centroid[1]<<" "
                     <<std::setprecision(precision)<<centroid[2]<<std::endl;
    }  
    TetCentroidFile.close();
  }
}



void Mesh::writeTris(std::string outputFileName)
{
  size_t nTri=this->nTri();
  if(consistentState && nTri)
  {
      std::string fname=outputFileName+".tris";
      std::ofstream TrisFile(fname.c_str());
      TrisFile<<nTri<<std::endl;
      for(size_t iTri=0; iTri<this->nTri(); iTri++)
      {
        TrisFile<<triangles[iTri].vertex[0]<<" "<<triangles[iTri].vertex[1]<<" "<<triangles[iTri].vertex[2]<<std::endl;
      }
      TrisFile.close();
  }
}

void Mesh::meshRescaling(double rescaling)
{
  if(consistentState)
  {
    for(size_t iPt=0; iPt<this->nPt(); iPt++)
    {
      for(short int ic=0; ic<3; ic++)
      {
        points[iPt].coord[ic]=rescaling*points[iPt].coord[ic];
      }
    }
  }
}


