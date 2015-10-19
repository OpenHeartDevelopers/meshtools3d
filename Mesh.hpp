#ifndef _MESH_HPP
#define _MESH_HPP

#include <fstream>
#include <string>
#include<vector>
#include <string>
//#include <values.h>
#include <cfloat>
#include<set>
#include<map>

struct MeshInfo{
  std::vector<double> baricenter;
  std::vector<std::vector<double> > bbox;
  std::vector<float> checksum;
  MeshInfo()
  :baricenter(3,0.0),
  bbox(3),
  checksum(3,0.0)
  {
    bbox[0].resize(2);
    bbox[1].resize(2);
    bbox[2].resize(2);
    for(short int jdim=0; jdim<3; jdim++)
    {
        (bbox[jdim])[0]= DBL_MAX;
        (bbox[jdim])[1]=-DBL_MAX;
    }
  }
  void clear()
  {
    for(short int jdim=0; jdim<3; jdim++)
    {
        (bbox[jdim])[0]= DBL_MAX;
        (bbox[jdim])[1]=-DBL_MAX;
        baricenter[jdim]=0.0;
        checksum[jdim]=0.0;
    }
  }
  inline double & bx(){return baricenter[0];};
  inline double & by(){return baricenter[1];};
  inline double & bz(){return baricenter[2];};

  inline float chksum_x(){return checksum[0];};
  inline float chksum_y(){return checksum[1];};
  inline float chksum_z(){return checksum[2];};
  inline float chksum_pts(){return(checksum[0]+checksum[1]+checksum[2]); };
};

struct Point{
  std::vector<double> coord;
  Point()
  :coord(3,0)
  {}

  inline double & x(){return coord[0];};
  inline double & y(){return coord[1];};
  inline double & z(){return coord[2];};
  
};


class Element
{
  public:
    Element()
      :regionLabel(0),
      vertex(0),
      _vtkCellType(0)
      {};
  ~Element()
    {
      vertex.clear();
      regionLabel=0;
      _vtkCellType=0;
    };

  inline const short int nbV() const {return static_cast<short int>(vertex.size());};
  int regionLabel;
  std::vector<size_t> vertex;
  int _vtkCellType;
};


class Triangle
:public Element
{
  public:
    Triangle()
    :Element()
    {
        _vtkCellType=5;
        vertex.resize(3);
        vertex[0]=0;
        vertex[1]=0;
        vertex[2]=0;
    }
    inline size_t & p0(){return vertex[0];};
    inline size_t & p1(){return vertex[1];};
    inline size_t & p2(){return vertex[2];};
};


class Tetrahedron
:public Element
{
  public:
  Tetrahedron()
    :Element()
    {
        _vtkCellType=10;
        vertex.resize(4);
        vertex[0]=0;
        vertex[1]=0;
        vertex[2]=0;
        vertex[3]=0;
    }
    inline size_t & p0(){return vertex[0];};
    inline size_t & p1(){return vertex[1];};
    inline size_t & p2(){return vertex[2];};
    inline size_t & p3(){return vertex[3];};
};



class Mesh
{
  typedef std::set<int>  regionSetType;
  typedef std::set<size_t> connectSetType;
  typedef connectSetType::iterator connectSetTypeIterator;
  typedef std::map<int, connectSetType> regionSubdivisionType;
  typedef regionSubdivisionType::iterator regionSubdivisionTypeIterator;
  typedef std::set<size_t> facetype;
  typedef facetype::iterator facetype_iter;
  typedef std::pair<size_t,facetype> faceLabtype;
  typedef std::multimap<size_t, faceLabtype > mapfacetype;
  typedef std::vector<connectSetType> pointToElemConnectionType;

  // type for elemBoundaryTris  
  typedef std::map<short int, size_t> boundaryFaceInTetraType;
  typedef std::vector<boundaryFaceInTetraType> ElemToBoundaryFaceConnectionType;



  public:
  
    Mesh();
    ~Mesh();
    Mesh(const  std::string & inputFileame);
    Mesh(const  Mesh & srcMesh);
    void readFromFile(const  std::string & inputFileame);
    void initializePtVector(size_t dim);
    void initializeTetraVector(size_t dim);
    void preprocessingOperations();
    void extractBoundary();
    void evalBoundaryLabels();
    void unsetBoundaryLabels();
    void clear();
    void initializeConnectivities();
    void writeCarpMesh(std::string outputFileName, bool binary=false, double rescaling=1.0);
    void writeVTKMesh(std::string outputFileName, bool binary=false, double rescaling=1.0);
    void writeBoundaryLabels(std::string & fileDir, std::string & FileName);
    void writeTetraCentroids(std::string outputFileName, double rescaling=1.0);
    void writeTris(std::string outputFileName);
    void meshRescaling(double rescaling=1.0);
     
    inline size_t nPt() const {return points.size();};
    inline size_t nTri() const {return triangles.size();};
    inline size_t nTet() const {return tetrahedra.size();};
    inline MeshInfo & Info(){return info;};
    inline const MeshInfo & Info() const {return info;};
    inline const facetype & endoTria() const {return _endoTris;};
    inline const facetype & epiTria() const {return _epiTris;};
    inline const ElemToBoundaryFaceConnectionType & elemBoundaryTris() const {return _elemBoundaryTris;};
    inline const ElemToBoundaryFaceConnectionType & faceToFace() const {return _faceToFace;};
    
    inline const std::vector<Point>  & Pt() const {return points;};
    inline const std::vector<Triangle> & Tri() const {return triangles;};
    inline const std::vector<Tetrahedron> & Tet() const {return tetrahedra;};

    
    std::vector<double> copyLabelVectorForVTKOutput();
    
#ifndef NDEBUG
    inline Point & Pt(size_t iPt){return points.at(iPt);};
    inline Triangle & Tri(size_t iTri){return triangles.at(iTri);};
    inline Tetrahedron & Tet(size_t iTet){return tetrahedra.at(iTet);};
    inline size_t & triaToTetMap(size_t iTri){return triaToTet.at(iTri);};
    inline const Point & Pt(size_t iPt) const {return points.at(iPt);};
    inline const Triangle & Tri(size_t iTri) const {return triangles.at(iTri);};
    inline const Tetrahedron & Tet(size_t iTet) const {return tetrahedra.at(iTet);};
    inline const size_t & triaToTetMap(size_t iTri) const {return triaToTet.at(iTri);};
    inline const boundaryFaceInTetraType & elemBoundaryTris(size_t iTet) const {return _elemBoundaryTris.at(iTet);};
    inline const boundaryFaceInTetraType & faceToFace(size_t iTet) const {return _faceToFace.at(iTet);};
    
#else
    inline Point & Pt(size_t iPt){return points[iPt];};
    inline Triangle & Tri(size_t iTri){return triangles[iTri];};
    inline Tetrahedron & Tet(size_t iTet){return tetrahedra[iTet];};
    inline size_t & triaToTetMap(size_t iTri){return triaToTet[iTri];};
    inline const Point & Pt(size_t iPt) const {return points[iPt];};
    inline const Triangle & Tri(size_t iTri) const {return triangles[iTri];};
    inline const Tetrahedron & Tet(size_t iTet) const{return tetrahedra[iTet];};
    inline const size_t & triaToTetMap(size_t iTri) const{return triaToTet[iTri];};
    inline const boundaryFaceInTetraType & elemBoundaryTris(size_t iTet) const {return _elemBoundaryTris[iTet];};
    inline const boundaryFaceInTetraType & faceToFace(size_t iTet) const {return _faceToFace[iTet];};
#endif
    inline const std::set<long int> & Endocardium() const {  return(_Endo);};
    inline const std::set<long int> & Epicardium() const {  return(_Epi);};;  

    double hTri(size_t iTri); // diameter
    double rhoTri(size_t iTri);  //radius (inscribed)
    double qTri(size_t iTri);  //quality of triangle
    double AreaTri(size_t iTri) const; //area
    
    double hTet(size_t iTet); // diameter (circumscribed)
    double rhoTet(size_t iTet);  //radius (inscribed)
    double AreaTet(size_t iTet);
    double VolTet(size_t iTet) const; //volume
    std::vector<double> TetJacobian(size_t iTet) const;
    std::vector<double> TetJacobianTransponse(size_t iTet) const;    
    std::vector<double> TetInvJacobian(size_t iTet) const;
    std::vector<double> TetInvJacobianTransponse(size_t iTet) const;

    std::vector<double> TetraCentroid(size_t iTet) const;
    std::vector<double> TriaCentroid(size_t iTri) const;


  private:
    void checkConnectivity();
    void evalTriangles(mapfacetype bound_faces, size_t & nbTri, bool outwardNormOnBoundary = false);
    void writePoints(std::string outputFileName,double rescaling=1.0, bool binary=false);
    void writeElements(std::string outputFileName, bool binary=false);
    void writeFakeFiebers(std::string outputFileName, bool binary=false);
    bool isLittleEndian();
    void SwapBytes(void *pv, size_t n);
    std::vector<double>InvertA3X3(const std::vector<double> & Mat0) const;
    short int RM3X3Ind(short int irow, short int jcol) const;
    bool consistentState;
    std::vector<Point> points;
    std::vector<Triangle> triangles;
    std::vector<Tetrahedron> tetrahedra;
    std::vector<size_t> triaToTet;
    MeshInfo info;
    // map between a region and the list of points belonging to it
    // two regions could share the same point if on boundary
    regionSubdivisionType pointRegions;
    regionSetType regionLabels;
    std::multimap<size_t,int> nbElToRegionLab;
    std::set<long int> _Endo;
    std::set<long int> _Epi;
    int _endoLabel;
    int _epiLabel;
    ElemToBoundaryFaceConnectionType _faceToFace;
    ElemToBoundaryFaceConnectionType _elemBoundaryTris;
    facetype _endoTris;
    facetype _epiTris;
};

#endif

