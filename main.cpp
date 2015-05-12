#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "Chrono.hpp"
#include "Mesh.hpp"
#include "GetPot.hpp"
#include<CGALDataType.hpp>
#include<set>

int main(int argc,char **argv)
{

  GetPot command_line(argc,argv);
	const char* data_file_name = command_line.follow("data", 2,  "-f","--file");
	GetPot param_file(data_file_name);
	if( command_line.search(2, "-i","--info") ) 
	{
		param_file.print();
		exit(0);
	}
  std::string seg_dir            = param_file("segmentation/seg_dir",".");
  std::string seg_name           = param_file("segmentation/seg_name","image.inr");

  FaceNumericalType fangle      = param_file("meshing/facet_angle",30);
  FaceNumericalType fsize       = param_file("meshing/facet_size",0.8);
  FaceNumericalType fapprox     = param_file("meshing/facet_distance",4);
  CellNumericalType cR_E_ratio  = param_file("meshing/cell_rad_edge_ratio",2.0);
  CellNumericalType csize       = param_file("meshing/cell_size",1.0);

  std::string out_dir        = param_file("output/outdir","."); 
  std::string out_name       = param_file("output/name","imgmesh"); 
  bool out_binary            = param_file("output/out_binary",false);
  bool out_medit             = param_file("output/out_medit",false);
  bool out_carp              = param_file("output/out_carp",false);  
  bool out_vtk               = param_file("output/out_vtk",false);
  double rescaling           = param_file("meshing/rescaleFactor",1.0);
  
  Facet_criteria facet_criteria(fangle, fsize, fapprox); 
  Cell_criteria cell_criteria(cR_E_ratio, csize); 
  
  Mesh CarpMesh;
  std::string imageName=seg_dir+"/"+seg_name;
  CGAL::Image_3 image1;
  if(image1.read(imageName.c_str()))
  {
    Mesh_domain domain(image1);
    // Set mesh criteria
    Mesh_criteria criteria(facet_criteria, cell_criteria);
    std::cout<<"MESHING...";
    
    Chrono chrono;
    chrono.start();

    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, 
                                           CGAL::parameters::lloyd(), CGAL::parameters::odt(), 
                                           CGAL::parameters::perturb(), CGAL::parameters::exude());
    chrono.stop();
    std::cout<<" done in "<<chrono<<std::endl;
    if(out_medit)
    {
      std::string mfileoutName=out_dir+"/"+out_name+".mesh";
      std::ofstream medit_file(mfileoutName.c_str());
      c3t3.output_to_medit(medit_file);
      medit_file.close();
    }
    //get the triangulation; 
    //Here I need a reference, otherwise when I go to take the tetra vertex
    // these will be not consistent (most of them 0)
    
    const Tr & triang= c3t3.triangulation();  
    size_t nbPt=static_cast<size_t>(triang.number_of_vertices());
    CarpMesh.initializePtVector(nbPt);
    size_t icount=0;
    std::map<Tr::Vertex_handle, size_t> Vertices;
    for (Tr::Finite_vertices_iterator it=triang.finite_vertices_begin(); it != triang.finite_vertices_end(); ++it)
    {

	    Vertices[it]=icount;
	    
	    CarpMesh.Pt(icount).x()=it->point().x();
	    CarpMesh.Pt(icount).y()=it->point().y();
	    CarpMesh.Pt(icount).z()=it->point().z();
     
      for(short int ic=0; ic<3; ic++)
      {
        CarpMesh.Info().baricenter[ic]=CarpMesh.Info().baricenter[ic]+(CarpMesh.Pt(icount).coord[ic]);

        if(CarpMesh.Info().bbox[ic][0]>CarpMesh.Pt(icount).coord[ic])
        {
          CarpMesh.Info().bbox[ic][0]=CarpMesh.Pt(icount).coord[ic];
        }
        if(CarpMesh.Info().bbox[ic][1]<CarpMesh.Pt(icount).coord[ic])
        {
          CarpMesh.Info().bbox[ic][1]=CarpMesh.Pt(icount).coord[ic];
        }
      }
      ++icount;
    }
    for(short int ic=0; ic<3; ic++)
    {
      CarpMesh.Info().baricenter[ic]=CarpMesh.Info().baricenter[ic]/nbPt;
    }

//TETRA
    size_t nbTet=static_cast<size_t>(c3t3.number_of_cells_in_complex());
    
    CarpMesh.initializeTetraVector(nbTet);
    icount=0;
    for(Cell_iterator itc = c3t3.cells_in_complex_begin(); itc != c3t3.cells_in_complex_end(); ++itc)
    {
      for(int iv=0; iv<4; iv++)
      {
        size_t vindex= Vertices.at(itc->vertex(iv)) ;
        
        CarpMesh.Tet(icount).vertex[iv]=vindex;
      }
      CarpMesh.Tet(icount).regionLabel=static_cast<int>(itc->subdomain_index());
      icount++;
    }
  }
  else
  {
    std::cerr<<"ERROR: IMAGE NOT READ"<<std::endl;
    exit(1);
  }
 
  CarpMesh.extractBoundary();
  //CarpMesh.evalBoundaryLabels();
  //CarpMesh.writeBoundaryLabels(out_dir, out_name);


  if(out_carp)
  {
      std::string cfileoutName=out_dir+"/"+out_name;
      CarpMesh.writeCarpMesh(cfileoutName,rescaling);
  }
  if(out_vtk)
  {
      std::string cfileoutName=out_dir+"/"+out_name;
      CarpMesh.writeVTKMesh(cfileoutName,rescaling,out_binary);
  }
  return 0;
  
  
}
