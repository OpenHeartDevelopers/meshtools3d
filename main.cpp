#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cfloat>
#include "Chrono.hpp"
#include "Mesh.hpp"
#include "LaplaceSolver.hpp"
#include "ThicknessEvaluation.hpp"
#include "VtkWriter.hpp"
#include "GetPot.hpp"
#include<CGALDataType.hpp>
#include<set>
#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/task_scheduler_init.h>
#endif


int main(int argc,char **argv)
{

  GetPot command_line(argc,argv);
	std::string data_file_name = command_line.follow("data", 2,  "-f","--file");
	GetPot param_file(data_file_name.c_str());
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
  bool eval_thickness        = param_file("others/eval_thickness",false); 
  bool out_medit             = param_file("output/out_medit",false);
  bool out_carp              = param_file("output/out_carp",false);  
  bool out_carp_binary       = param_file("output/out_carp_binary",false);
  bool out_vtk               = param_file("output/out_vtk",false);
  bool out_vtk_binary        = param_file("output/out_vtk_binary",false);
  bool out_potential         = param_file("output/out_potential",false);
  double rescaling           = param_file("meshing/rescaleFactor",1.0);
  

#ifdef CGAL_LINKED_WITH_TBB
  int numThreads = 1;
  char * nthr=NULL;
  nthr = getenv("TBB_NUM_THREADS");
  if (nthr != NULL)
    numThreads = atoi(nthr);
  if ((nthr == NULL) || (numThreads == 0)) {
    numThreads = 1;
    std::cout << "TBB_NUM_THREADS not set; nb of threads is: "
              << numThreads
              << " (default)"
              << std::endl;
  } else {
    std::cout << "nb of threads is: " << numThreads << std::endl;
  }
  tbb::task_scheduler_init init(numThreads);
#endif

  Mesh CarpMesh;
//create the output dir
  
  std::string command="mkdir -p "+out_dir;
  int ierr=system(command.c_str());
  if(!ierr==0)
  {
    std::cerr<<"Problem in creating the directory "<<out_dir<<std::endl;
    exit(1);
  }
  else
  {
    Facet_criteria facet_criteria(fangle, fsize, fapprox); 
    Cell_criteria cell_criteria(cR_E_ratio, csize); 
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
      chrono.reset();

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
      size_t nbTet=static_cast<size_t>(c3t3.number_of_cells_in_complex());
      if(nbPt<3 ||nbTet<1)
      {
        std::cerr<<"Problem with Triangulations, only "<<nbPt<<" Vertices and "<<nbTet<<"Thetraedra"<<std::endl;
        std::string mfileoutName=out_dir+"/"+out_name+".mesh";
        std::ofstream medit_file(mfileoutName.c_str());
        c3t3.output_to_medit(medit_file);
        medit_file.close();
        exit(1);
      }
      else
      {
        std::cout<<"Vertices: "<<nbPt<<std::endl;
        std::cout<<"Tetra:    "<<nbTet<<std::endl;
      }
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
          CarpMesh.Info().checksum[ic]=CarpMesh.Info().checksum[ic]+static_cast<float>(CarpMesh.Pt(icount).coord[ic]);

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
  }
  std::cout<<"pre-processing..."<<std::flush;
  Chrono chrono2;
  chrono2.start();
  CarpMesh.preprocessingOperations();
  chrono2.stop();
  std::cout<<" done in "<<chrono2<<std::endl;
  chrono2.reset();
  
  std::cout<<"boundary extraction..."<<std::flush;
  chrono2.start();
  CarpMesh.extractBoundary();
  chrono2.stop();
  std::cout<<" done in "<<chrono2<<std::endl;
  chrono2.reset();

  std::cout<<"boundary re-labeling..."<<std::flush;
  chrono2.start();
  CarpMesh.evalBoundaryLabels();
  chrono2.stop();
  std::cout<<" done in "<<chrono2<<std::endl;
  chrono2.reset();

  CarpMesh.writeBoundaryLabels(out_dir, out_name);
  CarpMesh.meshRescaling(rescaling); //rescale the whole mesh; not on output. So output of rescaling is 1 now
  if(out_carp)
  {
      std::string cfileoutName=out_dir+"/"+out_name;
      CarpMesh.writeCarpMesh(cfileoutName,out_carp_binary);
  }

  
  VtkWriter writerVTK(& CarpMesh, out_vtk_binary);
  if(out_vtk)
  {
    writerVTK.setOutputDir(out_dir);
    writerVTK.setPrefixName(out_name);
    writerVTK.openFileForOutput();
    std::vector<double> meshLabels = CarpMesh.copyLabelVectorForVTKOutput();
    writerVTK.writeVariable(meshLabels, "region_labels",VtkWriter::Scalar);
    meshLabels.clear();
  }

  //here I free some variables
  CarpMesh.unsetBoundaryLabels();
  
  if(eval_thickness)
  {
    
    if(!(param_file("others/thickalgo",1)==2))
    {
      CarpMesh.initializeConnectivities();
    }
    
    ThicknessEvaluation ThicknessCompute(param_file, &CarpMesh );
    ThicknessCompute.setBCValue(CarpMesh.Endocardium(), 0.0);  
    ThicknessCompute.setBCValue(CarpMesh.Epicardium(), 1.0);  
    ThicknessCompute.solve();
    std::string cfileoutName=out_dir+"/"+out_name;
    ThicknessCompute.writeElementGradient(cfileoutName);
    CarpMesh.writeTetraCentroids(cfileoutName);
    CarpMesh.writeTris(cfileoutName);
    
    if(out_potential)
    {
      if(out_vtk)
      {
        writerVTK.writeVariable(ThicknessCompute.sol(), "potential_func",VtkWriter::Scalar);
      }
      else
      {
        if(out_carp || !out_vtk) 
        {
          ThicknessCompute.writeSolution(out_dir+"/"+out_name);
        }
      }
    }
    
    if(out_vtk)
    {
      writerVTK.writeVariable(ThicknessCompute.thickness(), "Thickness",VtkWriter::Scalar);
    }
   
   
   if(out_carp || !out_vtk) 
   {
      ThicknessCompute.writeThickness(out_dir+"/"+out_name);
   }
    
  }
  
  
  if(out_vtk)
  {
    writerVTK.CloseFile();
  }
  
  return 0;
  
  
}
