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
#include "INRreader.hpp"
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
  
  if( command_line.search(2, "-h","--help") or (argc==1) ) 
  {
      std::cout<<"MeshTools3D\nUsage: ";
      std::cout<<"./meshtools3d -f <data_filename> -seg_dir <segmentation_dir> -seg_name <segmentation_name> -out_dir <output_dir> -out_name <output_name>"<<std::endl;
      std::cout<<"\t data_filename is the file name where parameters are"<<std::endl;
      std::cout<<"\t segmentation_dir, segmentation_name are the directory and the name of the input segmentation"<<std::endl;
      std::cout<<"\t output_dir, output_name are the output directory and the output suffix"<<std::endl;
      std::cout<<"NOTE:"<<std::endl;
      std::cout<<"\t If -seg_dir, -seg_name, -out_dir, -out_name are optional aruments,\n\t they overwrite those specified in the data file (used, for example, to use meshtools3d within a script)"<<std::endl;
      exit(0);
  }
  
  
  std::string seg_dir            = param_file("segmentation/seg_dir",".");
  std::string seg_name           = param_file("segmentation/seg_name","image.inr");
  bool mesh_from_segmentation    = param_file("segmentation/mesh_from_segmentation",true);

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
  bool debug_output          = param_file("output/debug_output",false);
  size_t debug_frequency     = param_file("output/debug_frequency",100000);
  
  double rescaling           = param_file("meshing/rescaleFactor",1.0);

  seg_dir  = command_line.follow(seg_dir, 2,  "-seg_dir","--segmentation_directory");
  seg_name = command_line.follow(seg_name, 2,  "-seg_name","--segmentation_name");

  out_dir  = command_line.follow(out_dir, 2,  "-out_dir","--output_directory");
  out_name = command_line.follow(out_name, 2,  "-out_name","--output_name");


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
  nthr=NULL;
#endif

  Mesh CarpMesh;
  INRreader Segmentation;

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
      std::string imageName=seg_dir+"/"+seg_name;
      Chrono chrono;    
      if(mesh_from_segmentation)   // CGAL reads segmentation with labels and mesh it
      {
          C3t3 c3t3;
          CGAL::Image_3 image1;
          if(image1.read(imageName.c_str()))
          {
              Mesh_domain domain(image1);
              Facet_criteria facet_criteria(fangle, fsize, fapprox); 
              Cell_criteria cell_criteria(cR_E_ratio, csize);
              Mesh_criteria criteria(facet_criteria, cell_criteria);
              std::cout<<"MESHING...";
              chrono.start();
              c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, 
                                               CGAL::parameters::lloyd(), CGAL::parameters::odt(), 
                                               CGAL::parameters::perturb(), CGAL::parameters::exude());
              chrono.stop();
              std::cout<<" done in "<<chrono<<std::endl;
              chrono.reset();
          }
          else
          {
                std::cerr<<"ERROR: IMAGE NOT READ"<<std::endl;
                exit(1);
          }
          if(out_medit) // print mesh in .mesh format
          {
              std::string mfileoutName=out_dir+"/"+out_name+".mesh";
              std::ofstream medit_file(mfileoutName.c_str());
              c3t3.output_to_medit(medit_file);
              medit_file.close();
          }
          //get the triangulation; 
          // Here I need a reference, otherwise when I go to take the tetra vertex
          // these will be not consistent (most of them 0)
          const Tr & triang= c3t3.triangulation();  
          size_t nbPt=static_cast<size_t>(triang.number_of_vertices());
          size_t nbTet=static_cast<size_t>(c3t3.number_of_cells_in_complex());
          if(nbPt<3 ||nbTet<1)
          {
              std::cerr<<"Problem with Triangulations, only "<<nbPt<<" Vertices and "<<nbTet<<" Thetraedra"<<std::endl;
              if(!out_medit)
              {
                  std::string mfileoutName=out_dir+"/"+out_name+".mesh";
                  std::ofstream medit_file(mfileoutName.c_str());
                  c3t3.output_to_medit(medit_file);
                  medit_file.close();
              }
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
      } // Segmentation is read by INRreader; mesh domain implemented through a function wrapper with no labels
      else
      {
          Segmentation.readSegmentation(imageName);
          Mesh_domain_manualseg domain(Segmentation);
          Facet_criteria_manualseg facet_criteria(fangle, fsize, fapprox); 
          Cell_criteria_manualseg  cell_criteria(cR_E_ratio, csize);
          Mesh_criteria_manualseg  criteria(facet_criteria, cell_criteria);
          std::cout<<"MESHING...";
          chrono.start();
          C3t3_manualseg c3t3= CGAL::make_mesh_3<C3t3_manualseg>(domain, criteria, 
                                              CGAL::parameters::lloyd(), CGAL::parameters::odt(), 
                                              CGAL::parameters::perturb(), CGAL::parameters::exude());
          chrono.stop();
          std::cout<<" done in "<<chrono<<std::endl;
          chrono.reset();
          if(out_medit) // print mesh in .mesh format
          {
              std::string mfileoutName=out_dir+"/"+out_name+".mesh";
              std::ofstream medit_file(mfileoutName.c_str());
              c3t3.output_to_medit(medit_file);
              medit_file.close();
          }
          //get the triangulation; 
          // Here I need a reference, otherwise when I go to take the tetra vertex
          // these will be not consistent (most of them 0)
          //get the triangulation; 
          // Here I need a reference, otherwise when I go to take the tetra vertex
          // these will be not consistent (most of them 0)
          const Tr_manualseg & triang= c3t3.triangulation();  
          size_t nbPt=static_cast<size_t>(triang.number_of_vertices());
          size_t nbTet=static_cast<size_t>(c3t3.number_of_cells_in_complex());
          if(nbPt<3 ||nbTet<1)
          {
              std::cerr<<"Problem with Triangulations, only "<<nbPt<<" Vertices and "<<nbTet<<" Thetraedra"<<std::endl;
              if(!out_medit)
              {
                  std::string mfileoutName=out_dir+"/"+out_name+".mesh";
                  std::ofstream medit_file(mfileoutName.c_str());
                  c3t3.output_to_medit(medit_file);
                  medit_file.close();
              }
              exit(1);
          }
          else
          {
              std::cout<<"Vertices: "<<nbPt<<std::endl;
              std::cout<<"Tetra:    "<<nbTet<<std::endl;
          }
          CarpMesh.initializePtVector(nbPt);
          size_t icount=0;
          std::map<Tr_manualseg::Vertex_handle, size_t> Vertices;
          for (Tr_manualseg::Finite_vertices_iterator it=triang.finite_vertices_begin(); it != triang.finite_vertices_end(); ++it)
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
          for(Cell_iterator_manualseg itc = c3t3.cells_in_complex_begin(); itc != c3t3.cells_in_complex_end(); ++itc)
          {
            for(int iv=0; iv<4; iv++)
            {
              size_t vindex= Vertices.at(itc->vertex(iv)) ;
              CarpMesh.Tet(icount).vertex[iv]=vindex;
            }
            //CarpMesh.Tet(icount).regionLabel=static_cast<int>(itc->subdomain_index());
            CarpMesh.Tet(icount).regionLabel=static_cast<int>(1); //label =1 by default: this is used as myocardium value
            icount++;
          }
      }
  }//end if on mkdir
  std::cout<<"pre-processing..."<<std::flush;
  Chrono chrono2;
  chrono2.start();
  CarpMesh.preprocessingOperations();
  chrono2.stop();
  std::cout<<" done in "<<chrono2<<std::endl;
  chrono2.reset();
  
  std::cout<<"boundary extraction..."<<std::flush;
  chrono2.start();
  if(!mesh_from_segmentation)
  {
    CarpMesh.extractBoundary(true);
  }
  else
  {
    CarpMesh.extractBoundary();
  }
  
  chrono2.stop();
  std::cout<<" done in "<<chrono2<<std::endl;
  chrono2.reset();

  if(!mesh_from_segmentation)
  {
    // HERE GOES RELABELING 
    // NB: AT this point algorthm expects that tria are defined
    // and that have the same label of the tetra they belongs to
    // If you are here, all of your tetra have a label = 1, since
    // you used a un-labeled segmentation
    // Perhaps it is worth to discard tetra label at this point
    std::cout<<"Re-apply labels on triangles"<<std::endl;
    chrono2.start();
    Segmentation.createBoundingBoxes();
    const INRreader::BboxMapType & bboxmap=Segmentation.bboxlabels();
    for(INRreader::BboxMapCIterType mapit=bboxmap.begin(); mapit!=bboxmap.end(); ++mapit)
    {
      std::cout<<"Relabeling region "<<mapit->first<<std::endl;
      Mesh::IDsetType candidateSet=CarpMesh.extractTrianglesFromBBOX((mapit->second).bbox());
      for(Mesh::IDiteratorType tr_it=candidateSet.begin(); tr_it!=candidateSet.end(); ++tr_it)
      {
          std::vector<double> cgT=CarpMesh.TriaCentroid(*tr_it);
          std::vector<double> voxval=Segmentation.interpolatedNonZeroVoxelValue(cgT);
          int regionLabel=static_cast<int>(voxval[0]);
          CarpMesh.Tri(*tr_it).regionLabel=regionLabel;
      }
    }
    chrono2.stop();
    std::cout<<" done in "<<chrono2<<std::endl;
  }
  CarpMesh.meshRescaling(rescaling); //rescale the whole mesh; not on output. So output of rescaling is 1 now
  std::cout<<"boundary re-labeling..."<<std::flush;
  chrono2.start();
  CarpMesh.evalBoundaryLabels(debug_output,out_dir,debug_frequency);
  chrono2.stop();
  std::cout<<" done in "<<chrono2<<std::endl;
  chrono2.reset();

  CarpMesh.writeBoundaryLabels(out_dir, out_name);
  
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
    
    ThicknessEvaluation ThicknessCompute(param_file, &CarpMesh );
    
    ThicknessCompute.setBCValue(CarpMesh.Endocardium(), 0.0);  
    ThicknessCompute.setBCValue(CarpMesh.Epicardium(), 1.0);  
    
    ThicknessCompute.solve();
    if(ThicknessCompute.algorithm()!=static_cast<unsigned char>(2)  )
    {
      CarpMesh.initializeConnectivities();
    }
    
    ThicknessCompute.evalThickness();
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
  }// end if on eval thickness
  
  if(out_vtk)
  {
    writerVTK.CloseFile();
  }
  
  return 0;
  
  
}
