#include "ThicknessEvaluation.hpp"
#include<iostream>
#include "Chrono.hpp"
#include "mgmres.hpp"
#include<cmath>

#ifndef MAX_PATH_LEN
#define MAX_PATH_LEN 50000
#endif

#ifndef IN_ELEM_TOLL
#define IN_ELEM_TOLL 1.0
#endif

#ifndef MAX_ITER_THICK
#define MAX_ITER_THICK 1000
#endif

ThicknessEvaluation::ThicknessEvaluation()
:LaplaceSolver()
{}

ThicknessEvaluation::ThicknessEvaluation(const Mesh *  _mesh)
:LaplaceSolver(_mesh)
{
  _thickness.resize(_ptrmesh->nPt(),0);
}

void ThicknessEvaluation::setMesh(const Mesh * _mesh)
{
  LaplaceSolver::setMesh(_mesh);
  _thickness.resize(_ptrmesh->nPt(),0);
}

ThicknessEvaluation::ThicknessEvaluation(const GetPot & dfile, const Mesh * _mesh)
:LaplaceSolver(dfile, _mesh)
{
  _thickness.resize(_ptrmesh->nPt(),0);
}

ThicknessEvaluation::~ThicknessEvaluation()
{
  _thickness.clear();
}

void ThicknessEvaluation::evalThickness()
{
  typedef std::set<size_t> facetype;
  typedef std::map<short int, size_t> boundaryFaceInTetraType;
  
  int interactions=0;
	int problems=0;
	std::vector<double> path(_ptrmesh->nTri(),0);
  
  const facetype & illumTris  = _ptrmesh->endoTria();
  const facetype & groundTris = _ptrmesh->epiTria();
  const std::vector<boundaryFaceInTetraType> & elemBoundaryTris =  _ptrmesh->elemBoundaryTris();
  const std::vector<boundaryFaceInTetraType> & faceToFace =  _ptrmesh->faceToFace();
  const std::vector<Tetrahedron> & elems = _ptrmesh->Tet();
  const std::vector<Triangle> & tris = _ptrmesh->Tri();
  const std::vector<Point> & coords = _ptrmesh->Pt();
  Chrono chrono;
  std::cout << "*******************************"<<std::endl;
  std::cout << "* Computing wall thickness... *"<<std::endl;
  std::cout << "*******************************"<<std::endl;
  chrono.start();
  for(facetype::iterator it=illumTris.begin();it!=illumTris.end();it++)
	{
    // Picks-out particular triangle
		size_t thisIllumTri = *it;

		// Picks-out the TISSUE element associated with this surface triangle from which we are illuminating
	  size_t surfaceElem = _ptrmesh->triaToTetMap(thisIllumTri);

		/////////////////////////////////////////////////////////////////////////////////////
		// Iterates over all nodes of selected triangle and derives centroid (of triangle) //
		/////////////////////////////////////////////////////////////////////////////////////
	  std::vector<double> x_coord_epi(3,0);
	  
	  for(short int jVertex=0;jVertex<3;jVertex++)
		{
			// One of the 3 triangle nodes
			size_t node = tris[thisIllumTri].vertex[jVertex];
			// Calculates centroid of triangle
			for(short int iCoord=0; iCoord<3; iCoord++)
			{
			  x_coord_epi[iCoord] +=coords[node].coord[iCoord]/3.0;
			}
		}
		
		// Set-ups list to store the 4 nodes in each element        	
		facetype elemNodes;
  	facetype::iterator i2;

		// Sets-up list to store the previous triangle nodes from th face we just passed through
		facetype faceRem;
		facetype::iterator it_fr;

		// Initialises this list with the nodes of the epicardial triangle surface face
		faceRem.clear();
  	for(short int jVertex=0;jVertex<3;jVertex++)
    {
      faceRem.insert(tris[thisIllumTri].vertex[jVertex]);
    }
		
		// Initial coordinates of photon packet set as centroid of triangle
		std::vector<double> x_coord = x_coord_epi;
		std::vector<double> x_coord_s(3,0);
		
		// Defines initial element we're in
		size_t element = surfaceElem;
		size_t oldElement = element;
    
    // Define initial direction as vector gradient direction
    std::vector<double> mu = ElementTetraGradient(element);
    std::vector<double> mu_new=mu;
    
    // Defines some counters and flags
    bool inTissue = true;
    int newCounter = 0;
    bool adjustMu = false;
    bool wrongSurface = false;
    while(inTissue)
    {
      // Counts total interactions
			interactions++;				
      // Checks path to see if it's too long
      if(path[thisIllumTri] > MAX_PATH_LEN)
      {
        std::cerr<<"WARNING: Triangle  "<<thisIllumTri<<" MAX LEN ("<<MAX_PATH_LEN<<") PASSED: len =  "<<path[thisIllumTri]<<std::endl;
        break;
      }

      /////////////////////////////////////////////////////////////////////
      // Quick check to see if we are not in the element we should be in //
      /////////////////////////////////////////////////////////////////////

      // this part is not used after (if needed better to re-implement with affine transfo)
      bool notInElement=!(isInElement(x_coord ,element));

      newCounter++;
      
      if(newCounter > MAX_ITER_THICK)
      {
        inTissue = false;
      }

      if( (adjustMu == false) && (wrongSurface == false))
			{
					// Initial directional cosines set from normal of triangle (calculating dot-product of each axis with normal)
					mu = ElementTetraGradient(element);
			}
			else 
			{
					mu = mu_new;
			}

      // Updates starting position of step
      x_coord_s = x_coord;
      wrongSurface = false;			
      
      //////////////////////////////////////////////////////////
      // Check for intersection with faces of current element //
      //////////////////////////////////////////////////////////
      
      // Defines things to define closest face
      //double aMin = DBL_MAX;
      double aMin = 1000000.0;
      
      short int intFace = -1;
      
      // Creates a set with all nodes defining element
      elemNodes.clear();
      for(short int iVertex=0;iVertex<4;iVertex++)
      {
        elemNodes.insert(elems[element].vertex[iVertex]);      
      }

      std::vector<double> Nmin(3,0);

      // Iterates over each triangle face
      for(short int iVertex = 0;iVertex < 4; iVertex++)
      {
					// Erases one node (which then defines a triangle)
					elemNodes.erase(elems[element].vertex[iVertex]);
					// Defines each of 3 remaining nodes which now define the triangle face
				  std::vector<size_t> triNodes(3,0);
					short int p = 0; 
					for(i2=elemNodes.begin();i2!=elemNodes.end();i2++)
					{
						triNodes[p] = *i2;
						p++;
					}
					// Checks to see if this is the same face and we've just pass through (if it is, don't include it)
					short int c = 0;
					for(short int m=0;m<3;m++)
					{
						it_fr = faceRem.find(triNodes[m]);
						if(it_fr != faceRem.end())
						{
						  c++;
						}
					}
					// Only if this isn't the face we've just interacted with do we now calculate an a value
					if(c < 3)
					{
						// Calculates inward-facing normal of this triangle
            
            std::vector<double> centroid = _ptrmesh->TetraCentroid(element);
            GrahmOperatorOutput	geoQuant = GrahmOperations( coords[triNodes[0]],  coords[triNodes[1]],  coords[triNodes[2]], centroid);
            centroid.clear();
            double a = distanceOfPointToPlane(geoQuant.N,  coords[triNodes[0]], x_coord_s,  mu);
						// Calculates intersection distance
						if((a >= 0) && (a<aMin))
						{
							// Updates smallest intersection distance
							aMin = a;
							// Defines closest intersecting face
							intFace = iVertex;
							Nmin=geoQuant.N;
						}
						// Calculates NORMAL distance of point to plane
            double a_n =normalDistanceOfPointToPlane(geoQuant.N,  coords[triNodes[0]], x_coord_s);
            //double dot = N[0]*grads[element][0] + N[1]*grads[element][1] + N[2]*grads[element][2];
            if(a_n > 0)
            {
              problems++;
            }
						while(a_n > 0)
						{
							  // Move along line towards centroid by a small amount and recompute a
							  double delta = 1.0;
							  centroid = _ptrmesh->TetraCentroid(element);
							  std::vector<double> v(3,0);
							  v =  waxpy(x_coord, centroid,-1.0);
							  centroid.clear();
							  double magV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
							  for(short int hCoord=0; hCoord<3; hCoord++)
							  {
							    v[hCoord] = v[hCoord]/magV;
							  }
							  //x=x+dlta*v
							  for(short int hCoord=0; hCoord<3; hCoord++)
							  {
							    x_coord[hCoord] =x_coord[hCoord] + delta*v[hCoord];	
							  }
							  
                a_n = normalDistanceOfPointToPlane(geoQuant.N,coords[triNodes[0]],x_coord);
                a   = distanceOfPointToPlane(geoQuant.N,coords[triNodes[0]],x_coord, mu);
							  if(a >= 0 && a<aMin)
                {
                  // Updates smallest intersection distance
                  aMin = a;
                  // Defines closest intersecting face
                  intFace = iVertex;
                  Nmin=geoQuant.N;
                }
						}
					}
					// Puts the erased node back
					elemNodes.insert(elems[element].vertex[iVertex]);
				}
				
				if(intFace != -1)
				{
          //lastintface = intFace;
					//lastaMin = aMin;
					// Add details of face we've just moved to
					elemNodes.erase(elems[element].vertex[intFace]);
          // Defines each of 3 remaining nodes which now define the triangle face
					faceRem.clear();
					for(i2=elemNodes.begin();i2!=elemNodes.end();i2++)
					{
					  faceRem.insert(*i2);
					}
					path[thisIllumTri] += aMin;	

					// Checks to see if this closest intersecting face is also a boundary face
					// This is for interior boundaries
					const boundaryFaceInTetraType & ElementBTris = elemBoundaryTris[element];
					boundaryFaceInTetraType::const_iterator iterBoundaryOnFace = ElementBTris.find(intFace);
					
          if(iterBoundaryOnFace != ElementBTris.end() )
          {
            // Checks to see if this is part of the ground electrode
            if(groundTris.find(iterBoundaryOnFace->second) != groundTris.end())
            {
                inTissue = false;
                break;
            }
            else
            {
              wrongSurface = true;
							// Move a small distance in towards the centroid
							std::vector<double> centroid = _ptrmesh->TetraCentroid(element);
							double delta = 1.0;
              std::vector<double> v(3,0);
              v =  waxpy(x_coord,centroid,-1.0);
              centroid.clear();
              double magV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
              for(short int hCoord=0; hCoord<3; hCoord++)
							{
							  v[hCoord]       = v[hCoord]/magV;
							  x_coord[hCoord] = x_coord[hCoord]+ delta*v[hCoord];	
							}

              // Calculate the normal to the current gradient and the face normal
              std::vector<double> w(3,0),u(3,0);
              u=ElementTetraGradient(element);
              
              std::vector<double> N1(3,0);
              N1[0] = (u[1]*Nmin[2] - u[2]*Nmin[1]);
              N1[1] = (u[2]*Nmin[0] - u[0]*Nmin[2]);
              N1[2] = (u[0]*Nmin[1] - u[1]*Nmin[0]);
							
							// Calculate the normal to the face normal and this normal to get the vector parallel to the face
							w[0] = (N1[1]*Nmin[2] - N1[2]*Nmin[1]);
							w[1] = (N1[2]*Nmin[0] - N1[0]*Nmin[2]);
							w[2] = (N1[0]*Nmin[1] - N1[1]*Nmin[0]);
              double magW = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
              for(short int hCoord=0; hCoord<3; hCoord++)
              {
                  w[hCoord] = w[hCoord]/magW;
              }
              
              // Check using the dot product that this new vector is in the same direction as the previous gradient
              double dot3 = w[0]*u[0] + w[1]*u[1] + w[2]*u[2];
              if(dot3 < 0)
              {
                for(short int hCoord=0; hCoord<3; hCoord++)
                {
                  w[hCoord] = -1.0*w[hCoord];
                }
              }
              
              mu_new = w;
            } //end else (face boundary but not ground)
          } //end if on face boundary

					// We also update our 'current element' as the one we're moving into as we cross this boundary
          const boundaryFaceInTetraType & faceToFaceMap = faceToFace[element];
          boundaryFaceInTetraType::const_iterator itfaceToFaceMap = faceToFaceMap.find(intFace);
          
          // Only if we haven't passed through a boundary
          if(itfaceToFaceMap != faceToFaceMap.end())
          {
            size_t newElement = itfaceToFaceMap->second;
						oldElement = element;
						element = newElement;
						
						mu_new = ElementTetraGradient(newElement);
            if( (Nmin[0]*mu_new[0] + Nmin[1]*mu_new[1] + Nmin[2]*mu_new[2]) >= 0)
						{	
							  adjustMu = true;
							  mu_new = ElementTetraGradient(newElement);
							  std::vector<double> mu_old = ElementTetraGradient(oldElement);
							  for(short int hCoord=0; hCoord<3; hCoord++)
							  {
							    mu_new[hCoord] = 0.5*(mu_new[hCoord] + mu_old[hCoord]);
							  }
							  
							  if((Nmin[0]*mu_new[0] + Nmin[1]*mu_new[1] + Nmin[2]*mu_new[2]) >= 0)
							  {
							  	  element = oldElement;
							  }
						}
						else
						{
						    adjustMu = false;
						}
            // Move-up to face of element and update position
            
            //x = x_s +aMin*mu
            x_coord=waxpy(mu, x_coord_s,aMin);
						
						
						////////////////////////////////////////////////////////////////
						// THE NEXT LINES ARE USELESS SINCE delta=0.0; I COMMENT THEM //
						////////////////////////////////////////////////////////////////
						/*
						// Instead of this, move a tiny amount into the new element
						std::vector<double> centroid = _ptrmesh->TetraCentroid(newElement);
						std::vector<double> v(3,0);
						v = waxpy(x_coord, centroid,-1.0);
						centroid.clear();
						double vMag = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2]  );
						double delta = 0.0;
						//x = x + (delta/vMag)*v;
						x_coord= waxpy(v, x_coord,delta/vMag);
						*/
						
          }
				}// end if on intface!=-1
    }//end while on inTissue
	}//end loop on endoTria

  chrono.stop();
  std::cout << "Done in "<<chrono<<std::endl;
  std::cout << "Total number of interactions = " << interactions << std::endl;
  double propProb = double(problems)/double(interactions);
  std::cout << "Total number of problems = " << problems << " ( or " << propProb << " ) "<<std::endl;

// I commented the following instructions since elemData is not used thereafter  
  /*
  std::vector<double> elemData(_ptrmesh->nTet(),0);
  for(facetype::iterator it=illumTris.begin();it!=illumTris.end();it++)
	{
	  size_t elem_n = _ptrmesh->triaToTetMap(*it);
	  size_t tri_n = *it;
	  elemData[elem_n]=path[tri_n];
  }
  */
  
  // Maps data over to node list keeping track of how many surface triangles are associated with each node 
  _thickness.clear();
  _thickness.resize(_ptrmesh->nPt(),0);
  
  std::vector<int> surfPlotterCounter(_ptrmesh->nPt(),0);
  
  
  for(facetype::iterator it=illumTris.begin();it!=illumTris.end();it++)
  {
	  size_t tri_n=*it;
	  for(short int jVertex=0;jVertex<3;jVertex++)
	  {
	  	size_t node_n = (tris[tri_n]).vertex[jVertex];
	  	_thickness[node_n] += path[tri_n];
	  	surfPlotterCounter[node_n]++;
	  }
  } 

  for(size_t iPt=0;iPt<_thickness.size();iPt++)
  {
  	if(surfPlotterCounter[iPt] != 0)
  	{
  	  _thickness[iPt] = _thickness[iPt]/double(surfPlotterCounter[iPt]);
  	}
  }

}

void ThicknessEvaluation::solve()
{
  LaplaceSolver::solve();
  evalThickness();
}



void ThicknessEvaluation::writeThickness(std::string filename)
{
  std::string fname=filename+"_thickness.dat";
  std::ofstream fsol(fname.c_str());
  if(!fsol)
  {
    std::cerr<<"ERROR: FILE "<<fname<<" NOT OPENED"<<std::endl;
    exit(1);
  }
  size_t nPt=_ptrmesh->nPt();
  for(size_t iPt=0; iPt<nPt; iPt++)
  {
    fsol<<_thickness[iPt]<<std::endl;
  }
  fsol.close();


}


bool ThicknessEvaluation::isInElement(const std::vector<double> & xc ,const size_t & iElem)
{
  // Tells if a point of coordinate xc is within element iElem
  // not in element is not the output of this function
  bool notInElement=false;
  const Tetrahedron & Tetra = _ptrmesh->Tet(iElem);
  const std::vector<Point> & coords = _ptrmesh->Pt();
  for(short int qCoord=0;qCoord<3;qCoord++)
  {
    int qA = 0;
    int qB = 0;
    for(short int qVertex=0;qVertex<4;qVertex++)
    {
  		if(  xc[qCoord] > coords[Tetra.vertex[qVertex]].coord[qCoord] + IN_ELEM_TOLL   )
  		{
    		qA++;
  		}
  		if(  xc[qCoord] < coords[Tetra.vertex[qVertex]].coord[qCoord] - IN_ELEM_TOLL    )
      {
  		  qB++;
      }
    }
  	if(qA > 3 || qB > 3)
  	{
  				notInElement = true;
  	}	
  }
  return(!notInElement);
}




GrahmOperatorOutput ThicknessEvaluation::GrahmOperations( const Point & p0, const Point & p1, const Point & p2, const std::vector<double> & x_c)
{
  std::vector<double> N(3,0);
  std::vector<double> p0p1(3,0), p0p2(3,0);
  std::vector<double> p0pxc(3,0);
  for(short int jCoord=0;jCoord<3;jCoord++)
  {
		p0p1[jCoord] = p0.coord[jCoord] - p1.coord[jCoord];
		p0p2[jCoord] = p0.coord[jCoord] - p2.coord[jCoord];
		p0pxc[jCoord] = x_c[jCoord] - p0.coord[jCoord]; 
  }
  // Calculates normal by cross-product
  N[0] = (p0p1[1]*p0p2[2] - p0p1[2]*p0p2[1]);
  N[1] = (p0p1[2]*p0p2[0] - p0p1[0]*p0p2[2]);
  N[2] = (p0p1[0]*p0p2[1] - p0p1[1]*p0p2[0]);
	
  // Normalises
  double NMag = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
  for(short int jCoord=0;jCoord<3;jCoord++)
  {
      N[jCoord] = N[jCoord]/NMag;
  }
		
  // Calculates constant for plane equation
  double d = -(N[0]*p0.coord[0] + N[1]*p0.coord[1] + N[2]*p0.coord[2]);

  // Checks if the centroid of the element is a positive or negative distance from the plane
  double dist = (N[0]*x_c[0] + N[1]*x_c[1] + N[2]*x_c[2] + d);

  // If positive, then the normal is pointing into the tissue and all is ok
  // If negative, then we need to reverse the direction of the triangle normal
  if(dist < 0.0)
  {
    for(short int jCoord=0;jCoord<3;jCoord++)
    {
      N[jCoord] = -N[jCoord];
    }
    dist=-1.0*dist;
  }
	std::vector<double> projPt(3,0);
	for(short int jCoord=0;jCoord<3;jCoord++)
	{
	  projPt[jCoord] = p0.coord[jCoord]+(p0pxc[jCoord]-dist*N[jCoord]);
	}
	GrahmOperatorOutput result;
	result.N=N;
	result.dist=dist;
	result.projPt=projPt;
	
	return(result);
}





double ThicknessEvaluation::distanceOfPointToPlane(const std::vector<double> & N, const Point & p0,  const std::vector<double> & x_c, const std::vector<double> & mu)
{

	// Plug the intersection point (P') into the plane equation to get the value of a (i.e. distance of P from plane along U)
	// N.(P + Ut) + d = 0
	// => t = -(NP + d)/(N.U)
	double NdotQ = (N[0]*p0.coord[0] + N[1]*p0.coord[1] + N[2]*p0.coord[2]);
	double d = -NdotQ;
	double NdotP = (N[0]*x_c[0] + N[1]*x_c[1] + N[2]*x_c[2]);
	double NdotU = (N[0]*mu[0] + N[1]*mu[1] + N[2]*mu[2]);
	double a = -(NdotP + d)/(NdotU);
	return(a);
}



double ThicknessEvaluation::normalDistanceOfPointToPlane(const std::vector<double> & N, const Point & p0, const std::vector<double> & x_c)
{
  double NdotQ = (N[0]*p0.coord[0] + N[1]*p0.coord[1] + N[2]*p0.coord[2]);
  double d = -NdotQ;
  double NdotP = (N[0]*x_c[0] + N[1]*x_c[1] + N[2]*x_c[2]);
  double a = -(NdotP + d);
  return(a);
}

std::vector<double> ThicknessEvaluation::waxpy(const std::vector<double> & x, const std::vector<double> & y, double  alpha)
{
  
  // implements w=ax+y
  size_t vecsize=x.size();
#ifndef NDEBUG
  if(vecsize != y.size())
  {
    std::cerr<<"ERROR: x and y have different sizes"<<std::endl;
    exit(1);
  }
#endif
  std::vector <double> w=y;
  for(size_t ir=0; ir<vecsize; ir++)
  {
    w[ir]=y[ir]+alpha*x[ir];
  }
  return(w);
}

