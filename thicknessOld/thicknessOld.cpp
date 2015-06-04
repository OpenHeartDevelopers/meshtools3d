/////////////////////////////////////////////////////////////////////////////////
// Monte Carlo Photon Scattering
// Simulates illumination with an unstructured finite element domain
//
// Martin Bishop, Department of Biomedical Engineering, King's College London
// 6th December 2011
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <list>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include "IO.h"
#include "loading.h"
#include "functions.h"

using namespace std;

int main( int argc, char *argv[] )
{
    	
	if (argc < 3)
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << " Mesh Basename " << std::endl;
		std::cerr << " endo surface list" << std::endl;
		std::cerr << " gradient vectors file name from Laplace solve " << std::endl;
		std::cerr << " epi (ground) surface list " << std::endl;
		std::cerr << " flag specifying migration direction (1 == endo-epi; 0 == epi-endo) " << std::endl;
		return 1;
	}
	cout.precision(15);	

	double tol = 1.0;

	//////////////////////////////
	// Mesh
	//////////////////////////////
	cout << "*******************************\n";
  cout << "Reading-in mesh...\n";
  cout << "*******************************\n";
	// Reads Nodes
	// Cesare: perhaps multiplication by 1000.0 is useless (input mesh is still in um)
	double** coords = NULL;
	unsigned int num_nodes;
	loadNodes(argv[1],coords,num_nodes);
	for(unsigned int i=0;i<num_nodes;i++)
		for(int j=0;j<3;j++)
			coords[i][j] = coords[i][j]*1000.0;

	// Reads Elements
	int** elems = NULL;
	unsigned int num_elems;
	loadElements(argv[1],elems,num_elems);

	// Does a check on element tags to switch to defaults (i.e. bath=0, tissue=2)			
	// Cesare: NOT CLEAR FOR ME
	for(unsigned int n=0;n<num_elems;n++)
		if(elems[n][4] == 1)
			elems[n][4] = 2;
	else if(elems[n][4] == 3 || elems[n][4] == 4)
		elems[n][4] = 0;

	// Reads Surface triangles (all surf boundary)
	int** tris = NULL;
	unsigned int num_tris;
	loadTris(argv[1],tris,num_tris);

	// Illumination surface - list of nodes which are illuminated
	// Cesare: Reads ENDOocardium
	int* illumSurf = NULL;
	unsigned int num_surfNodes;
	loadSurfaceList(argv[2],illumSurf,num_surfNodes);
	// Makes a set of illumination surface nodes
	set<int> illumNodesSet ;
	set<int>::iterator it_surf;
	for(unsigned int i=0;i<num_surfNodes;i++)
		illumNodesSet.insert(illumSurf[i]);

	// Ground surface 
	// Cesare: Reads EPIcardium
  int* groundSurf = NULL;
  unsigned int num_surfNodesGround;
  loadSurfaceList(argv[4],groundSurf,num_surfNodesGround);
	// Makes a set of illumination surface nodes
  set<int> groundNodesSet ;
  for(unsigned int i=0;i<num_surfNodesGround;i++)
    groundNodesSet.insert(groundSurf[i]);

	// Element centroid file
	// Cesare: Reads Centroids
	// Cesare: perhaps multiplication by 1000.0 is useless (input mesh is still in um)
	double** cents = NULL;
	loadCentroids(argv[1],cents,num_elems);
	for(unsigned int i=0;i<num_elems;i++)
                for(int j=0;j<3;j++)
                        cents[i][j] = cents[i][j]*1000.0;

	// Reads-in vectors from Laplace solve
	// Cesare: Reads Gradients
	double** grads = NULL;
	loadVectors(argv[3],grads,num_elems);
	int flag = atoi(argv[5]);
	if(flag == 1)
	{
		for(unsigned int i=0;i<num_elems;i++)
			for(int h=0;h<3;h++)
				grads[i][h] *= -1.0;
	}

	//Cesare: From now on algorithms starts
	
	//////////////////////////////
	// Coordinate Variables
	//////////////////////////////
	// Coordinates
	double x, y, z; 

	// Surface coordinates
	double x_epi, y_epi, z_epi;

	// Intermediate steps
	double x_s, y_s, z_s;

	// Angles of travel
	double mu_x, mu_y, mu_z;
	double mu_x_new, mu_y_new, mu_z_new;	
	mu_x_new = 0.0;
	mu_y_new = 0.0;
	mu_z_new = 0.0;
	//////////////////////////////
	// Storing Vectors
	//////////////////////////////
	// pathLength
	double* path;
	path = new double [num_tris];
	for(unsigned int i=0;i<num_tris;i++)
		path[i] = 0.0;

	double* illumElems;
	illumElems = new double [num_elems];
	for(unsigned int i=0;i<num_elems;i++)
		illumElems[i] = 0.0;

	//////////////////////////////
	// PRE-PROCESSING: Creates connectivity array which lists all elements attached to each node in mesh
	//////////////////////////////
	cout << "*******************************\n";
  cout << "Performing preprocessing...\n";
  cout << "*******************************\n";

	cout << " Producing node-element connectivity array..." << "\n";
	// Defines array
	int** connNodes = NULL;
	produceConnectivityArray(connNodes,elems,num_nodes,num_elems);

	cout << " Producing node-triangle connectivity array..." << "\n";
	// Defines array
	int** connNodesTris = NULL;
	produceConnectivityArrayTris(connNodesTris,tris,num_nodes,num_tris);


	//////////////////////////////
	// PRE-PROCESSING: Creates a list of which element faces connect which element faces 
	//////////////////////////////
	cout << "Constructing element face-to-face look-up table... \n";
	// Defines array
	int** faceToFace = NULL;
	faceToFace = new int*[num_elems];
	// Initialises to -1
	for(unsigned int i=0;i<num_elems;i++)
	{
		faceToFace[i] = new int[4];
		for(int j=0;j<4;j++)
			faceToFace[i][j] = -1;
	}

	for(unsigned int i=0;i<num_elems;i++)
	{
		// Makes a set listing all elements directly connected to this element
		set<int> elemConnElem;
		set<int>::iterator it_ce;
		elemConnElem.clear();
		// Iterates over each node in current element
		for(int j=0;j<4;j++)
		{
			int node = elems[i][j];
			// Iterates-over all possible entries in connNodes list
			for(int k=0;k<100;k++)
			{
				if(connNodes[node][k] != -1)
					elemConnElem.insert(connNodes[node][k]);
				else
					break;
			}
		}

		// Removes the current element from the list!
		elemConnElem.erase(i);

		// Creates a set with all nodes defining element
		set<int>elemNodes;
		set<int>::iterator i2;
		elemNodes.clear();
		for(int n=0;n<4;n++)
			elemNodes.insert(elems[i][n]);

		// Iterates over each triangle face
		for(int n=0;n<4;n++)
		{
			// Erases one node (which then defines a triangle)
			elemNodes.erase(elems[i][n]);

			// Iterates-over all elements in list of connected elements
			for(it_ce=elemConnElem.begin();it_ce!=elemConnElem.end();it_ce++)
			{
				int connElem = *it_ce;

				// Tries to find all 3 current triangle nodes (in turn) within this connected element
				int p = 0;
				for(int g=0;g<4;g++)
				{
					if(elemNodes.find(elems[connElem][g]) != elemNodes.end())
						p++;
				}
				// If p = 3 the we've found this triangle and thus this element must connect via the current element face
				if(p == 3)
				{
					faceToFace[i][n] = connElem;
					break;
				}
			}

			// Puts the erased node back before we try to move-on to the next triangle
			elemNodes.insert(elems[i][n]);
		}
	}
	//cout << "faceToFace = " << faceToFace[100][2] << "\n";

	//////////////////////////////
	// PRE-PROCESSING: Create look-up table to list which elements have boundary triangles (connected via FACES)
	//////////////////////////////
	cout << "Creating look-up table to list elements which have boundary triangles connected via FACES ... \n";
	// Array to store boundary triangle numbers
	int** elemBoundaryTris;
	elemBoundaryTris = new int*[num_elems];
	for(unsigned int i=0;i<num_elems;i++)
		elemBoundaryTris[i] = new int[4];

	// Initialises all entries to -1
	for(unsigned int i=0;i<num_elems;i++)
		for(int j=0;j<4;j++)
			elemBoundaryTris[i][j] = -1;


	for(unsigned int i=0;i<num_elems;i++)
	{
		// Makes a set listing all triangles directly connected to this element
		set<int> elemConnTris;
		set<int>::iterator it_ct;
		elemConnTris.clear();
		// Iterates over each node in current element
		for(int j=0;j<4;j++)
		{
			int node = elems[i][j];
			// Iterates-over all possible entries in connNodesTris list
			for(int k=0;k<50;k++)
			{
				// Adds connected triangle to list
				if(connNodesTris[node][k] != -1)
					elemConnTris.insert(connNodesTris[node][k]);
				else
					break;
			}
		}

		// Makes a set of all 4 element nodes
		set<int>elemNodes;
		set<int>::iterator it2;
		elemNodes.clear();
		for(int n=0;n<4;n++)
			elemNodes.insert(elems[i][n]);

		// Iterates over each triangle face
		for(int n=0;n<4;n++)
		{
			// Erases one node (which then defines a triangle)
			elemNodes.erase(elems[i][n]);

			// Iterates-over all triangles in list of connected elements
			for(it_ct=elemConnTris.begin();it_ct!=elemConnTris.end();it_ct++)
			{
				int connTri = *it_ct;
				int c = 0;
				// Tries to find each triangle node within the list of current element nodes
				for(int m=0;m<3;m++)
				{
					it2 = elemNodes.find(tris[connTri][m]);
					if(it2 != elemNodes.end())
						c++;
				}
				// If it finds more than 2 nodes (i.e. 3) then it inserts the current triangle number into the correct face number for this current element
				if(c > 2)
				{
					elemBoundaryTris[i][n] = connTri;
					break;
				}
			}
			// Puts the node we erased back before we move-on to the next triangle
			elemNodes.insert(elems[i][n]);
		}
	}

	int counterTris = 0;
	for(unsigned int q=0;q<num_elems;q++)
	{
		for(int q2=0;q2<4;q2++)
			if(elemBoundaryTris[q][q2] == -1)
				counterTris++;
	}
	cout << "Number of bounding Triangles = " << num_elems*4 - counterTris << "\n";

	//////////////////////////////
	// PRE-PROCESSING: Derives surface illumination triangles
	//////////////////////////////
	cout << "Making list of surface triangles to start migration... \n";
	// Defines a list of surface triangles which contain illuminating nodes
	set<int> illumTris;
	set<int>::iterator it;

	// Finds which surface triangle contains a surface illumination node and inserts the tri number into the set
	/*
	int surfNode;
	for(unsigned int i=0;i<num_surfNodes;i++)
	{
		surfNode = illumSurf[i];

		for(unsigned int j=0;j<num_tris;j++)
		{
			if(surfNode == tris[j][0] || surfNode == tris[j][1] || surfNode == tris[j][2])
				illumTris.insert(j);
		}	
	}
	*/

	for(unsigned int j=0;j<num_tris;j++)
        {
		for(int k=0;k<3;k++)
		{
			// Tries to find any triangle node within list of surface nodes
			it_surf = illumNodesSet.find(tris[j][k]);
			// If it does, add triangle to list of illuminated triangles
			if(it_surf != illumNodesSet.end() )
				illumTris.insert(j);	
		}
	}

	// Finds size (or total number of illuminating surface triangles)
	int numIllumTris = illumTris.size();
	cout << "Number of triangles which we are starting migration from  = " << numIllumTris << "\n";

	//////////////////////////////
        // PRE-PROCESSING: Defines list of surface triangles corresponding to other electrode
        //////////////////////////////
        
	// Defines a list of surface triangles which contain illuminating nodes
        set<int> groundTris;
        set<int>::iterator itg;

	/*
        for(unsigned int i=0;i<num_surfNodesGround;i++)
        {
                surfNode = groundSurf[i];

                for(unsigned int j=0;j<num_tris;j++)
                {
                        if(surfNode == tris[j][0] || surfNode == tris[j][1] || surfNode == tris[j][2])
                                groundTris.insert(j);
                }
        }
	*/
	for(unsigned int j=0;j<num_tris;j++)
        {
                for(int k=0;k<3;k++)
                {
                        // Tries to find any triangle node within list of surface nodes
                        it_surf = groundNodesSet.find(tris[j][k]);
                        // If it does, add triangle to list of illuminated triangles
                        if(it_surf != groundNodesSet.end() )
                                groundTris.insert(j);
                }
        }

        // Finds size (or total number of illuminating surface triangles)
        int numGroundTris = groundTris.size();
        cout << "Number of triangles which form part of ground = " << numGroundTris << "\n";

	//////////////////////////////
	// PRE-PROCESSING: Associating each surface illumination triangle with corresponding tissue element
	//////////////////////////////
	cout << "Making a look-up table to associate surface triangles from where migration starts with corresponding tissue elements... \n";

	int **illumElemsWithTris;
	illumElemsWithTris = new int*[numIllumTris];
	for(int i=0;i<numIllumTris;i++)
		illumElemsWithTris[i] = new int[2];

	int c8 = 0;
	for(it=illumTris.begin();it!=illumTris.end();it++)
	{
		// Picks-out particular triangle
		int thisIllumTri = *it;

		// Makes a set listing all elements directly connected to this triangle
		set<int> triConnElem;
		set<int>::iterator it_ct2;
		triConnElem.clear();
		// Iterates over each node in current tri
		for(int j=0;j<3;j++)
		{
			int node = tris[thisIllumTri][j];
			// Iterates-over all possible entries in connNodes list
			for(int k=0;k<100;k++)
			{
				if(connNodes[node][k] != -1)
					triConnElem.insert(connNodes[node][k]);
				else
					break;
			}
		}

		// Iterates-over all elements in list of connected elements
		for(it_ct2=triConnElem.begin();it_ct2!=triConnElem.end();it_ct2++)
		{
			int connElem = *it_ct2;
			if(elems[connElem][4] == 2)
			{
				// Makes a set of all 4 element nodes
				set<int> elemNodes;
				elemNodes.clear();
				for(int j=0;j<4;j++)
					elemNodes.insert(elems[connElem][j]);

				// Removes (or tries to remove) all 3 nodes defined by triangle 
				for(int j=0;j<3;j++)
					elemNodes.erase(tris[thisIllumTri][j]);

				// If we've successfully removed 3 nodes, then this trianlge is a face of this element
				if(elemNodes.size() < 2)
				{
					illumElemsWithTris[c8][1] = connElem;
					illumElemsWithTris[c8][0] = thisIllumTri;
				}
			}
		}
		c8++;
	}	
	//cout << " illumElemsWithTris = " << illumElemsWithTris[100][0]  << " " << illumElemsWithTris[100][1] << "\n";
	// Checks that faceToFace entries match elemBoundaryTris entires
	int badElement = 0;
	for(unsigned int i=0;i<num_elems;i++)
	{
		for(int h=0;h<4;h++)
		{
			if(faceToFace[i][h] == -1 && elemBoundaryTris[i][h] == -1)
				badElement++;
		}
	}
	cout << "Number of bad elements = " << badElement << "\n";

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	// MAIN - ILLUMINATES TISSUE WITH PHOTONS
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	cout << "*******************************\n";
        cout << "Computing wall thickness...\n";
        cout << "*******************************\n";
	// Iterates over all surface illumination triangles
	// Sets path going from each surface illuminating triangle
	int interactions=0;
	int problems=0;
	for(it=illumTris.begin();it!=illumTris.end();it++)
	{
		// Picks-out particular triangle
		int thisIllumTri = *it;
		// Picks-out the TISSUE element associated with this surface triangle from which we are illuminating
		int surfaceElem = 0;
		for(int j=0;j<numIllumTris;j++)
			if(illumElemsWithTris[j][0] == thisIllumTri)
		{
			surfaceElem = illumElemsWithTris[j][1];
			illumElems[surfaceElem] = 1.0;
		}

		//////////////////////////////
		// Iterates over all nodes of selected triangle and derives centroid (of triangle)
		//////////////////////////////
		x_epi = 0.0, y_epi = 0.0, z_epi = 0.0;
		int node;
		for(int j=0;j<3;j++)
		{
			// One of the 3 triangle nodes
			node = tris[thisIllumTri][j];
			// Calculates centroid of triangle
			x_epi += coords[node][0]/3.0;
			y_epi += coords[node][1]/3.0;
			z_epi += coords[node][2]/3.0;	
		}

		// Set-ups list to store the 4 nodes in each element        	
		set<int>elemNodes;
		set<int>::iterator i2;

		// Sets-up list to store the previous triangle nodes from th face we just passed through
		set<int> faceRem;
		set<int>::iterator it_fr;
		// Initialises this list with the nodes of the epicardial triangle surface face
		faceRem.clear();
		for(int n=0;n<3;n++)
			faceRem.insert(tris[thisIllumTri][n]);

		// Initial coordinates of photon packet set as centroid of triangle
		x = x_epi;
		y = y_epi;
		z = z_epi;

		// Defines initial element we're in
		int element = surfaceElem;
		int oldElement = 0;

		// Define initial direction as vector gradient direction
		mu_x = grads[element][0];
		mu_y = grads[element][1];
		mu_z = grads[element][2];

		// Defines some counters and flags
		int inTissue = 1;
		int newCounter = 0;
		int notInElement = 0;
		int adjustMu = 0;
		int wrongSurface = 0;
		double lastaMin = 0.0;
		int lastintface = 0;
		
		while(inTissue == 1)
		{
				//cout << "x,y,z = " << x << " " << y << " " << z << "\n";
				//cout << "element = " << element << "\n";
				//cout << "mus " << mu_x << " " << mu_y << " " << mu_z << "\n";
				// Counts total interactions
				interactions++;				
				
				// Checks path to see if it's too long
				if(path[thisIllumTri] > 50000)
					break;
					

				///////////////////////////////
				// Quick check to see if we are not in the element we should be in
				///////////////////////////////                 
				double * r;
				r = new double[3];
				r[0] = x;
				r[1] = y;
				r[2] = z;
				notInElement = 0;
				for(int q=0;q<3;q++)
				{
					int qA = 0;
					int qB = 0;
					for(int q2=0;q2<4;q2++)
					{
						if(r[q] > coords[elems[element][q2]][q] + tol)
							qA++;
						if(r[q] < coords[elems[element][q2]][q] - tol)
							qB++;
					}
					if(qA > 3 || qB > 3)
					{
						notInElement = 1;
					}	
				}

				newCounter++;
				if(newCounter > 1000)
					inTissue = 0;

				if(adjustMu == 0 && wrongSurface == 0)
				{
					// Initial directional cosines set from normal of triangle (calculating dot-product of each axis with normal)
					mu_x = grads[element][0];
					mu_y = grads[element][1];
					mu_z = grads[element][2];
				}
				else 
				{
					mu_x = mu_x_new;
					mu_y = mu_y_new;
					mu_z = mu_z_new; 
				}
				
				// Updates starting position of step
				x_s = x;
				y_s = y;
				z_s = z;
				
				wrongSurface = 0;			
				
				//////////////////////////////
				// Check for intersection with faces of current element
				//////////////////////////////
				// Defines things to define closest face
				double aMin = 1000000.0;

				int intFace = -1;

				// Creates a set with all nodes defining element
				elemNodes.clear();
				for(int n=0;n<4;n++)
					elemNodes.insert(elems[element][n]);

				double Nmin[3];
				for(int h=0;h<3;h++)
					Nmin[h] = 0.0;
				// Iterates over each triangle face
				for(int n=0;n<4;n++)
				{
					// Erases one node (which then defines a triangle)
					elemNodes.erase(elems[element][n]);
					// Defines each of 3 remaining nodes which now define the triangle face
					int triNodes[3];
					int p = 0;
					for(i2=elemNodes.begin();i2!=elemNodes.end();i2++)
					{
						triNodes[p] = *i2;
						p++;
					}

					// Checks to see if this is the same face and we've just pass through (if it is, don't include it)
					int c = 0;
					for(int m=0;m<3;m++)
					{
						it_fr = faceRem.find(triNodes[m]);
						if(it_fr != faceRem.end())
							c++;
					}
					// Only if this isn't the face we've just interacted with do we now calculate an a value
					if(c < 3)
					{
						// Calculates inward-facing normal of this triangle
						double* N = NULL;
						N = calculateInwardNormal(element,triNodes[0],triNodes[1],triNodes[2],coords,cents);
						// Calculates intersection distance
						double a = distanceOfPointToPlane(N,triNodes[0],coords,x_s,y_s,z_s,mu_x,mu_y,mu_z);
						if(a >= 0 && a<aMin)
						{
							// Updates smallest intersection distance
							aMin = a;
							// Defines closest intersecting face
							intFace = n;
							for(int h=0;h<3;h++)
								Nmin[h] = N[h];
						}
						// Calculates NORMAL distance of point to plane
						double a_n = normalDistanceOfPointToPlane(N,triNodes[0],coords,x_s,y_s,z_s);
						//double dot = N[0]*grads[element][0] + N[1]*grads[element][1] + N[2]*grads[element][2];
						if(a_n > 0)
							problems++;

						while(a_n > 0)
						{
							
							// Move along line towards centroid by a small amount and recompute a
							double delta = 1.0;
							double v[3];
							v[0] = cents[element][0] - x;
							v[1] = cents[element][1] - y;
							v[2] = cents[element][2] - z;
							double magV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
							for(int h=0;h<3;h++)
								v[h] = v[h]/magV;
							
							x = x + delta*v[0];
							y = y + delta*v[1];
							z = z + delta*v[2];
							a_n = normalDistanceOfPointToPlane(N,triNodes[0],coords,x,y,z);
							a = distanceOfPointToPlane(N,triNodes[0],coords,x,y,z,mu_x,mu_y,mu_z);
							if(a >= 0 && a<aMin)
                                                	{
                                                        	// Updates smallest intersection distance
                                                        	aMin = a;
                                                        	// Defines closest intersecting face
                                                        	intFace = n;
                                                        	for(int h=0;h<3;h++)
                                                                	Nmin[h] = N[h];
                                                	}
						}

						delete[] N;
					}

					// Puts the erased node back
					elemNodes.insert(elems[element][n]);
				}

				if(intFace != -1)
				{
					lastintface = intFace;
					lastaMin = aMin;
					// Add details of face we've just moved to
					elemNodes.erase(elems[element][intFace]);
					// Defines each of 3 remaining nodes which now define the triangle face
					faceRem.clear();
					for(i2=elemNodes.begin();i2!=elemNodes.end();i2++)
						faceRem.insert(*i2);


					path[thisIllumTri] += aMin;

					// Checks to see if this closest intersecting face is also a boundary face
					// This is for interior boundaries
					if(elemBoundaryTris[element][intFace] != -1)
                                        {
                                                // Checks to see if this is part of the ground electrode
                                                if(groundTris.find(elemBoundaryTris[element][intFace]) != groundTris.end())
                                                {
                                                        inTissue = 0;
                                                        break;
                                                }
                                                else
                                                {
							wrongSurface = 1;
							// Move a small distance in towards the centroid
							double delta = 1.0;
                                                        double v[3];
                                                        v[0] = cents[element][0] - x;
                                                        v[1] = cents[element][1] - y;
                                                        v[2] = cents[element][2] - z;
                                                        double magV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
                                                        for(int h=0;h<3;h++)
                                                                v[h] = v[h]/magV;
                                                        x = x + delta*v[0];
                                                        y = y + delta*v[1];
                                                        z = z + delta*v[2];

							// Calculate the normal to the current gradient and the face normal
							double w[3],u[3];
							for(int h=0;h<3;h++)
								u[h] = grads[element][h];
							double N1[3];
							N1[0] = (u[1]*Nmin[2] - u[2]*Nmin[1]);
						        N1[1] = (u[2]*Nmin[0] - u[0]*Nmin[2]);
						        N1[2] = (u[0]*Nmin[1] - u[1]*Nmin[0]);

							// Calculate the normal to the face normal and this normal to get the vector parallel to the face
							w[0] = (N1[1]*Nmin[2] - N1[2]*Nmin[1]);
							w[1] = (N1[2]*Nmin[0] - N1[0]*Nmin[2]);
							w[2] = (N1[0]*Nmin[1] - N1[1]*Nmin[0]);
							double magW = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
							for(int h=0;h<3;h++)
								w[h] = w[h]/magW;

							// Check using the dot product that this new vector is in the same direction as the previous gradient
							double dot3 = w[0]*u[0] + w[1]*u[1] + w[2]*u[2];
							if(dot3 < 0)
								for(int h=0;h<3;h++)
                                                                	w[h] = -1.0*w[h];

							mu_x_new = w[0];
                                                        mu_y_new = w[1];
                                                        mu_z_new = w[2];
                                                }
                                        }

					// We also update our 'current element' as the one we're moving into as we cross this boundary
					int newElement = faceToFace[element][intFace];

					// Only if we haven't passed through a boundary
					if(newElement != -1)
					{
						oldElement = element;
						element = newElement;

						mu_x_new = grads[newElement][0];
						mu_y_new = grads[newElement][1];
						mu_z_new = grads[newElement][2];
						if(Nmin[0]*mu_x_new + Nmin[1]*mu_y_new + Nmin[2]*mu_z_new >= 0)
						{	
							adjustMu = 1;
							mu_x_new = 0.5*grads[newElement][0] + 0.5*grads[oldElement][0];
							mu_y_new = 0.5*grads[newElement][1] + 0.5*grads[oldElement][1];
							mu_z_new = 0.5*grads[newElement][2] + 0.5*grads[oldElement][2];
							if(Nmin[0]*mu_x_new + Nmin[1]*mu_y_new + Nmin[2]*mu_z_new >= 0)
							{
								element = oldElement;
							}
						}
						else
							adjustMu = 0;

						// Move-up to face of element and update position
                                                x = x_s + mu_x*aMin;
                                                y = y_s + mu_y*aMin;
                                                z = z_s + mu_z*aMin;
						
						// Instead of this, move a tiny amount into the new element
						double v[3];
						v[0] = cents[newElement][0] - x;
						v[1] = cents[newElement][1] - y;
						v[2] = cents[newElement][2] - z;
						double vMag = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2]  );
						double delta = 0.0;
						x = x + delta*v[0]/vMag;
						y = y + delta*v[1]/vMag;
						z = z + delta*v[2]/vMag;
					}
				}
		}
	}

cout << "Total number of interactions = " << interactions << "\n";
double propProb = double(problems)/double(interactions);
cout << "Total number of problems = " << problems << " ( or " << propProb << " ) \n";
cout << "WRITING-OUT DATA FILES... \n";


// Converts data for all illum tris onto a full element list basis
double* elemData;
elemData = new double[num_elems];
for(unsigned int i=0;i<num_elems;i++)
	elemData[i] = 0.0;

int elem_n = 0;
int tri_n = 0;
for(int i=0;i<numIllumTris;i++)
{
	elem_n = illumElemsWithTris[i][1];
	tri_n = illumElemsWithTris[i][0];
	elemData[elem_n] = path[tri_n];
}

// Maps data over to node list keeping track of how many surface triangles are associated with each node 
double* surfPlotter;
surfPlotter = new double[num_nodes];
for(unsigned int i=0;i<num_nodes;i++)
	surfPlotter[i] = 0.0;

int* surfPlotterCounter;
surfPlotterCounter = new int[num_nodes];
for(unsigned int i=0;i<num_nodes;i++)
	surfPlotterCounter[i] = 0;

int node_n;
for(int i=0;i<numIllumTris;i++)
{
	tri_n = illumElemsWithTris[i][0];
	for(int j=0;j<3;j++)
	{
		node_n = tris[tri_n][j];
		surfPlotter[node_n] += path[tri_n];
		surfPlotterCounter[node_n]++;
	}
}

for(unsigned int i=0;i<num_nodes;i++)
{
	if(surfPlotterCounter[i] != 0)
		surfPlotter[i] = surfPlotter[i]/double(surfPlotterCounter[i]);
}

//string output = "output.dat";
//writeDataFile(output,num_elems,elemData);

//string outputElem = "outputElem.dat";
//writeDataFile(outputElem,num_elems,illumElems);

string outputNodes = "outputNodes.dat";
writeDataFile(outputNodes,num_nodes,surfPlotter);

delete[] tris;
delete[] elems;
delete[] cents;
delete[] faceToFace;
delete[] illumElemsWithTris;
delete[] illumSurf;
delete[] connNodes;
delete[] elemBoundaryTris;		
delete[] path;
delete[] elemData;
delete[] surfPlotterCounter;
delete[] surfPlotter;

return 0;

}














