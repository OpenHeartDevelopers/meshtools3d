#ifndef _LOADING_H_
#define _LOADING_H_

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// Function to load-in nodes
void loadNodes(const char* file, double** &coords, unsigned int &num_nodes);

// Function to load-in elements
void loadElements(const char* file, int** &elems, unsigned int &num_elems);

// Function to load-in triangles file
void loadTris(const char* file, int** &tris, unsigned int &num_tris);

// Function to load-in list of surface nodes
void loadSurfaceList(const char* file, int* &surface, unsigned int &num_surfnodes);

// Function to load-in nodes
void loadCentroids(const char* file, double** &cents, unsigned int &num_elems);

// Function to load-in vector file
void loadVectors(const char* file, double** &vectors, unsigned int num_elems);

// Function to define a connectivity array of which elements are connected to each node
void produceConnectivityArray(int** &conn_nodes,int** elems,unsigned int &num_nodes,unsigned int &num_elems);

// Function to define a connectivity array of which triangles are connected to each node
void produceConnectivityArrayTris(int** &conn_nodes,int** tris,unsigned int &num_nodes,unsigned int &num_tris);
///
// Function to load-in trianges and individual list of surface nodes
void loadTrisWithNtris(const char* file, int** &tris, unsigned int &num_tris, int* &n_tris, unsigned int &num_n_tris);

// Function to load-in triangles file
void loadTrisWithNoDummy(const char* file, int** &tris, unsigned int &num_tris);

// Function to load-in triangles file
void loadTrisWithTags(const char* file, int** &tris, unsigned int &num_tris);

// Function to load-in fibre orientation (lon) file
void loadFibres(const char* file, double** &fibres, unsigned int &num_elems);

// Function to load-in fibre orientation (lon) file
void loadFibresNew(const char* file, double** &fibres, unsigned int &num_elems);

// Function to load-in vector file
void loadVectors4(const char* file, double** &vectors, unsigned int num_elems);


#endif
