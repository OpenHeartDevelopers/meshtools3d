#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <set>
#include <list>

using namespace std;

// Function to load-in nodes
void loadNodes(const char* file, double** &coords, unsigned int &num_nodes)
{
	// Puts the .pts on end of mesh base-name and loads
	string pts_end = ".pts";
	string total_pts_file = file + pts_end;
	ifstream coords_file(total_pts_file.c_str());

	// Reads-in contents of file
	if (coords_file.is_open())
	{
		cout << "Reading " << total_pts_file << " " << endl << flush;

		// Reads-in line 1 as number of nodes
		coords_file >> num_nodes;

		cout << "Number of node = " << num_nodes << endl << flush;

		// Defines size of array to accommodate nodes
		coords = new double*[num_nodes];

		// Iterates over total number of nodes
		for(unsigned i=0;i<num_nodes;i++)
		{
			// Defines each row to have 3 columns and reads-in to array
			coords[i] = new double[3];
			coords_file >> coords[i][0] >> coords[i][1] >> coords[i][2];
		}
		coords_file.close();
		cout << "Done." << endl << flush;
	}
	else
	{
		cerr << "Can't Open " << file << "!" << endl << flush;
	}
}

// Function to load-in elements
void loadElements(const char* file, int** &elems, unsigned int &num_elems)
{
	// Puts the .elem on end of mesh base-name and loads
	string elem_end = ".elem";
	string total_elem_file = file + elem_end;
	ifstream elems_file(total_elem_file.c_str());

	// Reads-in contents of file
	if (elems_file.is_open())
	{
		cout << "Reading " << total_elem_file << " " << flush;

		// Reads-in line 1 as number of elements
		elems_file >> num_elems;

		// Outputs
		cout << "Number of elements = " << num_elems << endl << flush;

		// Defines size of array to accommodate elements
		elems = new int*[num_elems];

		// Defines a string to read the Tt flag to
		string Tt;

		// Interates over total number of elements
		for(unsigned i=0;i<num_elems;i++)
		{
			// Defines each row to have 5 columns and reads-in to array (including tag as 5th column)
			elems[i] = new int[5];
			elems_file >> Tt >> elems[i][0] >> elems[i][1] >> elems[i][2] >> elems[i][3] >> elems[i][4];
		}
		elems_file.close();

		cout << "Done." << endl << flush;
	}
	else
	{
		cerr << "Can't Open " << file << "!" << endl << flush;
	}
}

// Function to load-in list of triangles
void loadTris(const char* file, int** &tris, unsigned int &num_tris)
{
	// Puts the .tris on end of mesh base-name and loads
	string tris_end = ".tris";
	string total_tris_file = file + tris_end;
	ifstream tris_file(total_tris_file.c_str());

	// Reads-in file
	if (tris_file.is_open())
	{
		cout << "Reading " << total_tris_file << " " << endl << flush;

		// Reads-in line 1 as total number of triangles
		tris_file >> num_tris;
	
		//string dummy;
		
		// Outputs
		cout << "Number of Trianlges = " << num_tris << endl << flush;

		// Defines array to accommodate triangles
		tris = new int*[num_tris];

		// Defines tag (if triangles have it)
		//int tag;

		// Iterates over total number of triangles
		for(unsigned i=0;i<num_tris;i++)
		{
			// Defines array to have 3 rows
			tris[i] = new int[3];
			// Inputs triangles to array (but not tag)
     			//tris_file >> dummy >> tris[i][0] >> tris[i][1] >> tris[i][2] >> tag;
			//tris_file >> dummy >> tris[i][0] >> tris[i][1] >> tris[i][2];
			tris_file >> tris[i][0] >> tris[i][1] >> tris[i][2];
		}
		tris_file.close();
		cout << "Done." << endl << flush;
	}

	else
	{
		cerr << "Can't Open " << file << "!" << endl << flush;
	}
}

// Function to load-in list of (individual) surface nodes
void loadSurfaceList(const char* file, int* &surface, unsigned int &num_surfnodes)
{
	// Defines filename
	ifstream surf_file(file);

	// Reads-in contents of file
	if (surf_file.is_open())
	{
		cout << "Reading " << file << " "  << endl << flush;

		// Reads-in line 1 as number of surface nodes
		surf_file >> num_surfnodes;

		// Outputs
		cout << "Number of surface nodes = " << num_surfnodes  << endl << flush;

		// Defines size of array to accommodate surface nodes
		surface = new int [num_surfnodes];

		// Iterates over total number of surface nodes
		for(unsigned i=0;i<num_surfnodes;i++)
		{	
			// Reads-in to array
			surf_file >> surface[i];
		}
		surf_file.close();
		cout << "Done." << endl << flush;
	}
	else
	{
		cerr << "Can't Open " << file << "!" << endl << flush;
	}
}

// Function to load-in centroids
void loadCentroids(const char* file, double** &cents, unsigned int &num_elems)
{
	// Puts the .pts on end of mesh base-name and loads
	string cpts_end = ".cpts";
	string total_cpts_file = file + cpts_end;
	ifstream cents_file(total_cpts_file.c_str());

	// Reads-in contents of file
	if (cents_file.is_open())
	{
		cout << "Reading " << total_cpts_file << " " << endl << flush;

		// Reads-in line 1 as number of nodes
		cents_file >> num_elems;

		cout << "Number of centroids = " << num_elems << endl << flush;

		// Defines size of array to accommodate nodes
		cents = new double*[num_elems];

		// Iterates over total number of nodes
		for(unsigned i=0;i<num_elems;i++)
		{
			// Defines each row to have 3 columns and reads-in to array
			cents[i] = new double[3];
			cents_file >> cents[i][0] >> cents[i][1] >> cents[i][2];
		}
		cents_file.close();
		cout << "Done." << endl << flush;
	}
	else
	{
		cerr << "Can't Open " << file << "!" << endl << flush;
	}
}


void loadVectors(const char* file, double** &vectors, unsigned int num_elems)
{
	ifstream vec_file(file);

	if (vec_file.is_open())
	{
		cout << "Reading " << file << " " << endl << flush;

		vectors = new double*[num_elems];

		for(unsigned i=0;i<num_elems;i++)
		{
			vectors[i] = new double[3];
			vec_file >> vectors[i][0] >> vectors[i][1] >> vectors[i][2];
		}
		vec_file.close();
		cout << "Done." << endl << flush;
	}
	else
	{
		cerr << "Can't Open " << file << "!" << endl << flush;
	}

}

void produceConnectivityArray(int** &conn_nodes,int** elems,unsigned int &num_nodes,unsigned int &num_elems)
{

	conn_nodes = new int*[num_nodes];
	for(unsigned int i=0;i<num_nodes;i++)
	{
	  conn_nodes[i] = new int[100];
    // Initialises all entries in array to -1
	  for(int j=0;j<100;j++)
	  {
	    conn_nodes[i][j] = -1;
	  }
	}

// Loops-over all nodes of all elements
	for(unsigned int i=0;i<num_elems;i++)
	{
		for(int j=0;j<4;j++)
		{
		// Picks-out node
			int node_n = elems[i][j];
			int k = 0;

		// Finds entry in connectivity array which is not yet filled (i.e. = -1)
			while(conn_nodes[node_n][k] != -1)
			{
				k++;
			}
		// Adds element number to connectivity array
			conn_nodes[node_n][k] = i;
		}
	}
}

void produceConnectivityArrayTris(int** &conn_nodesTris,int** tris,unsigned int &num_nodes,unsigned int &num_tris)
{

  conn_nodesTris = new int*[num_nodes];
  for(unsigned int i=0;i<num_nodes;i++)
  {
    conn_nodesTris[i] = new int[50];
    // Initialises all entries in array to -1
    for(int j=0;j<50;j++)
    {
      conn_nodesTris[i][j] = -1;
    }
  }

  // Loops-over all nodes of all elements
  for(unsigned int i=0;i<num_tris;i++)
  {
      for(int j=0;j<3;j++)
      {
        // Picks-out node
        int node_n = tris[i][j];
        int k = 0;
        // Finds entry in connectivity array which is not yet filled (i.e. = -1)
        while(conn_nodesTris[node_n][k] != -1)
        {
          k++;
        }
        // Adds element number to connectivity array
        conn_nodesTris[node_n][k] = i;
      }
  }
}









// Function to load-in list of triangles
void loadTrisWithNoDummy(const char* file, int** &tris, unsigned int &num_tris)
{
        // Puts the .tris on end of mesh base-name and loads
        string tris_end = ".tris";
        string total_tris_file = file + tris_end;
        ifstream tris_file(total_tris_file.c_str());

        // Reads-in file
        if (tris_file.is_open())
        {
                cout << "Reading " << total_tris_file << " " << endl << flush;

                // Reads-in line 1 as total number of triangles
                tris_file >> num_tris;

                // Outputs
                cout << "Number of Trianlges = " << num_tris << endl << flush;

                // Defines array to accommodate triangles
                tris = new int*[num_tris];

                // Defines tag (if triangles have it)
		int tag;

                // Iterates over total number of triangles
                for(unsigned i=0;i<num_tris;i++)
                {
                        // Defines array to have 3 rows
                        tris[i] = new int[3];
                        // Inputs triangles to array (but not tag)
                        tris_file >> tris[i][0] >> tris[i][1] >> tris[i][2] >> tag;
                }
                tris_file.close();
                cout << "Done." << endl << flush;
        }

        else
        {
                cerr << "Can't Open " << file << "!" << endl << flush;
        }
}



// Function to load-in list of triangles
void loadTrisWithTags(const char* file, int** &tris, unsigned int &num_tris)
{
        // Puts the .tris on end of mesh base-name and loads
        string tris_end = ".tris";
        string total_tris_file = file + tris_end;
        ifstream tris_file(total_tris_file.c_str());

        // Reads-in file
        if (tris_file.is_open())
        {
                cout << "Reading " << total_tris_file << " " << endl << flush;

                // Reads-in line 1 as total number of triangles
                tris_file >> num_tris;

                // Outputs
                cout << "Number of Trianlges = " << num_tris << endl << flush;

                // Defines array to accommodate triangles
                tris = new int*[num_tris];

                // Defines tag (if triangles have it)
                // Iterates over total number of triangles
                for(unsigned i=0;i<num_tris;i++)
                {
                        // Defines array to have 3 rows
                        tris[i] = new int[4];
                        // Inputs triangles to array (but not tag)
                        tris_file >> tris[i][0] >> tris[i][1] >> tris[i][2] >> tris[i][3];
                }
                tris_file.close();
                cout << "Done." << endl << flush;
        }

        else
        {
                cerr << "Can't Open " << file << "!" << endl << flush;
        }
}



// Function to load-in list of triangles
void loadTrisWithNtris(const char* file, int** &tris, unsigned int &num_tris, int* &n_tris, unsigned int &num_n_tris)
{
	// Puts the .tris on end of mesh base-name and loads
	string tris_end = ".tris";
	string total_tris_file = file + tris_end;
	ifstream tris_file(total_tris_file.c_str());

	// Reads-in file
	if (tris_file.is_open())
	{
		cout << "Reading " << total_tris_file << " " << endl << flush;

		// Reads-in line 1 as total number of triangles
		tris_file >> num_tris;

		// Outputs
		cout << "Number of Trianlges = " << num_tris << endl << flush;

		// Defines array to accommodate triangles
		tris = new int*[num_tris];

		// Defines tag (if triangles have it)
		int tag;

		// Iterates over total number of triangles
		for(unsigned i=0;i<num_tris;i++)
		{
			// Defines array to have 3 rows
			tris[i] = new int[3];

			// Inputs triangles to array (but not tag)
			tris_file >> tris[i][0] >> tris[i][1] >> tris[i][2] >> tag;
		}
		tris_file.close();
		cout << "Done." << endl << flush;

		//////////////////////////////////////////////////////////////
		// Puts all triangle nodes into a set (doesn't allow repeats)
		//////////////////////////////////////////////////////////////
		// Defines a set
		set<int> n_trisSet;
		set<int>::iterator it;

		// Iterates over all triangle nodes
		for(unsigned int i=0;i<num_tris;i++)
		{	
			for(int j=0;j<3;j++)
			{
				// Inserts into the set (if not already in it)
				n_trisSet.insert(tris[i][j]);
			}  
		}

		cout << "Total number of unique surface nodes = " << n_trisSet.size() << "\n";

		// Converts the set into an array for use
		int i = 0;
		int node;
		n_tris = new int [n_trisSet.size()];

		for(it=n_trisSet.begin();it!=n_trisSet.end();it++)
		{
			node = *it;
			n_tris[i] = node;	
			i++;	
		}

		// Defines number of individual surface triangles nodes for use
		num_n_tris = n_trisSet.size();

		cout << "Done. \n";
	}

	else
	{
		cerr << "Can't Open " << file << "!" << endl << flush;
	}
}

void loadFibres(const char* file, double** &fibres, unsigned int &num_elems)
{
	string elem_end = ".elem";
	string total_elem_file = file + elem_end;
	ifstream elems_file(total_elem_file.c_str());

	if (elems_file.is_open())
	{
		cout << "Reading " << total_elem_file << " " << endl << flush;
		elems_file >> num_elems;
		cout << num_elems << "\n";
		elems_file.close();
		cout << "Done." << endl << flush;
	}
	else
	{
		cerr << "Can't Open " << total_elem_file << "!" << endl << flush;
	}

	string lon_end = ".lon";
	string total_lon_file = file + lon_end;
	ifstream lon_file(total_lon_file.c_str());

//	int dummy;
//	lon_file >> dummy;

	if (lon_file.is_open())
	{
		cout << "Reading " << total_lon_file << " " << flush;

		fibres = new double*[num_elems];

		for(unsigned i=0;i<num_elems;i++)
		{
			fibres[i] = new double[3];
			lon_file >> fibres[i][0] >> fibres[i][1] >> fibres[i][2];
		}
		lon_file.close();
		cout << "Done." << endl << flush;
	}
	else
	{
		cerr << "Can't Open " << lon_file << "!" << endl << flush;
	}

}

void loadFibresNew(const char* file, double** &fibres, unsigned int &num_elems)
{
        string elem_end = ".elem";
        string total_elem_file = file + elem_end;
        ifstream elems_file(total_elem_file.c_str());

        if (elems_file.is_open())
        {
                cout << "Reading " << total_elem_file << " " << endl << flush;
                elems_file >> num_elems;
                cout << num_elems << "\n";
                elems_file.close();
                cout << "Done." << endl << flush;
        }
        else
        {
                cerr << "Can't Open " << total_elem_file << "!" << endl << flush;
        }

        string lon_end = ".lon";
        string total_lon_file = file + lon_end;
        ifstream lon_file(total_lon_file.c_str());

//      int dummy;
//      lon_file >> dummy;

        if (lon_file.is_open())
        {
                cout << "Reading " << total_lon_file << " " << flush;

                fibres = new double*[num_elems];
		int dummy ;
		lon_file >> dummy;

                for(unsigned i=0;i<num_elems;i++)
                {
                        fibres[i] = new double[3];
                        lon_file >> fibres[i][0] >> fibres[i][1] >> fibres[i][2];
                }
                lon_file.close();
                cout << "Done." << endl << flush;
        }
        else
        {
                cerr << "Can't Open " << lon_file << "!" << endl << flush;
        }

}


void loadVectors4(const char* file, double** &vectors, unsigned int num_elems)
{
        ifstream vec_file(file);

        if (vec_file.is_open())
        {
                cout << "Reading " << file << " " << endl << flush;

                vectors = new double*[num_elems];

                for(unsigned i=0;i<num_elems;i++)
                {
                        vectors[i] = new double[4];
                        vec_file >> vectors[i][0] >> vectors[i][1] >> vectors[i][2] >> vectors[i][3]; 
                }
                vec_file.close();
                cout << "Done." << endl << flush;
        }
        else
        {
                cerr << "Can't Open " << file << "!" << endl << flush;
        }

}





