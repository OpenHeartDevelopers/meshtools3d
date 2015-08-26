#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <set>
#include <list>

using namespace std;


double* calculateInwardNormal(int element, int node0, int node1, int node2, double** &coords, double** &cents)
{
	double *N;
	N = new double [3];
	double node0node1[3], node0node2[3];
	for(int j=0;j<3;j++)
	{
		node0node1[j] = coords[node0][j] - coords[node1][j];
		node0node2[j] = coords[node0][j] - coords[node2][j];
	}
	// Calculates normal by cross-product
	N[0] = (node0node1[1]*node0node2[2] - node0node1[2]*node0node2[1]);
	N[1] = (node0node1[2]*node0node2[0] - node0node1[0]*node0node2[2]);
	N[2] = (node0node1[0]*node0node2[1] - node0node1[1]*node0node2[0]);
	
	// Normalises
	double NMag = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
	for(int j=0;j<3;j++)
		N[j] = N[j]/NMag;
		
		// Calculates constant for plane equation
		double d = -(N[0]*coords[node0][0] + N[1]*coords[node0][1] + N[2]*coords[node0][2]);
		
		// Checks if the centroid of the element is a positive or negative distance from the plane
		double dist = (N[0]*cents[element][0] + N[1]*cents[element][1] + N[2]*cents[element][2] + d);
		
		// If positive, then the normal is pointing into the tissue and all is ok
		// If negative, then we need to reverse the direction of the triangle normal
		if(dist < 0)
			for(int j=0;j<3;j++)
				N[j] = -N[j];
	return(N);
				
}
	
double distanceOfPointToPlane(double* &N, int node, double** &coords, double x, double y, double z, double mu_x, double mu_y, double mu_z)
{
	double NdotQ = (N[0]*coords[node][0] + N[1]*coords[node][1] + N[2]*coords[node][2]);
	double d = -NdotQ;

	// Plug the intersection point (P') into the plane equation to get the value of a (i.e. distance of P from plane along U)
	// N.(P + Ut) + d = 0
	// => t = -(NP + d)/(N.U)
	
	double NdotP = (N[0]*x + N[1]*y + N[2]*z);
	
	double NdotU = (N[0]*mu_x + N[1]*mu_y + N[2]*mu_z);
	
	double a = -(NdotP + d)/(NdotU);
	
	return(a);
}
	
double normalDistanceOfPointToPlane(double* &N, int node, double** &coords, double x, double y, double z)
{
  double NdotQ = (N[0]*coords[node][0] + N[1]*coords[node][1] + N[2]*coords[node][2]);
  double d = -NdotQ;
  double NdotP = (N[0]*x + N[1]*y + N[2]*z);
  
  double a = -(NdotP + d);
  return(a);
}
	
	

//////////////////////////////
// Function to obtain a random number
//////////////////////////////
double GetRandomNumber()
{
  double random_number;
  int random_number_int;

  random_number_int = (rand()%100000001);
  
  random_number = (double)((random_number_int)/100000000.0);

  if(random_number < 0.00000001)
    random_number = 0.000000001;

  return (random_number);
}


//////////////////////////////
// Calculates angles of travel following scattering
//////////////////////////////
void GetScatteringAngles(double mu_x, double mu_y, double mu_z, double &mu_x_new, double &mu_y_new, double &mu_z_new, double g)
{
	double costheta, sintheta, cospsi, sinpsi, psi;
	double random_number;

	double pi = M_PI;

	// Rolls a dice for the scattering
	random_number = GetRandomNumber();

	// Finds the scattering angle theta
	costheta = (1/(2*g))*(1 + g*g - pow(((1 - g*g)/(1-g+(2*g*random_number))),2));
	sintheta = sqrt(1 - costheta*costheta);

	// Rolls a dice for the scattering
	random_number = GetRandomNumber();

	// Finds the azimuthal angle
	psi = 2*pi*random_number;

	cospsi = cos(psi);
	if(psi < pi)
		sinpsi = sqrt(1 - cospsi*cospsi);
	else
		sinpsi = -sqrt(1 - cospsi*cospsi);
	
	// Finds the new angles which we are travelling in
	if( sqrt(mu_z*mu_z) > 0.99999)
	{
		mu_x_new = sintheta*cospsi;
		mu_y_new = sintheta*sinpsi;

		if(mu_z > 0)
			mu_z_new = costheta;
		else
			mu_z_new = -costheta;      
	}
	else
	{
		mu_x_new = (sintheta/(sqrt(1 - mu_z*mu_z)))*(mu_x*mu_z*cospsi - mu_y*sinpsi) + mu_x*costheta;

		mu_y_new = (sintheta/(sqrt(1 - mu_z*mu_z)))*(mu_y*mu_z*cospsi + mu_x*sinpsi) + mu_y*costheta;

		mu_z_new = -sintheta*cospsi*sqrt((1 - mu_z*mu_z)) + mu_z*costheta;
	}
	
	double UMag = sqrt((mu_x_new*mu_x_new) + (mu_y_new*mu_y_new) + (mu_z_new*mu_z_new));

	mu_x_new = mu_x_new/UMag;
	mu_y_new = mu_y_new/UMag;
	mu_z_new = mu_z_new/UMag;

}

//////////////////////////////
// Function to map x,y,z coordinates onto voxel
//////////////////////////////
int GetCurrentPixel(double x, double y, double z, double x_res, double y_res, double z_res, int sizex, int sizey, int sizez)
{
	int p_x, p_y, p_z, p;
	
	// Translates x,y,z-coordinate into voxel by rounding-off
	p_x = floor( (x / x_res) );
	p_y = floor( (y / y_res) );
	p_z = floor( (z / z_res) );

	// Converts the x,y,z voxel positions to a single iterator
	p = p_z*sizex*sizey + p_y*sizex + p_x;
	
	return(p);
}

//////////////////////////////
// Function to Determine Reflection Coefficient R, given an Incident Angle
//////////////////////////////
double GetReflectionCoefficient(double alpha_i, double n_i, double n_t)
{
  double alpha_t, Ref;

  if(sin(alpha_i) < n_t/n_i)
    {
      alpha_t = asin(n_i*sin(alpha_i)/n_t);

      // Calculates the reflectance, R
      Ref = 0.5*(((sin(alpha_i - alpha_t))*(sin(alpha_i - alpha_t)))/((sin(alpha_i + alpha_t))*(sin(alpha_i + alpha_t))) + ((tan(alpha_i - alpha_t))*(tan(alpha_i - alpha_t)))/((tan(alpha_i + alpha_t))*(tan(alpha_i + alpha_t))));
    }
  else 
    Ref = 1.0;

  return(Ref);
}

double DeterminantCalc(double m00, double m01, double m02, double m03, double m10, double m11, double m12, double m13, double m20, double m21, double m22, double m23, double m30, double m31, double m32, double m33)
{	
	double det;
	det = m03 * m12 * m21 * m30-m02 * m13 * m21 * m30-m03 * m11 * m22 * m30+m01 * m13 * m22 * m30+ m02 * m11 * m23 * m30-m01 * m12 * m23 * m30-m03 * m12 * m20 * m31+m02 * m13 * m20 * m31 + m03 * m10 * m22 * m31-m00 * m13 * m22 * m31-m02 * m10 * m23 * m31+m00 * m12 * m23 * m31 + m03 * m11 * m20 * m32-m01 * m13 * m20 * m32-m03 * m10 * m21 * m32+m00 * m13 * m21 * m32 + m01 * m10 * m23 * m32-m00 * m11 * m23 * m32-m02 * m11 * m20 * m33+m01 * m12 * m20 * m33 + m02 * m10 * m21 * m33-m00 * m12 * m21 * m33-m01 * m10 * m22 * m33+m00 * m11 * m22 * m33;
	
	return(det);	
}

int GetCurrentElement(double x, double y, double z, double x_res, double y_res, double z_res, int sizex, int sizey, int sizez, double** &coords, int** &elems, int* &elementPixels)
{
	int p_x, p_y, p_z, p;
	int element=0, insideElement=0;

	//cout << "ELEMS = " << elems[10][0] << " " << elems[10][1] << "\n";
	//cout << "COORDS = " << coords[10][0] << " " << coords[10][1] << "\n";
	
	// Translates x,y,z-coordinate into voxel by rounding-off
	p_x = floor( (x / x_res) );
	p_y = floor( (y / y_res) );
	p_z = floor( (z / z_res) );
	//cout << "res " << x_res << " " << y_res << " " << z_res << "\n";
	//cout << "size " << sizex << " " << sizey << " " << sizez << "\n";	

	//cout << "ps " << p_x << " " << p_y << " " << p_z << "\n";

	//for(int i=0;i<1000;i++)
	//	cout << "elementPixels = " << elementPixels[i] << "\n";

	// Defines position vectors of element vertices
	double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
	int n0,n1,n2,n3;

	int d = 10;
	
	// Iterates 5 pixels either side of current pixel
	for(int k=-d;k<d;k++)
	{
		for(int j=-d;j<d;j++)
		{
			for(int i=-d;i<d;i++)
			{
				if(p_x-i >= 0 && p_y-j >= 0 && p_z-k >= 0 && p_x + i <= sizex && p_y+j <= sizey && p_z+k <= sizez)
				{
					// Converts 'pixel' coordinates into single pixel number
					p = (p_z + k)*sizex*sizey + (p_y + j)*sizex + (p_x + i);
					//cout << "p = " << p << "\n";		
					// Pixels-out which element lies in this pixel (if any)
					element = elementPixels[p];
					//cout << "element = " << element << "\n";
					
					if(element != 0)
					{
						// Picks-out node numbers of this element
						n0 = elems[element][0];
						n1 = elems[element][1];
						n2 = elems[element][2];
						n3 = elems[element][3];
					
						// Picks-out vertex coordinates of this element
						x0 = coords[n0][0];
						y0 = coords[n0][1];
						z0 = coords[n0][2];
						x1 = coords[n1][0];
						y1 = coords[n1][1];
						z1 = coords[n1][2];
						x2 = coords[n2][0];
						y2 = coords[n2][1];
						z2 = coords[n2][2];
						x3 = coords[n3][0];
						y3 = coords[n3][1];
						z3 = coords[n3][2];
						
						// Calculate 4 determinants
						double det0 = DeterminantCalc(x,y,z,1,x1,y1,z1,1,x2,y2,z2,1,x3,y3,z3,1);
						double det1 = DeterminantCalc(x0,y0,z0,1,x,y,z,1,x2,y2,z2,1,x3,y3,z3,1);
						double det2 = DeterminantCalc(x0,y0,z0,1,x1,y1,z1,1,x,y,z,1,x3,y3,z3,1);
						double det3 = DeterminantCalc(x0,y0,z0,1,x1,y1,z1,1,x2,y2,z2,1,x,y,z,1);
						
						// If all of the deterinants have the same sign, then the point lies within the tetrahedra!
						if(det0 <= 0 && det1 <= 0 && det2 <= 0 && det3 <=0)
						{
							insideElement = element;
							break;
						}
						else if(det0 >= 0 && det1 >= 0 && det2 >= 0 && det3 >=0)
						{
							insideElement = element;
							break;
						}
					}
				}
			}
		}
	}
	
	return(insideElement);
	
}	




int cornerEdgeChecker(double x, double y, double z, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
{
	int outOfCornerEdge = 0;
	int counter = 0;	

	if(x < xmin)
		counter++;
	if(x > xmax)
		counter++;
	if(y < ymin)
		counter++;
	if(y > ymax)
		counter++;
	if(z < zmin)
		counter++;
	if(z > zmax)
		counter++;
	
	if(counter > 1)
		outOfCornerEdge = 1;

	return(outOfCornerEdge);	
}

int faceChecker(double x, double y, double z, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
{
	int outOfFace = 0;

        if(x < xmin)
                outOfFace = 1;
        if(x > xmax)
                outOfFace = 2;
        if(y < ymin)
                outOfFace = 3;
        if(y > ymax)
                outOfFace = 4;
        if(z < zmin)
                outOfFace = 5;
        if(z > zmax)
                outOfFace = 6;

        return(outOfFace);

}	
