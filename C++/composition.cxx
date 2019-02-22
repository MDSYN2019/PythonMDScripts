/* --------------------------------------------------------------------------------
  ---------------------------------------------------------------------------------
  |Calculate the Clustering of the Mixed Lipid Bilayer through a DBSCAN Algorithm | 
  ---------------------------------------------------------------------------------

   How does the DBSCAN algorithm work?
   
   => Point P in a cluster is a 'core' point if there are a critical number of the same type of points within a distance E.
   => Point Q is a point in the cluster if there is a clear vector towards that from any of the elements of the P vector 
   => All points not reachable from any other points are outliers
   
   Two main algorithms are required, where there are 
   
   ----------------------------------
   | DBSCAN algorithm -  Pseudocode:|
   ----------------------------------

   DBSCAN(D, eps, MinPts) {
   C = 0
   for each point P in dataset D {
   if P is visited 
   continue next point 
   mark P as visited
   NeighbourPts = regionQuery(P, eps)
   if Sizeof(NeighbourPts) < MinPts
   mark P as NOISE
   else {
   C = next cluster
   expandCluster(P, NeighbourPts, C, eps, MinPits)
   }
   }
}

*/

/*

 -------------------------------------------------------------------------- 
 |                      Array of the Center of Masses                     | 
 --------------------------------------------------------------------------

 We are implementing a pseudo DBSCAN algorithm
    
 To measure the distance between the points, we measure the distances between 
 the centers of masses (COM). The DBSCAN algorithm originally requires to     
 calculate the `anchor' COMs, which are decided through a nearest neighbour   
 calculation. In the schematic below, 

 O  O  O ---> X  -|
  \/ \/           | -- the O's show the 'anchor', and the X's show the increasing cluster
  O  O ---- X   -|
    
 Calculating the COM of the C12E2
 
 The COM components are divided into the x, y, and z coordinates

*/

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <bits/stdc++.h> 
#include <utility>
#include <cerrno>
#include <cstdlib> 
#include <cmath>
#include <tuple>
#include <boost/progress.hpp>

//#include <Eigen>
//#include <Eigen/Dense>
//using namespace Eigen;

//double CenterOfMass(double* H7, double* H6_1, double* H6_2, double* T3_1, double* T3_2, double* T3_3, double* T4);
//double trueDist(double* COM1x, double* COM1y, double* COM1z, double* COM2x, double* COM2y, double* COM2z);

/* Function to calculate the centre of mass of each molecule */

const int numberOfPolymers = 998; // The number of polymers of each type - C12E2 or mimic 
const int numberofatoms = 71313; // Total number of beads in the simulation
const int indexCG = 7;
int numberofSS = 100; /*The number of screenshots in the dump file*/
const int boxdim = 3;

typedef struct {                                                                                                                              
  int index[7];
} C12E2_skeleton;

typedef struct {
  int a;
  int b;
  double x;
  double y;
  double z;
} inputCoord;

double trueDist(double* COM1x, double* COM1y, double* COM1z, double* COM2x, double* COM2y, double* COM2z) {
  double dist = pow((pow(COM1x-COM2x,2.0) + pow(COM1y-COM2y,2.0) + pow(COM1z-COM2z,2.0)),0.5);
  return dist;
}
 
bool compareByIndex(const inputCoord &a, const inputCoord &b) {
    return a.a < b.a;
}

/*
double CenterOfMass(C12E2_skeleton* input, int ind, std::vector<std::pair<int, double> >* vec1) {  
  // Need to update
  double COM; 
  // COM = (vec[index][input->index[0]]->second + vec[index][input->index[1]]->second + vec[index][input->index[2]]->second + vec[index][input->index[3]]->second + vec[index][input->index[4]]->second + vec[index][input->index[5]]->second + vec[index][input->index[6]]->second) / 7; 
  // With this COM definition we now know the COM in each cartesian coordinate  
  for (std::vector<std::pair<int, double> >::const_iterator it = vec1->begin() ; it != vec1->end(); it++) {
    std::cout << it->first << " " << it->second << std::endl;    
  }
  //return COM; 
}  
*/

bool sortbysec_int(const std::pair<int,int> &a, 
		   const std::pair<int,int> &b) 
{ 
    return (a.second < b.second); 
}

bool sortbysec_double(const std::pair<int,double> &a, 
		      const  std::pair<int,double> &b) 
{ 
    return (a.second < b.second); 
} 

//pointCluster* regionQuery (int , double*, double*, double*, double*, double*, double*);
//#define minClust 5 // the minimum sizes cluster we want to count as a cluster
class compute {
public:
  compute() {}; // Default constructor
  compute(std::string file) {}; // constructor reading file

  void storeFile() {    
    boost::progress_display show_progress(numberofSS);
    // open file for reading 
    ipf = fopen("dump.final", "r");  // Needs correction 
    // check for error opening file */

    if (ipf == NULL) {  
      std::cout << "Error opening file\n";
      exit(1);
    }
    // loop over the values  
    for (int SSno = 0; SSno < numberofSS; SSno++) {  

      l = 0;
      n = 0;
      index = 0;
      // Make sure to clear all vectors before going on to the next snapshot
      inputVector.clear();
      for (int k = 0; k < nlines; k++) { 
	// get a line from the file 
	// fgets() returns NULL when it reaches an error or end of file  

	fgets(line,sizeof(line),ipf);
	//while (fgets(line,sizeof(line),ipf) != NULL) {

	if (l < 5) {	
	  /* We are doing nothing */
	  //	printf("%s l= %d", line,l);
	  l++;
	}
	else if ((l > 4 && l < 8)) {
	  /*We are scanning the bit with just the box parameters*/
	  sscanf(line, "%lf %lf", &box1, &box2);
	  //	printf("%lf %lf \n",box1,box2 );
	  boxlength[l-5] = box2-box1; 
	  l++;
	}
	else if (l == 8) {
	  // printf(" **** l= %d \n",l);
	  /* We are doing nothing */
	  l++;
	}
	else {
	  //	printf(" ***** l = %d \n",l );
	  /* convert the text to numbers */
	  
	  sscanf(line,"%d %d %lf %lf %lf",&index,&atomtype,&x,&y,&z);
	  //std::cout << index << " " << atomtype << " " << x << " " << y << " " << z; 
	  //a.push_back(std::make_pair(index, index)); // Push back indices 
	  //b.push_back(std::make_pair(index, atomtype)); // Push back atomtypes
	  //xco.push_back(std::make_pair(index, x*boxlength[0])); // Push back boxlengths - x coordinates
	  //yco.push_back(std::make_pair(index, y*boxlength[1])); // Push back boxlengths - y coordinates
	  //zco.push_back(std::make_pair(index, z*boxlength[2])); // Push back boxlengths - z coordinates
	  inputline.a = index;
	  inputline.b = atomtype;
	  inputline.x = x*boxlength[0];
	  inputline.y = y*boxlength[1];
	  inputline.z = z*boxlength[2];
	  inputVector.push_back(inputline);
	  n++;
	  l++;
	  //printf("%d %d %lf %lf %lf\n",index,atomtype,x,y,z);
	}
      }
      inputTotal.push_back(inputVector);
      ++show_progress;
    }
  }

  void sortVectors () { // Sort vector values based on indices of the dump files 
    for (unsigned int i = 0; i < inputTotal.size(); ++i) {
      sort(inputTotal[i].begin(), inputTotal[i].end(), compareByIndex);
    }
  }
  
  void printVectorElements() {
     for (unsigned int i = 0; i < inputTotal.size(); ++i) {
	for (unsigned int j = 0; j < inputTotal[1].size(); ++j) {
	  //std::cout << i << " " << j << " " << " " << inputTotal[i][j].a << " " << inputTotal[i][j].b << " " <<  inputTotal[i][j].x << " " << inputTotal[i][j].y << " " << inputTotal[i][j].z << std::endl;

	}
     }
  }
 
  /*
  void CompositionProfile() {        
    //Calculating the composition profile 
    for (int i = 0; i <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1])-1; i++) { // TODO

      //double truedist(double COM1x, double COM1y, double COM1z, double COM2x, double COM2y, double COM2z) {
      // See if the distance between the nanopparticle is smaller than 15 
      //std::cout << zco[71312] << std::endl;

      if (trueDist(xco[71312], yco[71312], zco[71312], headGroupC12_E2xCOM[i], headGroupC12_E2yCOM[i], headGroupC12_E2zCOM[i]) <= 25.000) { 
	// Check if the c12e2 molecule is on the top layer or the bottom layer
	if (headGroupC12_E2zCOM[i] > zco[71312]) { // If on the top layer
	  tophead++;
	}
	else if (headGroupC12_E2zCOM[i] < zco[71312]) { // If on the bottom layer
	  downhead++;
	}
      }
      
      // This time, working with mimics 
       if(trueDist(xco[71312], yco[71312], zco[71312], headGroupMIMCxCOM[i], headGroupMIMICyCOM[i], headGroupMIMICzCOM[i] ) <= 25.000) { 
	 if (headGroupMIMICzCOM[i] > zco[71312]) {
	   mimictophead++;
	 }
	 else if (headGroupMIMICzCOM[i] < zco[71312]) {
	   mimicdownhead++;
	   }
	 }  
       }
    
    
    //std::cout << SSno << " " << tophead/(tophead + mimictophead) << " " << mimictophead/(tophead + mimictophead) << " " << downhead/(downhead + mimicdownhead) << " " << mimicdownhead/(downhead + mimicdownhead) << " " << (tophead-mimictophead)/(tophead + mimictophead) << " " << (downhead-mimicdownhead)/(downhead + mimicdownhead) << " " << ((tophead-mimictophead)/(tophead + mimictophead) + (downhead-mimicdownhead)/(downhead + mimicdownhead))/2 << " " <<  ((tophead-mimictophead)/(tophead + mimictophead) - (downhead-mimicdownhead)/(downhead + mimicdownhead))/2 <<   std::endl; 
    }
  */
  
  void check() {  
    for (unsigned int i = 0; i < inputTotal.size(); ++i) {	
      for (unsigned int j = 0; j <= inputTotal[1].size()/sizeof(inputTotal[1][1]); j++) {
	if (inputTotal[i][j].b == 7) {
	  std::cout<< inputTotal[i][j].a << " " << inputTotal[i][j].b << std::endl;
	  std::cout<< inputTotal[i][j+1].a << " " << inputTotal[i][j+1].b << std::endl;
	  std::cout<< inputTotal[i][j+2].a << " " << inputTotal[i][j+2].b << std::endl;
	  std::cout<< inputTotal[i][j+3].a << " " << inputTotal[i][j+3].b << std::endl;
	  std::cout<< inputTotal[i][j+4].a << " " << inputTotal[i][j+4].b << std::endl;
	  std::cout<< inputTotal[i][j+5].a << " " << inputTotal[i][j+5].b << std::endl;
	  std::cout<< inputTotal[i][j+6].a << " " << inputTotal[i][j+6].b << std::endl;
	} else if (inputTotal[i][j].b == 13) {
	  std::cout<< inputTotal[i][j].a << " " << inputTotal[i][j].b << std::endl;
	  std::cout<< inputTotal[i][j+1].a << " " << inputTotal[i][j+1].b << std::endl;
	  std::cout<< inputTotal[i][j+2].a << " " << inputTotal[i][j+2].b << std::endl;
	  std::cout<< inputTotal[i][j+3].a << " " << inputTotal[i][j+3].b << std::endl;
	  std::cout<< inputTotal[i][j+4].a << " " << inputTotal[i][j+4].b << std::endl;
	  std::cout<< inputTotal[i][j+5].a << " " << inputTotal[i][j+5].b << std::endl;
	  std::cout<< inputTotal[i][j+6].a << " " << inputTotal[i][j+6].b << std::endl;
	}
      }
    }
  }
  
  
  /*  
  void computeOrderphobic() {
    double dist;
    for (unsigned int i = 0; i < xcoTotal.size(); ++i) {	
      for (unsigned int j = 0; j <= sizeof(C12E2_struct)/sizeof(C12E2_struct[1]); j++) {
	CenterOfMass(&C12E2_struct[j], i, &xcoTotal[i]);	
	// std::cout << dist <<  i << " "  << j << std::endl;
      }
    }
  }
  */
  
  /*
void newProcess() {
  for (unsigned int i = 0; i < xcoTotal.size(); ++i) {
    aFormat.clear();
    bFormat.clear();
    xcoFormat.clear(); 
    ycoFormat.clear(); 
    zcoFormat.clear(); 
    for (unsigned int j = 0; j < xcoTotal[1].size(); ++j) {
      //std::cout << xcoTotal.size() << " " << xcoTotal[1].size() << "\n";
      aFormat.insert(aFormat.begin() + aTotal[i][j].first, aTotal[i][j].second);    
      bFormat.insert(bFormat.begin() + bTotal[i][j].first, bTotal[i][j].second);    
      xcoFormat.insert(xcoFormat.begin() + xcoTotal[i][j].first, xcoTotal[i][j].second);    
      ycoFormat.insert(ycoFormat.begin() + ycoTotal[i][j].first, ycoTotal[i][j].second);    
      zcoFormat.insert(zcoFormat.begin() + zcoTotal[i][j].first, zcoTotal[i][j].second);    
      std::cout << i << " " << " " << j << " " <<  aTotal[i][j].first << std::endl;
    }
    
    aTotalFormat.push_back(aFormat);
    bTotalFormat.push_back(bFormat);
    xcoTotalFormat.push_back(xcoFormat); 
    ycoTotalFormat.push_back(ycoFormat); 
    zcoTotalFormat.push_back(zcoFormat); 
  }
 
  // for (unsigned int i = 0; i < aTotalFormat.size(); ++i) {	
  // for (unsigned int j = 0; j <= aTotalFormat[1].size(); j++) {
  //    std::cout << i << " " << j << " " << aTotalFormat[i][j] << " " << bTotalFormat[i][j] << " " << xcoTotalFormat[i][j] <<  "\n";
  // }
  // }
}

  */

  /*
    std::vector<std::vector<std::vector< std::tuple<int,int, double> > > > vecOftup; // Damn ugly code!!! 
    std::vector<std::vector< std::tuple<int,int, double> > > closestDistanceVector;     
    std::tuple<int, int , double> foo; 
    std::vector< std::pair<double,int> > closestDistanceVector;
    
    for (unsigned int i = 0; i < xcoTotal.size(); ++i) {      
      for (unsigned int j = 0; j <= sizeof(C12E2M_struct)/sizeof(C12E2M_struct[1]); j++) { 	
	for (unsigned int k = 0; k <= sizeof(C12E2M_struct)/sizeof(C12E2M_struct[1]); k++) { 
	  dist = trueDist(&xcoTotal[i][C12E2_struct[j].index[3]], &ycoTotal[i][C12E2_struct[j].index[3]], &zcoTotal[i][C12E2_struct[j].index[3]],
			  &xcoTotal[i][C12E2_struct[k].index[3]], &ycoTotal[i][C12E2_struct[k].index[3]], &zcoTotal[i][C12E2_struct[k].index[3]]);
	  foo = std::make_tuple(j, k, dist);
	  closestDistanceVector.append(foo);	  
	}
      }
      vecOftup.append(closestDistanceVector);
      closestDistanceVector.clear(); // Clear after accumulating
    }
    */
  /*
  double ChainOrderAnalysis(double* H7_x, double* H7_y, double* H7_z, double* H6_1_x, double* H6_1_y, double* H6_1_z,  double* H6_2_x, double* H6_2_y, double* H6_2_z, double* T3_1_x, double* T3_1_y, double* T3_1_z,  double* T3_2_x, double* T3_2_y, double* T3_2_z, double* T3_3_x, double* T3_3_y, double* T3_3_z, double* T4_x, double* T4_y, double* T4_z) {

    // This bit is obsolete - needs to be modified

    double H7coord = H7;
    double H6coord_1 = H6_1;
    double H6coord_2 = H6_2;
    double T3coord_1 = T3_1; 
    double T3coord_2 = T3_2; 
    double T3coord_3 = T3_3; 
    double T4coord = T4; 
    double COM;
  }
  */
  
  /*
  void clusterAnalysis(std::vector<double>* C12E2xCOM, vector<double>* C12E2yCOM, vector<double>* C12E2zCOM, std::vector<C12E2_skeleton>* values, double* distance, int* PolymerCounter) { 
      for (unsigned i : indices(C12E2xCOM)) { // 
	for (unsigned j : indices(C12E2xCOM)) { // 
	  if (C12E2xCOM[i] != C12E2xCOM[j] && C12E2yCOM[i] != C12E2yCOM[j] && C12E2zCOM[i] != C12E2zCOM[j]) {	    
	    if (C12E2xCOM[i] - C12E2xCOM[j] !=0 && sqrt(pow(C12E2xCOM[i] - C12E2xCOM[j],2)) <= distance) {
	      if (C12E2yCOM[i] - C12E2yCOM[j] !=0 && sqrt(pow(C12E2yCOM[i] - C12E2xCOMy[j],2)) <= distance) {
		if (C12E2zCOM[i] - C12E2zCOM[j] !=0 && sqrt(pow(C12E2zCOM[i] - C12E2zCOM[j],2)) <= distance) {
		  PolymerCounter++; // Pointer as there will be two counters 
		}
	      }
	    }
	  }
	}
      }
    }
  */
private:    
  // Vectors to store trajectory values 
  std::vector<inputCoord> inputVector;
  std::vector<std::vector<inputCoord> > inputTotal;
  
  FILE *ipf; /* input file */  
  // COM vectors - C12E2
  std::vector<double> headGroupC12_E2xCOM; // COM C12E2 X
  std::vector<double> headGroupC12_E2yCOM; // COM C12E2 Y
  std::vector<double> headGroupC12_E2zCOM; // COM C12E2 Z
  // COM vectors - C12E2-M
  std::vector<double> headGroupMIMICxCOM; // COM MIMIC X
  std::vector<double> headGroupMIMICyCOM; // COM MIMIC Y
  std::vector<double> headGroupMIMICzCOM; // COM MIMIC Z
  
  // -----------------------------------------------------------------------//
  // Defining normals (A7 -- A6_1 -- A6_2 -- A3_1 -- A3_2 -- A3_3 -- A4     //
  // -----------------------------------------------------------------------//
  int A7, A6_1, A6_2, A3_1, A3_2, A3_3, A4; // First Batch
  int A7_2, A6_1_2, A6_2_2, A3_1_2, A3_2_2, A3_3_2, A4_2; // Second Batch
  int A7_3, A6_1_3, A6_2_3, A3_1_3, A3_2_3, A3_3_3, A4_3; // Third Batch 
  int A7_4, A6_1_4, A6_2_4, A3_1_4, A3_2_4, A3_3_4, A4_4; // Fourth Batch
  
  // ----------------------------------------------------------------------//
  // Defining mimics (A13 -- A12_1 -- A12_2 -- A9_1 -- A9_2 -- A9_3 -- A10 //
  // ----------------------------------------------------------------------//  
  int A13, A12_1, A12_2, A9_1, A9_2, A9_3, A10; // First Batch
  int A13_2, A12_1_2, A12_2_2, A9_1_2, A9_2_2, A9_3_2, A10_2; // Second Batch
  int A13_3, A12_1_3, A12_2_3, A9_1_3, A9_2_3, A9_3_3, A10_3; // Third Batch
  int A13_4, A12_1_4, A12_2_4, A9_1_4, A9_2_4, A9_3_4, A10_4; // Fourth Batch

  C12E2_skeleton C12E2_struct[numberOfPolymers];                                                                                              
  C12E2_skeleton C12E2M_struct[numberOfPolymers];    
 
  double tophead = 0;
  double downhead = 0;
  double mimictophead = 0;
  double mimicdownhead = 0;
  
  int MimicCounter = 0;    
  int PolymerCounter = 0;    
  int atomtype; 
  int index, l, n;
  int nlines = numberofatoms + 9;      
  double x,y,z; /*coordinates for the atoms in the box*/
  double box1;
  double box2; 
  double boxlength[boxdim];
  char line[100];  
  inputCoord inputline;

};
 
compute A;
  
int main (int argc, char *argv[])  {
  A.storeFile();
  A.sortVectors();
  //  A.printVectorElements();
  A.check();
  //  A.sortVectors();
  //A.AllocationVec();
  //  A.computeOrderphobic();
  // A.newProcess();
  return 0;    
}



