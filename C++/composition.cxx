// ------------------------------------------------------------------------- //
//                       Array of the center of masses                       //
// ------------------------------------------------------------------------- //

/*
 We are implementing a pseudo DBSCAN algorithm
    
 To measure the distance between the points, we measure the distances between //
 the centers of masses (COM). The DBSCAN algorithm originally requires to     //
 calculate the `anchor' COMs, which are decided through a nearest neighbour   //
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
#include <utility>
#include <cerrno>
#include <cstdlib> 
#include <cmath>
#include <string>
#include <boost/progress.hpp>

#include "compositon.h"


//#include <Eigen>
//#include <Eigen/Dense>
//using namespace Eigen;

//double CenterOfMass(double* H7, double* H6_1, double* H6_2, double* T3_1, double* T3_2, double* T3_3, double* T4);
//double trueDist(double* COM1x, double* COM1y, double* COM1z, double* COM2x, double* COM2y, double* COM2z);

/* Function to calculate the centre of mass of each molecule */

int numberofSS = 100; /*The number of screenshots in the dump file*/
const int numberOfPolymers = 1000; // The number of polymers of each type - C12E2 or mimic 
const int numberofatoms = 71313; // Total number of beads in the simulation
const int indexCG = 7;
const int boxdim = 3;

struct C12E2_skeleton {                                                                                                                              
  int index[7];
};


double CenterOfMass(double* H7, double* H6_1, double* H6_2, double* T3_1, double* T3_2, double* T3_3, double* T4) {  
  // Need to update
  double COM; 
  COM = (*H7 + *H6_1 + *H6_2 + *T3_1 + *T3_2 + *T3_3 + *T4)/7; 
  /* With this COM definition we now know the COM in each cartesian coordinate */ 
  return COM; 
}  

double trueDist(double* COM1x, double* COM1y, double* COM1z, double* COM2x, double* COM2y, double* COM2z) {
  double dist = pow((pow(COM1x-COM2x,2.0) + pow(COM1y-COM2y,2.0) + pow(COM1z-COM2z,2.0)),0.5);
  return dist;
}

void die (const char *message) {
  if (errno) {
    perror(message);
  } else {
    std::cout << message << std::endl;
  }
  exit(1);
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
      //  printf("This is the data for trajectory no %d \n", SSno); 
      l = 0;
      n = 0;
      index = 0;
      a.clear();
      b.clear();
      xco.clear();
      yco.clear();
      zco.clear();
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
	  a.push_back(index); // Push back indices 
	  b.push_back(atomtype); // Push back atomtypes
	  xco.push_back(x*boxlength[0]); // Push back boxlengths - x coordinates
	  yco.push_back(y*boxlength[1]); // Push back boxlengths - y coordinates
	  zco.push_back(z*boxlength[2]); // Push back boxlengths - z coordinates
	  n++;
	  l++;
	  //printf("%d %d %lf %lf %lf\n",index,atomtype,x,y,z);
	}
      }
      aTotal.push_back(a);
      bTotal.push_back(b);
      xcoTotal.push_back(xco);
      ycoTotal.push_back(yco);
      zcoTotal.push_back(zco);
      ++show_progress;
    }
  }

  void printVectorElements() {
    for (unsigned int i = 0; i < xcoTotal.size(); ++i)
      {
	for (unsigned int j = 0; j < xcoTotal[i].size(); ++j)
	  {
	    std::cout << i << " " << j << " " << " " << xcoTotal[i][j] << std::endl;
	  }
	
	//std::cout << std::endl;
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
 
  void AllocationVec() {
    /*
      The way that the C12E2 and C12E2-M batches are allocated in the simulation follows
      the fact that this was directly replicated in LAMMPS; hence, there isn't a single 
      block of list, but rather 4 batches of C12E2 and C12E2-M portions that we need to 
      allocate 
      
      Each batch has 250 molecules, which is why we loop from 0 to 249 
    */
    
    for (int i = 0; i <= 249; i++) {

      // C12E2 batches

      int batch1 = 0;
      int batch2 = 3500;
      int batch3 = 7000;
      int batch4 = 10500;

      // C12E2-M batches

      int batch1M = 1750;
      int batch2M = 5250;	
      int batch3M = 8750;
      int batch4M = 12250;
      
      // Allocate index

      A7 = 7*(i);
      A6_1= 7*(i) + 1;
      A6_2= 7*(i) + 2;
      A3_1 = 7*(i) + 3;
      A3_2 = 7*(i) + 4;
      A3_3 = 7*(i) + 5;
      A4 = 7*(i) + 6;
      
      // Second batch
      
      A7_2 = 7*(i) + batch2 + 1;
      A6_1_2 = 7*(i) + batch2 + 2;
      A6_2_2 = 7*(i) + batch2 + 3;
      A3_1_2 = 7*(i) + batch2 + 4;
      A3_2_2 = 7*(i) + batch2 + 5;
      A3_3_2 = 7*(i) + batch2 + 6;
      A4_2 = 7*(i) + batch2 + 7;
      
      // Third batch
      
      A7_3 = 7*(i) + batch3 + 1;
      A6_1_3 = 7*(i) + batch3 + 2;
      A6_2_3 = 7*(i) + batch3 + 3;
      A3_1_3 = 7*(i) + batch3 + 4;
      A3_2_3 = 7*(i) + batch3 + 5;
      A3_3_3 = 7*(i) + batch3 + 6;
      A4_3 = 7*(i) + batch3 + 7;
      
      // Fourth batch
      
      A7_4 = 7*(i)  + batch4 + 1;
      A6_1_4 = 7*(i) + batch4 + 2;
      A6_2_4 = 7*(i) + batch4 + 3;
      A3_1_4 = 7*(i) + batch4 + 4;
      A3_2_4 = 7*(i) + batch4 + 5;
      A3_3_4 = 7*(i) + batch4 + 6;
      A4_4 = 7*(i) + batch4 + 7;
      
      //  --- C12E2-M --- //
    
      A13 = 7*(i) + batch1M + 1;
      A12_1= 7*(i) + batch1M + 2;
      A12_2= 7*(i) + batch1M + 3;
      A9_1 = 7*(i) + batch1M + 4;
      A9_2 = 7*(i) + batch1M + 5;
      A9_3 = 7*(i) + batch1M + 6;
      A10 = 7*(i) + batch1M + 7;
    
      /* Mimic atoms - Second batch */
      
      A13_2 = 7*(i) + batch2M + 1;
      A12_1_2 = 7*(i) + batch2M + 2;
      A12_2_2 = 7*(i) + batch2M + 3;
      A9_1_2 = 7*(i) + batch2M + 4;
      A9_2_2 = 7*(i) + batch2M + 5;
      A9_3_2 = 7*(i) + batch2M + 6;
      A10_2 = 7*(i) + batch2M + 7;
      
      /* Mimic atoms - Third batch */
    
      A13_3 = 7*(i) + batch3M + 1;
      A12_1_3 = 7*(i) + batch3M + 2;
      A12_2_3 = 7*(i) + batch3M + 3;
      A9_1_3 = 7*(i) + batch3M + 4;
      A9_2_3 = 7*(i) + batch3M + 5;
      A9_3_3 = 7*(i) + batch3M + 6;
      A10_3 = 7*(i) + batch3M + 7;
      
      /* Mimic atoms - Fourth batch */
      
      A13_4 = 7*(i) + batch4M + 1;
      A12_1_4 = 7*(i) + batch4M + 2;
      A12_2_4 = 7*(i) + batch4M + 3;
      A9_1_4 = 7*(i) + batch4M + 4;
      A9_2_4 = 7*(i) + batch4M + 5;
      A9_3_4 = 7*(i) + batch4M + 6;
      A10_4 = 7*(i) + batch4M + 7;   
      
      
      // ----------- C12E2 indicies ----------- //
      
      /*
	Now that we know the indiices, allocate each index of the chain into the
	struct, so that we can pick out what is needed. 
      */
	
      // First batch
      C12E2_struct[i].index[0] = A7;
      C12E2_struct[i].index[1] = A6_1;
      C12E2_struct[i].index[2] = A6_2;
      C12E2_struct[i].index[3] = A3_1;
      C12E2_struct[i].index[4] = A3_2;
      C12E2_struct[i].index[5] = A3_3;
      C12E2_struct[i].index[6] = A4;

      // Second batch
      C12E2_struct[i+250].index[0] = A7_2;
      C12E2_struct[i+250].index[1] = A6_1_2;
      C12E2_struct[i+250].index[2] = A6_2_2;
      C12E2_struct[i+250].index[3] = A3_1_2;
      C12E2_struct[i+250].index[4] = A3_2_2;
      C12E2_struct[i+250].index[5] = A3_3_2;
      C12E2_struct[i+250].index[6] = A4_2;
      
      // Third batch 
      C12E2_struct[i+500].index[0] = A7_3;
      C12E2_struct[i+500].index[1] = A6_1_3;
      C12E2_struct[i+500].index[2] = A6_2_3;
      C12E2_struct[i+500].index[3] = A3_1_3;
      C12E2_struct[i+500].index[4] = A3_2_3;
      C12E2_struct[i+500].index[5] = A3_3_3;
      C12E2_struct[i+500].index[6] = A4_3;
      
      // Fourth batch
      C12E2_struct[i+750].index[0] = A7_4;
      C12E2_struct[i+750].index[1] = A6_1_4;
      C12E2_struct[i+750].index[2] = A6_2_4;
      C12E2_struct[i+750].index[3] = A3_1_4;
      C12E2_struct[i+750].index[4] = A3_2_4;
      C12E2_struct[i+750].index[5] = A3_3_4;
      C12E2_struct[i+750].index[6] = A4_4;
      
      // Mimic assignment
      
      
      // First batch
      C12E2M_struct[i].index[0] = A13;
      C12E2M_struct[i].index[1] = A12_1;
      C12E2M_struct[i].index[2] = A12_2;
      C12E2M_struct[i].index[3] = A9_1;
      C12E2M_struct[i].index[4] = A9_2;
      C12E2M_struct[i].index[5] = A9_3;
      C12E2M_struct[i].index[6] = A10;
      
      // Second batch
      C12E2M_struct[i+250].index[0] = A13_2;
      C12E2M_struct[i+250].index[1] = A12_1_2;
      C12E2M_struct[i+250].index[2] = A12_2_2;
      C12E2M_struct[i+250].index[3] = A9_1_2;
      C12E2M_struct[i+250].index[4] = A9_2_2;
      C12E2M_struct[i+250].index[5] = A9_3_2;
      C12E2M_struct[i+250].index[6] = A10_2;
      
      // Third batch 
      C12E2M_struct[i+500].index[0] = A13_3;
      C12E2M_struct[i+500].index[1] = A12_1_3;
      C12E2M_struct[i+500].index[2] = A12_2_3;
      C12E2M_struct[i+500].index[3] = A9_1_3;
      C12E2M_struct[i+500].index[4] = A9_2_3;
      C12E2M_struct[i+500].index[5] = A9_3_3;
      C12E2M_struct[i+500].index[6] = A10_3;
      
      // Fourth batch
      C12E2M_struct[i+750].index[0] = A13_4;
      C12E2M_struct[i+750].index[1] = A12_1_4;
      C12E2M_struct[i+750].index[2] = A12_2_4;
      C12E2M_struct[i+750].index[3] = A9_1_4;
      C12E2M_struct[i+750].index[4] = A9_2_4;
      C12E2M_struct[i+750].index[5] = A9_3_4;
      C12E2M_struct[i+750].index[6] = A10_4; 	
    }
  }

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

  std::vector<double> xco; 
  std::vector<double> yco; 
  std::vector<double> zco; 
  std::vector<int> a;
  std::vector<int> b;

  std::vector< std::vector<double>> xcoTotal; 
  std::vector< std::vector<double>> ycoTotal; 
  std::vector< std::vector<double>> zcoTotal; 
  std::vector< std::vector<int>> aTotal;
  std::vector< std::vector<int>> bTotal;
  //double xco[numberofatoms],yco[numberofatoms],zco[numberofatoms]; 
  /* ditto for a and b, which represent index and atomtype respectively */
  //int a[numberofatoms],b[numberofatoms];
  // double xco[numberofatoms], yco[numberofatoms], zco[numberofatoms]; 
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
  struct C12E2_skeleton C12E2_struct[numberOfPolymers];                                                                                              
  struct C12E2_skeleton C12E2M_struct[numberOfPolymers];    
 

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
};

compute A;
  
int main (int argc, char *argv[])  {
  
  A.storeFile();
  A.printVectorElements();
  
  return 0;    
}


