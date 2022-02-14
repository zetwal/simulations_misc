#include <cstdlib>
#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdexcept>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <fstream>
#include <algorithm>    // std::sort
#include <sstream>
#include <omp.h>
#include <sys/stat.h>
#include "json.hpp"

// Generic IO
#include "GenericIO.h"
#include "Partition.h"

#include "Halos_test.h"

//#include "MurmurHashNeutral2.cpp" 
#include "routines.h"

/// Assumes ids and  are consistent but ordering may have changed
/// Redistributes aongst ranks and sorts the values for a one-to-one check of the changes in 
/// several halo outputs


// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

int main( int argc, char** argv ) {

  MPI_Init( &argc, &argv );
  Partition::initialize();
  GenericIO::setNaturalDefaultPartition();


  //
	// Input
	nlohmann::json jsonInput;
	std::ifstream jsonFile(argv[1]);
	jsonFile >> jsonInput;



  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  string path1, path2;
  string haloProperties1, haloProperties2;
  string sodFile1, sodFile2;
  string partFile1, partFile2;
  string partFile1a, partFile2a;
  string partFile1b, partFile2b;

  int halo_opt, sod_opt, part_opt;
  float lim, mass_min, mass_max, box_size;


  halo_opt = jsonInput["halo-opt"];
  sod_opt = jsonInput["sod-opt"];
  part_opt = jsonInput["part-opt"];
  
  box_size = jsonInput["box-size"];
  lim = jsonInput["threshold-fraction"];

  mass_min = jsonInput["minimum-mass-for-position-matching"];
  mass_max = jsonInput["maximum-mass-for-position-matching"];



  path1 = jsonInput["file-paths"]["orig"];
  path2 = jsonInput["file-paths"]["orig"];

  haloProperties1 = jsonInput["haloproperties"]["orig"];
  haloProperties2 = jsonInput["haloproperties"]["comp"];

  sodFile1 = jsonInput["sodfile"]["orig"];
  sodFile2 = jsonInput["sodfile"]["comp"];

  partFile1 = jsonInput["part-file"]["orig"];
  partFile2 = jsonInput["part-file"]["comp"];

  partFile1a = jsonInput["part-file-a"]["orig"];
  partFile2a = jsonInput["part-file-a"]["comp"];

  partFile1b = jsonInput["part-file-b"]["orig"];
  partFile2b = jsonInput["part-file-b"]["comp"];


  if (rank==0){
    cout << "Parameter list"<< endl;
    cout << "path1: " << path1 << endl;
    cout << "path2: " << path2 << endl;
    cout << "halo_opt: " << halo_opt << endl;
    cout << "sod_opt: " << sod_opt << endl;
    cout << "part_opt: " << part_opt << endl;
    cout << "box_size: " << box_size << endl;
    cout << "lim: " << lim << endl;
    cout << "mass_min: " << mass_min <<endl;
    cout << "mass_max: " << mass_max << endl;
    cout << endl;
  }


  string fof_file, fof_file2;
  fof_file = path1 + haloProperties1;   // "/orig-499.haloproperties";
  fof_file2 = path2 + haloProperties2; // "/sz_001_velNoComp-499.haloproperties";

  string sod_file, sod_file2;
  sod_file = path1 + sodFile1;   // "/orig-delta200.0-499.sodpropertybins";
  sod_file2 = path2 + sodFile2; // "/sz_001_velNoComp-delta200.0-499.sodpropertybins";

  map<int64_t,int> tag_map_main;
  map<int64_t,int>* tag_map = &tag_map_main; // halos passing the mass thresholds

  string part_file, part_file2;
  part_file = path1 + partFile1;   // "/orig-delta200.0-499.bighaloparticles#104";
  part_file2 = path2 + partFile2; // "/sz_001_velNoComp-499.bighaloparticles#104";

  int err=0;


  if (halo_opt==0){
    err = perform_halo_check (fof_file, fof_file2, lim, mass_min, mass_max,tag_map);
  }
  else if (halo_opt==1)
  {
    err = match_pos(fof_file, fof_file2, lim, box_size, mass_min, mass_max);
  }
  else if (halo_opt==2)
    err = compare_dist(fof_file, fof_file2, lim);
  if ((sod_opt==0)&&(halo_opt==0))
    err = sodbin_check(sod_file,sod_file2,lim,tag_map);
  if (part_opt==0)
    err = part_check(part_file, part_file2);
  if (part_opt==1){
    err = part_check(part_file, part_file2);
    part_file = path1 + partFile1a;   
    part_file2 = path2 + partFile2a;  
    err = part_check(part_file, part_file2);
    part_file = path1 + partFile1b;   
    part_file2 = path2 + partFile2b;  
    err = part_check(part_file, part_file2);
  }


  Partition::finalize();
  MPI_Finalize();
  return 0;
}


// source ../env.sh
// export OMP_NUM_THREADS=1
// make clean
// make run_test
// mpirun -n 8 ./run_test inputs/json_file.json > log_name.log
