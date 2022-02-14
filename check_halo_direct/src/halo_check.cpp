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
#include <algorithm> // std::sort
#include <sstream>
#include <omp.h>
#include <sys/stat.h>

// Generic IO
#include "GenericIO.h"
#include "Partition.h"
#include "Halos_test.h"
#include "MurmurHashNeutral2.h"
#include "routines.h"

/// Assumes ids and  are consistent but ordering may have changed
/// Redistributes aongst ranks and sorts the values for a one-to-one check of the changes in
/// several halo outputs

// Cosmotools
using namespace std;
using namespace gio;
using namespace cosmotk;

bool comp_by_halo_dest(const halo_properties_test &a, const halo_properties_test &b)
{
  return a.rank < b.rank;
}

bool comp_by_halo_tag(const halo_properties_test &a, const halo_properties_test &b)
{
  return a.fof_halo_tag < b.fof_halo_tag;
}

inline void write_to_file(std::string filename, std::string text)
{
  std::fstream file;

  //std::cout << "@@ " << text << std::endl;

  file.open(filename, std::ios_base::app | std::ios_base::in);
  if (file.is_open())
    file << text;
}

int compute_mean_float(vector<float> *val1, vector<float> *val2, int num_halos, string var_name, string var_name2, float lim)
{
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  double diff = 0;
  double diff_frac = 0;
  double frac_max = 0;


  double mean;
  double stddev = 0;
  double stddevq = 0;
  double meanq = 0;
  double meanq_tot;
  int n_tot;

  for (int i=0; i<num_halos; i++){
      diff += (double)(val1->at(i)-val2->at(i));
      meanq += (double)(val1->at(i));
      if (fabs(val1->at(i))>1.e-6){
	      double frac = 0;
	      frac = (double)(fabs(val1->at(i)-val2->at(i))/fabs(val1->at(i)));
        diff_frac = (diff_frac<frac)?frac:diff_frac;
      }
   }

      MPI_Barrier(Partition::getComm());
      MPI_Allreduce(&diff, &mean, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
      MPI_Allreduce(&diff_frac, &frac_max, 1, MPI_DOUBLE, MPI_MAX,  Partition::getComm());
      MPI_Allreduce(&meanq, &meanq_tot, 1, MPI_DOUBLE, MPI_SUM, Partition::getComm());
      MPI_Allreduce(&num_halos, &n_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());

   mean = mean/n_tot;
   meanq_tot = meanq_tot/n_tot;

   double diffq;
   double diff_tmp;
   for (int i=0; i< num_halos; i++){
      diff_tmp = (double)(val1->at(i)-val2->at(i))-mean;
      diffq = (double)(val1->at(i))- meanq_tot;
      stddev += diff_tmp*diff_tmp/(n_tot-1);
      stddevq += diffq*diffq/(n_tot-1);
   }
   double stddev_tot;
   double stddevq_tot;
   MPI_Allreduce(&stddev, &stddev_tot, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   MPI_Allreduce(&stddevq, &stddevq_tot, 1, MPI_DOUBLE, MPI_SUM,  Partition::getComm());
   stddev_tot = sqrt(stddev_tot);
   stddevq_tot = sqrt(stddevq_tot);

   bool print_out = true;
   if (rank==0){

      if ((frac_max<lim)||((fabs(stddev_tot/stddevq_tot)<lim)&&(fabs(mean/meanq_tot)<lim))) // no values change by more than lim
        print_out=false;
      if (isnan(mean)) // if nan is present, then this is higher than threshold
        print_out=true;
        
      if (print_out){
        cout << " " << var_name << endl;
        cout << " " << var_name2 << endl;
        cout << " ______________________________________" <<endl;

        cout << " mean difference = " << mean << endl;
        write_to_file("outputs/test_mean.yaml", " " + var_name + ": " + std::to_string(mean) + "\n");

        cout << " maximum fractional difference = " << frac_max << endl;
        write_to_file("outputs/test_fracDiff.yaml", " " + var_name + ": " + std::to_string(frac_max) + "\n");

        cout << " standard deviation of difference = " << stddev_tot << endl;
        write_to_file("outputs/test_stdDevDiff.yaml", " " + var_name + ": " + std::to_string(stddev_tot) + "\n");

        cout << " mean of quantity = " << meanq_tot << endl;
        write_to_file("outputs/test_meanQty.yaml", " " + var_name + ": " + std::to_string(meanq_tot) + "\n");

        cout << " standard deviation of quantity = " << stddevq_tot << endl;
        write_to_file("outputs/test_stdDevQty.yaml", " " + var_name + ": " + std::to_string(stddevq_tot) + "\n");
        
        cout << endl;
        return 1;
      }
   }
   return 0;
}


int compute_mean_std_dist_ellipticity(Halos_test H_1 , Halos_test H_2, float lim ){
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = H_1.num_halos;


  int err=0;
  string e1_name = "axis_ratio_13";
  string e2_name = "axis_ratio_23";

  err += compute_mean_float(H_1.ellipticity_data[0],H_2.ellipticity_data[0],count,e1_name,e1_name,lim);
  err += compute_mean_float(H_1.ellipticity_data[1],H_2.ellipticity_data[1],count,e2_name,e2_name,lim);
  return err;
}



int  compute_mean_std_dist_halo(Halos_test H_1 , Halos_test H_2, float lim ){
  // compute the mean and std of the differences and make a histogram to look more closely
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  int count = H_1.num_halos;

  if (rank == 0)
  {
  write_to_file("outputs/test_mean.yaml",      "-\n");
  write_to_file("outputs/test_fracDiff.yaml",  "-\n");
  write_to_file("outputs/test_stdDevDiff.yaml","-\n");
  write_to_file("outputs/test_meanQty.yaml",   "-\n");
  write_to_file("outputs/test_stdDevQty.yaml", "-\n");
  }

  int err = 0;
  for (int i = 0; i < N_HALO_FLOATS; i++)
  {
    string var_name = float_var_names_test[i];
    string var_name2 = float_var_names_test2[i];
    err += compute_mean_float(H_1.float_data[i], H_2.float_data[i], count, var_name, var_name2, lim);
  }
  return err;
}

inline unsigned int tag_to_rank(int64_t fof_tag, int n_ranks)
{
  return MurmurHashNeutral2((void *)(&fof_tag), sizeof(int64_t), 0) % n_ranks;
}

halo_properties_test alter_ellipticities(halo_properties_test tmp){
    

    float s1x = tmp.ellipticity_data[0];
    float s1y = tmp.ellipticity_data[1];
    float s1z = tmp.ellipticity_data[2];
    float s2x = tmp.ellipticity_data[3];
    float s2y = tmp.ellipticity_data[4];
    float s2z = tmp.ellipticity_data[5];
    float s3x = tmp.ellipticity_data[6];
    float s3y = tmp.ellipticity_data[7];
    float s3z = tmp.ellipticity_data[8];

    float e1 = sqrt(s1x*s1x+s1y*s1y+s1z*s1z);
    float e2 = sqrt(s2x*s2x+s2y*s2y+s2z*s2z);
    float e3 = sqrt(s3x*s3x+s3y*s3y+s3z*s3z);

    for (int i = 0; i<N_HALO_FLOATS_E; i++){
        tmp.ellipticity_data[i]  = 0.0;
    }
    tmp.ellipticity_data[0] = e1/e3;
    tmp.ellipticity_data[1] = e2/e3;
	    

    if (N_HALO_FLOATS_E==18){
    float r1x = tmp.ellipticity_data[9];
    float r1y = tmp.ellipticity_data[10];
    float r1z = tmp.ellipticity_data[11];
    float r2x = tmp.ellipticity_data[12];
    float r2y = tmp.ellipticity_data[13];
    float r2z = tmp.ellipticity_data[14];
    float r3x = tmp.ellipticity_data[15];
    float r3y = tmp.ellipticity_data[16];
    float r3z = tmp.ellipticity_data[17];
    }
    return tmp;
}

void read_halos(Halos_test &H0, string file_name, int file_opt) {
 // Read halo files into a buffer
  GenericIO GIO(Partition::getComm(),file_name,GenericIO::FileIOMPI);
  GIO.openAndReadHeader(GenericIO::MismatchRedistribute);
  size_t num_elems = GIO.readNumElems();

  H0.Resize(num_elems + GIO.requestedExtraSpace()); 

  GIO.addVariable("fof_halo_tag",   *(H0.fof_halo_tag),true);
  GIO.addVariable("fof_halo_count", *(H0.fof_halo_count), true);
  if (H0.has_sod)
    GIO.addVariable("sod_halo_count", *(H0.sod_halo_count), true);

 if (file_opt==1){
  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names_test[i], *(H0.float_data[i]), true);
   }
 else{
  for (int i=0; i<N_HALO_FLOATS; ++i)
    GIO.addVariable((const string)float_var_names_test2[i], *(H0.float_data[i]), true);
 }
  for (int i=0; i<N_HALO_FLOATS_E; ++i)
    GIO.addVariable((const string)float_var_names_ellipticity[i], *(H0.ellipticity_data[i]), true);
   

  GIO.readData();
  H0.Resize(num_elems);
}  



int perform_halo_check(string fof_file, string fof_file2, float lim, float min_mass, float max_mass, map<int64_t, int> *tag_map){
  int rank, n_ranks;
  rank = Partition::getMyProc();
  n_ranks = Partition::getNumProc();

  if (rank==0){
	  cout << " Performing halo check with ID matching "<< endl;
          cout << " ______________________________________" <<endl;
          cout << endl;
  }
  // create halo buffers
  Halos_test H_1;
  Halos_test H_2;
  Halos_test H_3;

  H_1.Allocate();
  H_1.has_sod = true;
  H_2.Allocate();
  H_2.has_sod = true;
  H_1.Set_MPIType();
  H_2.Set_MPIType();
  H_3.Allocate();
  H_3.has_sod = true;
  H_3.Set_MPIType();



  read_halos(H_1, fof_file, 2);
  read_halos(H_2, fof_file2, 2);
  if (rank==0)
      cout << endl;
  // determine destination ranks
  vector<halo_properties_test> fof_halo_send;
  vector<int> fof_halo_send_cnt(n_ranks,0);
  vector<halo_properties_test> fof_halo_send2;
  vector<int> fof_halo_send_cnt2(n_ranks,0);
  vector<int64_t> fof_halo_tags_masked;
  vector<int64_t> fof_halo_tags_masked_full;

  //fof_halo_tags_masked.resize(0);

  int tag_send_count=0;
  for (int i=0; i<H_1.num_halos; ++i) {
    halo_properties_test tmp = H_1.GetProperties(i);
    if ((tmp.float_data[0]>min_mass)&&(tmp.float_data[0]<max_mass)){
    tmp.rank = tag_to_rank(tmp.fof_halo_tag, n_ranks);
    tmp = alter_ellipticities(tmp);
    fof_halo_send.push_back(tmp);
    fof_halo_tags_masked.push_back(tmp.fof_halo_tag);
    ++fof_halo_send_cnt[tmp.rank];
    ++tag_send_count;
    }
  }
  for (int i=0; i<H_2.num_halos; ++i) {
    halo_properties_test tmp = H_2.GetProperties(i);
    if ((tmp.float_data[0]>min_mass)&&(tmp.float_data[0]<max_mass)){
    tmp.rank = tag_to_rank(tmp.fof_halo_tag, n_ranks);
    tmp = alter_ellipticities(tmp);
    fof_halo_send2.push_back(tmp);
    ++fof_halo_send_cnt2[tmp.rank];
    }
  }

  
  // sort by destination rank
  sort(fof_halo_send.begin(),fof_halo_send.end(), comp_by_halo_dest);
  sort(fof_halo_send2.begin(),fof_halo_send2.end(), comp_by_halo_dest);
  MPI_Barrier(Partition::getComm());



  H_1.Resize(0);
  H_2.Resize(0);
  H_3.Resize(0);
  // create send and receive buffers and offsets
  vector<int> fof_halo_recv_cnt;
  vector<int> fof_halo_recv_cnt2;
  fof_halo_recv_cnt.resize(n_ranks,0);
  fof_halo_recv_cnt2.resize(n_ranks,0);
  
  int sizes[n_ranks];
  vector<int> tag_disp(n_ranks,0);
  MPI_Allgather(&tag_send_count,1,MPI_INT,sizes,1,MPI_INT,Partition::getComm());
  int sum=0;
  for(int i=0; i < n_ranks; i++){
      sum+=sizes[i];
      if(i==0){
         tag_disp[i]=0;
      }else{
         tag_disp[i]=tag_disp[i-1]+sizes[i-1];
      }
  }

  fof_halo_tags_masked_full.resize(sum);

  int size = (int)fof_halo_tags_masked.size();
  MPI_Allgatherv(&fof_halo_tags_masked[0], size, MPI_INT64_T, &fof_halo_tags_masked_full[0],sizes, &tag_disp[0], MPI_INT64_T,Partition::getComm());



  for (int i=0;i<sum;i++){
  (*tag_map).insert(pair<int64_t,int>(fof_halo_tags_masked_full[i],i) );
  }

  MPI_Alltoall(&fof_halo_send_cnt[0],1,MPI_INT,&fof_halo_recv_cnt[0],1,MPI_INT,Partition::getComm());
  MPI_Alltoall(&fof_halo_send_cnt2[0],1,MPI_INT,&fof_halo_recv_cnt2[0],1,MPI_INT,Partition::getComm());

  vector<int> fof_halo_send_off;
  fof_halo_send_off.resize(n_ranks,0);
  vector<int> fof_halo_recv_off;
  fof_halo_recv_off.resize(n_ranks,0);
  vector<int> fof_halo_send_off2;
  fof_halo_send_off2.resize(n_ranks,0);
  vector<int> fof_halo_recv_off2;
  fof_halo_recv_off2.resize(n_ranks,0);


  fof_halo_send_off[0] = fof_halo_recv_off[0] = 0;
  fof_halo_send_off2[0] = fof_halo_recv_off2[0] = 0;

  for (int i=1; i<n_ranks; ++i) {
    fof_halo_send_off[i] = fof_halo_send_off[i-1] + fof_halo_send_cnt[i-1];
    fof_halo_recv_off[i] = fof_halo_recv_off[i-1] + fof_halo_recv_cnt[i-1];
    fof_halo_send_off2[i] = fof_halo_send_off2[i-1] + fof_halo_send_cnt2[i-1];
    fof_halo_recv_off2[i] = fof_halo_recv_off2[i-1] + fof_halo_recv_cnt2[i-1];
  }

  int fof_halo_recv_total = 0;
  int fof_halo_recv_total2 = 0;
  for (int i=0; i<n_ranks; ++i) {
    fof_halo_recv_total += fof_halo_recv_cnt[i];
    fof_halo_recv_total2 += fof_halo_recv_cnt2[i];
  }

  vector<halo_properties_test> fof_halo_recv, fof_halo_recv2;
  fof_halo_recv.resize(fof_halo_recv_total);
  fof_halo_recv2.resize(fof_halo_recv_total2);
  MPI_Barrier(Partition::getComm());


 // send data 
 
  MPI_Alltoallv(&fof_halo_send[0],&fof_halo_send_cnt[0],&fof_halo_send_off[0], H_1.halo_properties_MPI_Type,\
                   &fof_halo_recv[0],&fof_halo_recv_cnt[0],&fof_halo_recv_off[0], H_1.halo_properties_MPI_Type, Partition::getComm());

  MPI_Alltoallv(&fof_halo_send2[0],&fof_halo_send_cnt2[0],&fof_halo_send_off2[0], H_2.halo_properties_MPI_Type,\
                   &fof_halo_recv2[0],&fof_halo_recv_cnt2[0],&fof_halo_recv_off2[0], H_2.halo_properties_MPI_Type, Partition::getComm());


  std::sort(fof_halo_recv.begin(),fof_halo_recv.end(),comp_by_halo_tag);
  std::sort(fof_halo_recv2.begin(),fof_halo_recv2.end(),comp_by_halo_tag);

  H_1.Resize(0);
  H_2.Resize(0);
  for (int i=0; i<fof_halo_recv_total; ++i) {
    halo_properties_test tmp =  fof_halo_recv[i];
    H_1.PushBack(tmp);
  }
  for (int i=0; i<fof_halo_recv_total2; ++i){
     halo_properties_test tmp = fof_halo_recv2[i];
     H_2.PushBack(tmp);
  }
  fof_halo_recv.resize(0);
  fof_halo_recv2.resize(0);


  int err = 0;
  bool skip_err = true;
  vector<int64_t> dup_idx;
  vector<int64_t> dup_idx2;

  int dup_1 = 0;
  int dup_2 = 0;
  for (int i=0; i<H_1.num_halos-1;i++){
    if (H_1.fof_halo_tag->at(i)==H_1.fof_halo_tag->at(i+1)){
        dup_idx.push_back(i);
	dup_idx.push_back(i+1);
    }
  }
    for (int i=0; i<H_2.num_halos-1;i++){
    if (H_2.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i+1)){
        dup_idx2.push_back(i);
        dup_idx2.push_back(i+1);
    }
  }
 
  for (int i=0; i<dup_idx.size();i++){
      H_1.Erase(dup_idx[i]-dup_1);
      dup_1 ++;
  }
    for (int i=0; i<dup_idx2.size();i++){
      H_2.Erase(dup_idx2[i]-dup_2);
      dup_2 ++;
  }
  cout << "duplicates "<< dup_1 << " and "<< dup_2<<endl;
  // remove both sets of duplicates amongst one rank.

  if (H_1.num_halos != H_2.num_halos){
      err += 1;
    }
   else{
    for (int i=0;i<H_1.num_halos;i++){
       if (H_1.fof_halo_tag->at(i)!=H_2.fof_halo_tag->at(i))
         err+=1;
     }
   }
   int numh1 = H_1.num_halos;
   int numh2 = H_2.num_halos;
   int dn_1 = 0;
   int dn_2 = 0;

   // if they don't all match up then run correction
   if (skip_err&&(err>0)){
       err = 0;
       int i=0;
       while ((i < H_1.num_halos)&&(i<H_2.num_halos)){
          if (H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i)){
	     halo_properties_test tmp = H_2.GetProperties(i);
	     H_3.PushBack(tmp);
	     assert(H_1.fof_halo_tag->at(i)==H_3.fof_halo_tag->at(i));
             i += 1;
           }
          else {
           bool not_found = true;
           for (int j=-200;j<200;j++){
	       if (((i+j)<H_2.num_halos)&&(i+j)>=0){
               if (H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i+j)){
		  halo_properties_test tmp = H_2.GetProperties(i+j);
                  H_3.PushBack(tmp);
		  assert(H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i+j));
		  assert(H_1.fof_halo_tag->at(i)==H_3.fof_halo_tag->at(i));

                not_found = false;
                i+=1; // iterate once found
                }
              }
	   }
           if (not_found){
             H_1.Erase(i);
             --numh1;
           }
         }
     }
     while (i<H_1.num_halos){
           bool not_found = true;
           for (int j=-200;j<200;j++){
               if (((i+j)<H_2.num_halos)&&(i+j)>=0){
               if (H_1.fof_halo_tag->at(i)==H_2.fof_halo_tag->at(i+j)){
                  halo_properties_test tmp = H_2.GetProperties(i+j);
                  H_3.PushBack(tmp);
                  assert(H_1.fof_halo_tag->at(i)==H_3.fof_halo_tag->at(i));

                not_found = false;
                i+=1; // iterate once found
                }
              }
           }
           if (not_found){
             H_1.Erase(i);
             --numh1;
           }
         }

   cout << "number of halos"<<endl;
   cout << H_1.num_halos << " " << H_3.num_halos<< endl;

   if (H_1.num_halos>H_3.num_halos){
   int diff_halos = H_1.num_halos-H_3.num_halos;
   for (int kh = 0; kh<diff_halos;kh++){
         H_1.Erase(H_1.num_halos+kh-1);
         --numh1;
   }
   }


   for (int kh = 0 ; kh<H_1.num_halos; kh++){
   	 assert(H_1.fof_halo_tag->at(kh)==H_3.fof_halo_tag->at(kh));

   }

   // look at number of erased halos
   dn_1 =  fof_halo_recv_total - H_1.num_halos;
   dn_2 =  fof_halo_recv_total2 - H_2.num_halos;

  }
  else{
	  if (rank==0){
		  cout << "catalogs line up perfectly "<< endl;
		  cout << "number of halos = " << H_1.num_halos << endl;
	  }
   for (int kh=0; kh<H_1.num_halos; kh++){
          halo_properties_test tmp = H_2.GetProperties(kh);
          H_3.PushBack(tmp);

   }

  }
  int ndiff_tot = 0;
  int ndiff_tot2 = 0;

  int ndups1=0;
  int ndups2 =0;
  MPI_Allreduce(&dn_1, &ndiff_tot, 1, MPI_INT, MPI_SUM,  Partition::getComm());
  MPI_Allreduce(&dn_2, &ndiff_tot2, 1, MPI_INT, MPI_SUM,  Partition::getComm());

  MPI_Allreduce(&dup_1, &ndups1, 1, MPI_INT, MPI_SUM,  Partition::getComm());
  MPI_Allreduce(&dup_2, &ndups2, 1, MPI_INT, MPI_SUM,  Partition::getComm());
  assert(H_1.num_halos == H_3.num_halos);
  err = compute_mean_std_dist_halo(H_1 ,H_3, lim);
  err = compute_mean_std_dist_ellipticity(H_1 ,H_3, lim);

    if ((rank==0)&&(err==0)){
      cout << " Results " << endl;
      cout << " _______ " << endl;
      cout << endl;
      cout << " Comparison test passed! " << endl;
      cout << " All variables within threshold of "  << lim << endl;
      cout << " Total number of non-matching halos = "<< ndiff_tot+ndiff_tot2 << endl;
      cout << " Number of duplicate pairs removed = "<< ndups1 << " in file 1 and "  << ndups2 << " in file 2." << endl;
      cout << endl;
  }
  if ((rank==0)&&(err>0)){
      cout << " Results " << endl;
      cout << " _______ " << endl;
      cout << endl;
      cout << " Comparison exceeded threshold of " << lim << " for " << err << " variables" << endl;
      cout << " out of a total of " <<  N_HALO_FLOATS << " variables " << endl;
      cout << " See above outputs for details  "<< endl;
      cout << " Total number of non-matching halos = "<< ndiff_tot+ndiff_tot2 << endl;
      cout << " Number of duplicate pairs removed = "<< ndups1 << " in file 1 and "  << ndups2 << " in file 2." << endl;
      cout << endl;
  }



  MPI_Barrier(Partition::getComm());
  H_1.Deallocate();
  H_2.Deallocate();
  H_3.Deallocate();
  return 0;

}