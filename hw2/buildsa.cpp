#include <string>
#include <filesystem>
#include <iostream>
#include <cctype>
#include <fstream>
#include <cereal/archives/json.hpp>
#include <cereal/archives/adapters.hpp>
#include <cereal/archives/binary.hpp>
#include <getopt.h>
#include <kseq++/seqio.hpp>
#include <ass2.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct_sa.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <unordered_map>

using namespace klibpp;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::milliseconds;

// Converts lowercase to uppercase and psedo-randomly converts anything other than ATCG to ATCG
void correct_seq(KSeq& record)
{
  std::string valids = "ATCG";
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distr(0, 3);
  for (long int i=0; i < record.seq.size(); i++)
  {
    record.seq[i] = toupper(record.seq[i]);
    if (valids.find(record.seq[i]) == std::string::npos)
    {
      record.seq[i] = valids[distr(gen)];
    }
  }
}


KSeq read_seq(const char* filename)
{
  KSeq record;
  SeqStreamIn iss(filename);
  while (iss >> record) {
    break;  // Assuming that there is only one sequence in the FASTA file
  }
  if (record.seq.empty()){
    std::cout << "The specified file doesn't have a sequence" << std::endl;
    std::exit (EXIT_FAILURE);
  }
  correct_seq(record);
  return record;
}


// template<uint8_t fixedIntWidth>
void my_serialize(KSeq& record, sdsl::int_vector<0>& temp,  std::unordered_map<std::string, interval>& preftab, long int k, std::string ofile)
{
  std::fstream obj;
  obj.open(ofile, std::ios::binary | std::ios::trunc | std::ios::out);
  cereal::BinaryOutputArchive oarchive(obj);
  
  std::stringstream temp_stream;
  sdsl::serialize(temp, temp_stream);
  long int temp_size = temp.size();
  // Serializing only the sequence part from the record.
  oarchive(record.seq, temp_size, temp_stream.str(), preftab, k);
  obj.close();
}


void build_preftab(sdsl::int_vector<0>& SA, std::string& seq,  std::unordered_map<std::string, interval>& pref_Tab, long int k){
  long int curr_pos;
  std::string prefix;
  std::string curr_prefix;
  interval rng;
  for(long int i=0; i < SA.size(); i++)
  {
    if (SA[i] + k < seq.size())
    {
      if (!prefix.empty())
      {
        curr_prefix = seq.substr(SA[i], k);
        if (prefix.compare(curr_prefix) != 0)
        {
          rng.begin = curr_pos;
          rng.end = i;
          pref_Tab[prefix] = rng;
          prefix = curr_prefix;
          curr_pos = i;
        } 
      }else{
        prefix = seq.substr(SA[i], k);
        curr_pos = i;
      }

    }else{
      if (!prefix.empty()){
        rng.begin = curr_pos;
        rng.end = i;
        pref_Tab[prefix] = rng;
        prefix = "";
      }
    }
  }
  if (!prefix.empty()){
    rng.begin = curr_pos;
    rng.end = SA.size() - 1; 
    pref_Tab[prefix] = rng;
    prefix = "";
  }
  // for( const auto& [key, value] : pref_Tab ) {
  //   std::cout << key <<":" << value.begin << "," << value.end << std::endl;
  // }
}


int main(int argc, char *argv[])
{
  std::string filename;
  std::string outfile;
  long int preftab = 0;
  bool if_profile = 0;
  option readopts[] = {
        {"reference", required_argument, NULL, 'r'}, 
        {"output", required_argument, NULL, 'o'}, 
        {"preftab", optional_argument, NULL, 'p'},
        {"profile", optional_argument, NULL, 'e'}, {0}};
  while (1) {
        const int opt = getopt_long(argc, argv, "proe:", readopts, 0);
        if (opt == -1) {
            break;
        }
        switch (opt) {
            case 'r':
              filename = optarg;
              break;
            case 'o':
              outfile = optarg;
              break;
            case 'p':
              preftab = atoi(optarg);
              break;
            case 'e':
              if_profile = bool(atoi(optarg));
              break;
        }
    }
    if (filename.empty()){
      std::cout << "You did not specify the --reference" << std::endl;
      std::exit (EXIT_FAILURE);
    }
    if (outfile.empty()){
      std::cout << "You did not specify the --output" << std::endl;
      std::exit (EXIT_FAILURE);
    }
  
  if (if_profile)
  {
    // By default profiling information will be saved in create_profiling.csv
    std::ofstream obj;
    obj.open("./create_profiling.csv");

    obj << "seq_size,preftab_k,serialized_size,create_time\n";

    duration<double, std::milli> ms_double;
  
    // profiling on smaller file. DEafult location expected is CMSC858D_S22_Project2_sample/ecoli.fa
    filename = "./CMSC858D_S22_Project2_sample/ecoli.fa";
    KSeq record = read_seq(filename.c_str());
    int SA_width = ceil(std::log2(record.seq.length()));
    int k = 0;
    while (k <= 15)
    {
      sdsl::int_vector SA(record.seq.length(), 0);
      SA.width(SA_width);
      auto t1 = high_resolution_clock::now();
      sdsl::algorithm::calculate_sa((const unsigned char*) record.seq.c_str(), record.seq.length(), SA);
      std::unordered_map<std::string, interval> pref_Tab;
      if (k != 0){
        build_preftab(SA, record.seq, pref_Tab, k);
      }
      auto t2 = high_resolution_clock::now();
      ms_double = t2 - t1;
      std::string index_file =  "indexes/ecoli_" + std::to_string(k);
      my_serialize(record, SA, pref_Tab, k, index_file);
      std::filesystem::path p(index_file);
      obj << record.seq.size() << "," << k << "," << std::filesystem::file_size(p) << "," << ms_double.count() << "\n";
      std::cout << "SA creation time at i=" << filename << ",k=" << k << ":" << ms_double.count() << std::endl;
      if (k == 0)
      {
        k = 1;
      }
      else{
        k += 2;
      }
    }

    // profiling on larger file. DEafult location expected is CMSC858D_S22_Project2_sample/human_chr20.fa
    filename = "./CMSC858D_S22_Project2_sample/human_chr20_mod.fa";
    record = read_seq(filename.c_str());
    SA_width = ceil(std::log2(record.seq.length()));
    k = 0;
    while (k <= 15)
    {
      sdsl::int_vector SA(record.seq.length(), 0);
      SA.width(SA_width);
      auto t1 = high_resolution_clock::now();
      sdsl::algorithm::calculate_sa((const unsigned char*) record.seq.c_str(), record.seq.length(), SA);
      std::unordered_map<std::string, interval> pref_Tab;
      if (k != 0){
        build_preftab(SA, record.seq, pref_Tab, k);
      }
      auto t2 = high_resolution_clock::now();
      ms_double = t2 - t1;
      std::string index_file =  "indexes/human_chr20_" + std::to_string(k);
      my_serialize(record, SA, pref_Tab, k, index_file);
      std::filesystem::path p(index_file);
      obj << record.seq.size() << "," << k << "," << std::filesystem::file_size(p) << "," << ms_double.count() << "\n";
      std::cout << "SA creation time at i=" << filename << ",k=" << k << ":" << ms_double.count() << std::endl;
      if (k == 0)
      {
        k = 1;
      }
      else{
        k += 2;
      }
    }
    obj.close();
  }
  else
  {
    KSeq record = read_seq(filename.c_str());
    int SA_width = ceil(std::log2(record.seq.length()));
    sdsl::int_vector SA(record.seq.length(), 0);
    SA.width(SA_width);
    std::unordered_map<std::string, interval> pref_Tab;
    sdsl::algorithm::calculate_sa((const unsigned char*) record.seq.c_str(), record.seq.length(), SA);
    if (preftab != 0)
    {
      build_preftab(SA, record.seq, pref_Tab, preftab);
    }
    my_serialize(record, SA, pref_Tab, preftab, outfile);
  }
}