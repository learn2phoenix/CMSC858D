#include <string>
#include <iostream>
#include <cctype>
#include <fstream>
#include <math.h>
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


sdsl::int_vector<0> my_deserialize(KSeq& record,  std::unordered_map<std::string, interval>& preftab, long int& k, std::string ifile)
{
  std::fstream obj;
  obj.open(ifile, std::ios::binary | std::ios::in);
  cereal::BinaryInputArchive iarchive(obj);
  
  std::string temp_str;
  long int temp_size;
  // Serializing only the sequence part from the record.
  iarchive(record.seq, temp_size, temp_str, preftab, k);

  std::stringstream temp_stream(temp_str);
  sdsl::int_vector<0> temp;
  sdsl::load(temp, temp_stream);
  obj.close();
  return temp;
}


long int LCP(std::string X, std::string Y)
{
    long int i = 0;
    while (i < X.length() && i < Y.length())
    {
        if (X[i] != Y[i]) {
            break;
        }
        i++;
    }
    return i;
}


std::string get_seq_str(sdsl::int_vector<0>& SA, std::string& seq, std::string& query, long int pos)
{
  std::string seq_str;
  if (SA[pos] + query.size() < seq.size())
  {
    seq_str = seq.substr(SA[pos], query.size());
  }else{
    seq_str = seq.substr(SA[pos], seq.size() - SA[pos] - 1);
  }
  return seq_str;
}


void search_naive(sdsl::int_vector<0>& SA, std::string& seq, std::string& query, std::vector<long int>& pos, long int start, long int end, long int lower_limit, long int upper_limit)
{
  long int mid = floor((start + end) /  2);
  std::string seq_str = get_seq_str(SA, seq, query, mid);

  if (SA[mid] + query.size() < seq.size())
  {
    seq_str = seq.substr(SA[mid], query.size());
  }else{
    seq_str = seq.substr(SA[mid], seq.size() - SA[mid] - 1);
  }
  // Now comparing
  if (seq_str == query){
    pos.push_back(SA[mid]);
    if((mid >= lower_limit) && (mid > start))
    {
      search_naive(SA, seq, query, pos, start, mid - 1, lower_limit, upper_limit);
    }
    if((mid <= upper_limit) && (mid < end))
    {
      search_naive(SA, seq, query, pos, mid + 1, end, lower_limit, upper_limit);
    }
  }else{
    if (seq_str > query)
    {
      if((mid >= lower_limit) && (mid > start))
      {
        search_naive(SA, seq, query, pos, start, mid - 1, lower_limit, upper_limit);
      }
    }else{
      if((mid <= upper_limit) && (mid < end))
      {
        search_naive(SA, seq, query, pos, mid + 1, end, lower_limit, upper_limit);
      }
    }
  }
}


void search_accel(sdsl::int_vector<0>& SA, std::string& seq, std::string query, std::vector<long int>& pos, long int start, long int end, long int lower_limit, long int upper_limit, long int start_lcp, 
long int end_lcp)
{
  long int mid = floor((start + end) /  2);
  std::string seq_str;
  long int lcp = std::min(start_lcp, end_lcp);
  std::string curr_query = query.substr(lcp, query.length() - lcp);
  long int new_lcp;
  if (SA[mid] + query.size() < seq.size())
  {
    seq_str = seq.substr(SA[mid] + lcp, curr_query.size());
  }else{
    seq_str = seq.substr(SA[mid] + lcp, seq.size() - SA[mid] - lcp - 1);
  }
  // Now comparing
  if (seq_str == curr_query){
    pos.push_back(SA[mid]);
    if((mid >= lower_limit) && (mid > start))
    {
      new_lcp = std::min(start_lcp, (long int)query.length());
      search_accel(SA, seq, query, pos, start, mid - 1, lower_limit, upper_limit, start_lcp, new_lcp);
    }
    if((mid <= upper_limit) && (mid < end))
    {
      new_lcp = std::min(end_lcp, (long int)query.length());
      search_accel(SA, seq, query, pos, mid + 1, end, lower_limit, upper_limit, new_lcp, end_lcp);
    }
  }else{
    if (seq_str > curr_query)
    {
      if((mid >= lower_limit) && (mid > start))
      {
        new_lcp = std::min(start_lcp, LCP(seq_str, curr_query) + lcp);
        search_accel(SA, seq, query, pos, start, mid - 1, lower_limit, upper_limit, start_lcp, new_lcp);
      }
    }else{
      if((mid <= upper_limit) && (mid < end))
      {
        new_lcp = std::min(end_lcp, LCP(seq_str, curr_query) + lcp);
        search_accel(SA, seq, query, pos, mid + 1, end, lower_limit, upper_limit, new_lcp, end_lcp);
      }
    }
  }
}


void search_queries(sdsl::int_vector<0>& SA, std::string& seq, std::string qfile, std::string mode, std::string ofile,std::unordered_map<std::string, interval> preftab, long int k, bool if_profile)
{
  KSeq record;
  SeqStreamIn iss(qfile.c_str());
  while (iss >> record) {
    if (record.seq.empty()){
        std::cout << "The query is empty. It will not be processed" << std::endl;
    }
    if (if_profile)
    {
      // By default profiling information will be saved in query_profiling.csv
      std::ofstream obj;
      obj.open("./query_profiling.csv");
      obj << "seq_size,query_size,preftab_k,mode,search_time\n";

      break;
    }
    else
    {
      std::ofstream obj;
      obj.open(ofile);
      std::vector<long int> positions;
      if (mode.compare("naive") == 0)
      {
        if (k == 0)
        {
          search_naive(SA, seq, record.seq, positions, 0, SA.size() - 1, 0, SA.size() - 1);
        }
        else{
          // We consult the preftab
          interval search_interval;
          if (preftab.count(record.seq.substr(0, k)) != 0)
          {
            search_interval = preftab[record.seq.substr(0, k)];
            search_naive(SA, seq, record.seq, positions, search_interval.begin, search_interval.end, search_interval.begin, search_interval.end);
          }
        }
      }
      else{
        long int lcp = 0;
        long int start_lcp;
        long int end_lcp;
        if (k == 0)
        {
          start_lcp = LCP(seq.substr(SA[0], record.seq.length()), record.seq);
          end_lcp = LCP(seq.substr(SA[SA.size() - 1], record.seq.length()), record.seq);
          search_accel(SA, seq, record.seq, positions, 0, SA.size() - 1, 0, SA.size() - 1, start_lcp, end_lcp);
        }
        else{
          interval search_interval;
          if (preftab.count(record.seq.substr(0, k)) != 0)
          {
            search_interval = preftab[record.seq.substr(0, k)];
            start_lcp = LCP(seq.substr(SA[search_interval.begin], record.seq.length()), record.seq);
            end_lcp = LCP(seq.substr(SA[search_interval.end], record.seq.length()), record.seq);
            search_accel(SA, seq, record.seq, positions, search_interval.begin, search_interval.end, search_interval.begin, search_interval.end, start_lcp, end_lcp);
          }
        }
      }
      obj << record.name << "\t" << positions.size();
      for(auto i: positions)
      {
        obj << "\t" << i;
      }
      obj << "\n";
      obj.close();
    } 
  }
}


int main(int argc, char *argv[])
{
  std::string ifile;
  std::string qfile;
  std::string outfile;
  std::string mode;
  bool if_profile = 0;
  option readopts[] = {
        {"index", required_argument, NULL, 'i'}, 
        {"queries", required_argument, NULL, 'q'}, 
        {"query_mode", required_argument, NULL, 'm'},
        {"output", required_argument, NULL, 'o'},
        {"profile", optional_argument, NULL, 'e'},{0}};
  while (1) {
        const int opt = getopt_long(argc, argv, "iqmoe:", readopts, 0);
        if (opt == -1) {
            break;
        }
        switch (opt) {
            case 'i':
              ifile = optarg;
              break;
            case 'q':
              qfile = optarg;
              break;
            case 'o':
              outfile = optarg;
              break;
            case 'm':
              mode = optarg;
              break;
            case 'e':
              if_profile = bool(atoi(optarg));
              break;
        }
    }
  if (ifile.empty()){
    std::cout << "You did not specify the --index" << std::endl;
    std::exit (EXIT_FAILURE);
  }
  if (qfile.empty()){
    std::cout << "You did not specify the --queries" << std::endl;
    std::exit (EXIT_FAILURE);
  }
  if (mode.empty()){
    std::cout << "You did not specify the --query_mode" << std::endl;
    std::exit (EXIT_FAILURE);
  }
  if (outfile.empty()){
    std::cout << "You did not specify the --output" << std::endl;
    std::exit (EXIT_FAILURE);
  }
  if ((mode.compare("naive") != 0) && (mode.compare("simpaccel") != 0))
  {
    std::cout << "Accepted options for --query_mode are <naive> and <simpaccel>" << std::endl;
    std::exit (EXIT_FAILURE);
  }
  KSeq record;
  std::unordered_map<std::string, interval> pref_Tab;
  sdsl::int_vector SA;
  long int k;
  SA = my_deserialize(record, pref_Tab, k, ifile);
  search_queries(SA, record.seq, qfile, mode, outfile, pref_Tab, k, if_profile);
}