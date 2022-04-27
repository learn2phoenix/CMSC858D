#include <iostream>
#include <fstream>
#include <getopt.h>
#include <kseq++/seqio.hpp>
#include <ass2.hpp>
#include <algorithm>

using namespace klibpp;


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



int main(int argc, char *argv[]){
    std::string filename;
    std::string outfile;
    long int preftab = 0;
    bool if_profile = 0;
    option readopts[] = {
        {"reference", required_argument, NULL, 'r'}, 
        {"output", required_argument, NULL, 'o'}, {0}};
    while (1) {
        const int opt = getopt_long(argc, argv, "ro:", readopts, 0);
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
        }
    }

    KSeq record = read_seq(filename.c_str());
    SeqStreamOut oss("./CMSC858D_S22_Project2_sample/human_chr20_mod.fa");
    oss << record;


    std::random_device rd;
    std::mt19937 gen(rd());
    std::ofstream obj;
    obj.open(outfile);

    int i = 0;
    std::uniform_int_distribution<> start_picker(0, record.seq.size() - 61);
    std::uniform_int_distribution<> length_picker(10, 60);
    std::uniform_int_distribution<> fudger_picker(0, 9);
    long int start_pos;
    int length_seq ;
    int fudger;
    while(i < 10000)
    {
        start_pos = start_picker(gen);
        length_seq = length_picker(gen);
        obj << ">"<< i << ":" << start_pos << ":";

        fudger = fudger_picker(gen);
        if (fudger == 9)
        {
            obj << "M\n"; 
            std::string temp = record.seq.substr(start_pos, length_seq);
            std::random_shuffle(temp.begin(), temp.end());
            obj << temp << "\n";
        }
        else{
            obj << "R\n";
            obj << record.seq.substr(start_pos, length_seq) << "\n";
        }
        i++;
    }

}