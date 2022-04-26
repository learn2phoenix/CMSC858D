# include <math.h>
# include <iostream>
# include <fstream>
# include <stdlib.h>
# include <unistd.h>
# include <sdsl/vectors.hpp>
# include <sdsl/int_vector.hpp>
# include <kseq++/seqio.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>


class interval{
  public:
    long int begin;
    long int end;

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(begin, end);
    }
};
