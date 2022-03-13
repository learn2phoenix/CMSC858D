# include <math.h>
# include <iostream>
# include <fstream>
# include <stdlib.h>
# include <unistd.h>
# include <sdsl/vectors.hpp>
# include <sdsl/int_vector.hpp>

typedef sdsl::bit_vector::size_type size_type;

// Function has been taken from https://www.geeksforgeeks.org/smallest-power-of-2-greater-than-or-equal-to-n/
unsigned long int nextPowerOf2(unsigned long int n)
{
  unsigned long int p = 1;
  if (n && !(n & (n - 1)))
  {
    return n;
  }

  while (p < n)
  {
    p <<= 1;
  }
  return p;
}

sdsl::bit_vector int_to_bit(unsigned int size, unsigned int value){
  sdsl::bit_vector temp(size, 0);
  temp.set_int(0, value);
  return temp;
}


namespace ass1{

class sparse_array{
  public:
    sparse_array():
    ret_var(1){

    }
    int ret_var;

    unsigned long int overhead()
    {
      // We have not included the constant overhead of creating a c++ vector
      return  (super_block.size() * super_block[0].size()) + (block.size() * block[0].size());
    }
    
    void create(unsigned long int x)
    {
      orig_size = x;
      sparse.resize(nextPowerOf2(x));
      sdsl::util::set_to_value(sparse, 0);
    }

    void append(std::string elem, unsigned long int pos)
    {
      dense.push_back(elem);
      sparse[pos] = 1;
      build_rank_dynamic(pos);
    }

    unsigned long int size()
    {
      return sparse.size();
    }

    bool get_at_rank(unsigned long int r, std::string& elem)
    {
      if (r < dense.size())
      {
        elem = dense[r];
        return true;
      }
      return false;
    }

    bool get_at_index(unsigned long int r, std::string& elem)
    {
      if (sparse[r] == 0)
      {
        return false;
      }
      elem = dense[num_elem_at(r)];
      return true;
    }

    unsigned long int num_elem()
    {
      return std::popcount(*sparse.data());
    }

    unsigned long int num_elem_at(unsigned long int ind)
    {
      long int rank;
      rank = super_block[ind / supblock_size].get_int(0) + block[ind / block_size].get_int(0);
      if (rank > dense.size())
      {
        std::cout << ind << " " << rank << " " << dense.size() << " sparse size " << sparse.size() << " " << supblock_size << " " << block_size << std::endl;
        std::cout << sparse << std::endl;
      }
      if (ind % block_size)
      {
        auto first = sparse.begin() + (block_size * (ind / block_size));
        auto last = first + (ind % block_size);
        sdsl::bit_vector temp(ind % block_size, 0);
        copy(first, last, temp.begin());
        rank += std::popcount(*temp.data());
      }

      return rank;
    }
    void save(std::string& fname)
    {
      std::fstream obj;
      obj.open(fname, std::ios::binary | std::ios::trunc | std::ios::out);
      unsigned long int dense_size = dense.size();
      obj.write(reinterpret_cast<char*>(&dense_size), sizeof(dense_size));
      obj.write(reinterpret_cast<char*>(&orig_size), sizeof(orig_size));
      obj.write(reinterpret_cast<char*>(&supblock_width), sizeof(supblock_width));
      obj.write(reinterpret_cast<char*>(&supblock_size), sizeof(supblock_size));
      obj.write(reinterpret_cast<char*>(&block_width), sizeof(block_width));
      obj.write(reinterpret_cast<char*>(&block_size), sizeof(block_size));
      supblock_cnt = super_block.size();
      block_cnt = block.size();
      obj.write(reinterpret_cast<char*>(&supblock_cnt), sizeof(supblock_cnt));
      obj.write(reinterpret_cast<char*>(&block_cnt), sizeof(block_cnt));

      for(int i=0; i < super_block.size(); i++)
      {
        sdsl::serialize(super_block[i], obj);
      }
      for(int i=0; i < block.size(); i++)
      {
        sdsl::serialize(block[i], obj);
      }
      sdsl::serialize(sparse, obj);
      for(int i=0; i < dense_size; i++)
      {
        obj << dense[i] << std::endl;
      }
      obj.close();
    }

    void load(std::string& fname)
    {
      std::ifstream obj;
      obj.open(fname, std::ios::binary | std::ios::in);
      unsigned long int dense_size;
      obj.read((char*)&dense_size, sizeof(dense_size));
      obj.read((char*)&orig_size, sizeof(orig_size));
      obj.read((char*)&supblock_width, sizeof(supblock_width));
      obj.read((char*)&supblock_size, sizeof(supblock_size));
      obj.read((char*)&block_width, sizeof(block_width));
      obj.read((char*)&block_size, sizeof(block_size));
      obj.read((char*)&supblock_cnt, sizeof(supblock_cnt));
      obj.read((char*)&block_cnt, sizeof(block_cnt));
      rank_exists = true;
      // First reading superblocks
      super_block.clear();
      block.clear();
      for(int i=0; i < supblock_cnt; i++)
      {
        sdsl::bit_vector temp(supblock_width, 0);
        sdsl::load(temp, obj);
        super_block.push_back(temp);
      }
      // Now reading blocks
      for(int i=0; i < block_cnt; i++)
      {
        sdsl::bit_vector temp(block_width, 0);
        sdsl::load(temp, obj);
        block.push_back(temp);
      }
      sdsl::util::set_to_value(sparse, 0);
      sdsl::load(sparse, obj);
      dense.clear();
      std::string line;
      while(getline(obj, line)){
        dense.push_back(line);
        // std::cout << line << std::endl;
      }
      obj.close();
    }

  private:
    sdsl::bit_vector sparse;
    std::vector<std::string> dense;
    std::vector<sdsl::bit_vector> super_block;
    std::vector<sdsl::bit_vector> block;
    long int supblock_width;
    long int block_width;
    long int supblock_size;
    long int block_size;
    long int supblock_cnt;
    long int block_cnt;
    long int orig_size;
    bool rank_exists = false;

    void build_rank_dynamic(long int pos)
    {
      if (!rank_exists)
      {
        rank_exists = true;
        long int cnt = 0;
        long int cnt_sub = 0;
        supblock_size = pow(std::log2(sparse.size()), 2);
        supblock_width = ceil(std::log2(orig_size));
        if (!(orig_size % 2)){
          supblock_width++;
        }
        super_block.clear();
        block.clear();
        block_size = std::log2(sparse.size());
        block_width = ceil(std::log2(supblock_size));
        if (!(supblock_size % 2)){
          block_width++;
        }
        for(long int i=0; i < sparse.size(); i++)
        {
          if (!(i % supblock_size))
          {
            super_block.push_back(int_to_bit(supblock_width, cnt));
            cnt_sub = 0;
          }
          if (!(i % block_size))
          {
            block.push_back(int_to_bit(block_width, cnt_sub));
          }
          if (sparse[i]){
            cnt++;
            cnt_sub++;
          }
        }
      }else{
        long int supblock_pos = pos / supblock_size;
        long int block_pos = pos / block_size;
        if (!(supblock_pos == super_block.size() - 1))
        {
          super_block[supblock_pos + 1] = int_to_bit(supblock_width, super_block[supblock_pos + 1].get_int(0) + 1);
        }
        if (!(block_pos == block.size() - 1))
        {
          block[block_pos + 1] = int_to_bit(block_width, block[block_pos + 1].get_int(0) + 1);
        }
      }
      
    } 
};

class select_support{
  public:
    select_support(sdsl::bit_vector* x):
    
    ret_var(1){
      padToPow2(x);
    }
    int ret_var;

    unsigned long int overhead()
    {
      // We have not included the constant overhead of creating a c++ vector
      return  (super_block.size() * super_block[0].size()) + (block.size() * block[0].size());
    }

    unsigned long int select1(long int select)
    {
      long int left = 0;
      long int right = origBit.size() - 1;
      while ((left != right) && (right - left > 1))
      {
        unsigned long int curr_rank =rank1((right + left)/2);
        if (curr_rank >= select)
        {
          right = (right + left)/2;
        }else{
          left = (right + left)/2;
        }
      }
      if (rank1(left) == select)
      {
        return left;
      }else{
        return right;
      }
    }

    unsigned long int rank1(long int ind)
    {
      long int rank;
      rank = super_block[ind / supblock_size].get_int(0) + block[ind / block_size].get_int(0);
      if (ind % block_size)
      {
        auto first = origBit.begin() + (block_size * (ind / block_size));
        auto last = first + (ind % block_size);
        sdsl::bit_vector temp(ind % block_size, 0);
        copy(first, last, temp.begin());
        rank += std::popcount(*temp.data());
      }

      return rank;
    }

    void save(std::string& fname)
    {
      std::fstream obj;
      obj.open(fname, std::ios::binary | std::ios::trunc | std::ios::out);
      obj.write(reinterpret_cast<char*>(&supblock_width), sizeof(supblock_width));
      obj.write(reinterpret_cast<char*>(&supblock_size), sizeof(supblock_size));
      obj.write(reinterpret_cast<char*>(&block_width), sizeof(block_width));
      obj.write(reinterpret_cast<char*>(&block_size), sizeof(block_size));
      supblock_cnt = super_block.size();
      block_cnt = block.size();
      obj.write(reinterpret_cast<char*>(&supblock_cnt), sizeof(supblock_cnt));
      obj.write(reinterpret_cast<char*>(&block_cnt), sizeof(block_cnt));

      for(int i=0; i < super_block.size(); i++)
      {
        sdsl::serialize(super_block[i], obj);
      }
      for(int i=0; i < block.size(); i++)
      {
        sdsl::serialize(block[i], obj);
      }
      sdsl::serialize(origBit, obj);
      obj.close();
    }

    void load(std::string& fname)
    {
      std::ifstream obj;
      obj.open(fname, std::ios::binary | std::ios::in);
      obj.read((char*)&supblock_width, sizeof(supblock_width));
      obj.read((char*)&supblock_size, sizeof(supblock_size));
      obj.read((char*)&block_width, sizeof(block_width));
      obj.read((char*)&block_size, sizeof(block_size));
      obj.read((char*)&supblock_cnt, sizeof(supblock_cnt));
      obj.read((char*)&block_cnt, sizeof(block_cnt));
      // First reading superblocks
      super_block.clear();
      block.clear();
      for(int i=0; i < supblock_cnt; i++)
      {
        sdsl::bit_vector temp(supblock_width, 0);
        sdsl::load(temp, obj);
        super_block.push_back(temp);
      }
      // Now reading blocks
      for(int i=0; i < block_cnt; i++)
      {
        sdsl::bit_vector temp(block_width, 0);
        sdsl::load(temp, obj);
        block.push_back(temp);
      }
      sdsl::util::set_to_value(origBit, 0);
      sdsl::load(origBit, obj);
      obj.close();
    }

  private:
    std::vector<sdsl::bit_vector> super_block;
    std::vector<sdsl::bit_vector> block;
    long int supblock_width;
    long int block_width;
    long int supblock_size;
    long int block_size;
    long int supblock_cnt;
    long int block_cnt;
    sdsl::bit_vector origBit;

    // Function to pad to the nearest pow of 2
    void padToPow2(sdsl::bit_vector* x)
    {
      origBit = *x;
      unsigned long int pow2 = nextPowerOf2((*x).size());
      origBit.resize(pow2);
      long int cnt = 0;
      long int cnt_sub = 0;
      supblock_size = pow(std::log2(pow2), 2);
      supblock_width = ceil(std::log2((*x).size()));
      if (!((*x).size() % 2)){
        supblock_width++;
      }
      block_size = std::log2(pow2);
      block_width = ceil(std::log2(supblock_size));
      if (!(supblock_size % 2)){
        block_width++;
      }
      for(long int i=0; i < pow2; i++)
      {
        if (!(i % supblock_size))
        {
          super_block.push_back(int_to_bit(supblock_width, cnt));
          cnt_sub = 0;
        }
        if (!(i % block_size))
        {
          block.push_back(int_to_bit(block_width, cnt_sub));
        }
        if (origBit[i]){
          cnt++;
          cnt_sub++;
        }
      }
    }  
};

class rank_support {       // The class
  public:             // Access specifier
    rank_support(sdsl::bit_vector* x):
    
    ret_var(1){
      padToPow2(x);
    }
    int ret_var;

    unsigned long int overhead()
    {
      // We have not included the constant overhead of creating a c++ vector
      return  (super_block.size() * super_block[0].size()) + (block.size() * block[0].size());
    }

    unsigned long int rank1(long int ind)
    {
      long int rank;
      rank = super_block[ind / supblock_size].get_int(0) + block[ind / block_size].get_int(0);
      if (ind % block_size)
      {
        auto first = origBit.begin() + (block_size * (ind / block_size));
        auto last = first + (ind % block_size);
        sdsl::bit_vector temp(ind % block_size, 0);
        copy(first, last, temp.begin());
        rank += std::popcount(*temp.data());
      }

      return rank;
    }

    void save(std::string& fname)
    {
      std::fstream obj;
      obj.open(fname, std::ios::binary | std::ios::trunc | std::ios::out);
      obj.write(reinterpret_cast<char*>(&supblock_width), sizeof(supblock_width));
      obj.write(reinterpret_cast<char*>(&supblock_size), sizeof(supblock_size));
      obj.write(reinterpret_cast<char*>(&block_width), sizeof(block_width));
      obj.write(reinterpret_cast<char*>(&block_size), sizeof(block_size));
      supblock_cnt = super_block.size();
      block_cnt = block.size();
      obj.write(reinterpret_cast<char*>(&supblock_cnt), sizeof(supblock_cnt));
      obj.write(reinterpret_cast<char*>(&block_cnt), sizeof(block_cnt));

      for(int i=0; i < super_block.size(); i++)
      {
        sdsl::serialize(super_block[i], obj);
      }
      for(int i=0; i < block.size(); i++)
      {
        sdsl::serialize(block[i], obj);
      }
      sdsl::serialize(origBit, obj);
      obj.close();
    }

    void load(std::string& fname)
    {
      std::ifstream obj;
      obj.open(fname, std::ios::binary | std::ios::in);
      obj.read((char*)&supblock_width, sizeof(supblock_width));
      obj.read((char*)&supblock_size, sizeof(supblock_size));
      obj.read((char*)&block_width, sizeof(block_width));
      obj.read((char*)&block_size, sizeof(block_size));
      obj.read((char*)&supblock_cnt, sizeof(supblock_cnt));
      obj.read((char*)&block_cnt, sizeof(block_cnt));
      // First reading superblocks
      super_block.clear();
      block.clear();
      for(int i=0; i < supblock_cnt; i++)
      {
        sdsl::bit_vector temp(supblock_width, 0);
        sdsl::load(temp, obj);
        super_block.push_back(temp);
      }
      // Now reading blocks
      for(int i=0; i < block_cnt; i++)
      {
        sdsl::bit_vector temp(block_width, 0);
        sdsl::load(temp, obj);
        block.push_back(temp);
      }
      sdsl::util::set_to_value(origBit, 0);
      sdsl::load(origBit, obj);
      obj.close();
    }

  private:
    std::vector<sdsl::bit_vector> super_block;
    std::vector<sdsl::bit_vector> block;
    long int supblock_width;
    long int block_width;
    long int supblock_size;
    long int block_size;
    long int supblock_cnt;
    long int block_cnt;
    sdsl::bit_vector origBit;

    // Function to pad to the nearest pow of 2
    void padToPow2(sdsl::bit_vector* x)
    {
      origBit = *x;
      unsigned long int pow2 = nextPowerOf2((*x).size());
      origBit.resize(pow2);
      long int cnt = 0;
      long int cnt_sub = 0;
      supblock_size = pow(std::log2(pow2), 2);
      supblock_width = ceil(std::log2((*x).size()));
      if (!((*x).size() % 2)){
        supblock_width++;
      }
      block_size = std::log2(pow2);
      block_width = ceil(std::log2(supblock_size));
      if (!(supblock_size % 2)){
        block_width++;
      }
      for(long int i=0; i < pow2; i++)
      {
        if (!(i % supblock_size))
        {
          super_block.push_back(int_to_bit(supblock_width, cnt));
          cnt_sub = 0;
        }
        if (!(i % block_size))
        {
          block.push_back(int_to_bit(block_width, cnt_sub));
        }
        if (origBit[i]){
          cnt++;
          cnt_sub++;
        }
      }
    }   
};
}
