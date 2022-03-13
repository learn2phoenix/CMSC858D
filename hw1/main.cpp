# include <chrono>
# include <iostream>
# include <fstream>
# include <stdlib.h>
# include <bits/stdc++.h>
# include <sdsl/vectors.hpp>
# include <sdsl/util.hpp>
# include <sdsl/int_vector.hpp>
# include <sdsl/suffix_arrays.hpp>
# include <ass1.hpp>

using namespace std;
using namespace sdsl;
using namespace ass1;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::milliseconds;


long int gen_random(long int start, long int end)
{
    std::random_device seed;
    std::mt19937 gen{seed()}; // seed the generator
    std::uniform_int_distribution dist{start, end}; // set min and max
    long int guess = dist(gen); 
    return guess;
}

// This function has been copied from https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
string gen_random_str(const int len) {
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    string tmp_s;
    tmp_s.reserve(len);

    for (int i = 0; i < len; ++i) {
        tmp_s += alphanum[rand() % (sizeof(alphanum) - 1)];
    }
    
    return tmp_s;
}

// This function has been copied from https://stackoverflow.com/questions/28574346/find-average-of-input-to-vector-c
double average(std::vector<double> const& v){
    if(v.empty()){
        return 0;
    }

    auto const count = static_cast<double>(v.size());
    return std::reduce(v.begin(), v.end()) / count;
}


void test_sparse_array()
{
    std::vector<long int> sizes = {1000, 10000, 100000, 1000000};
    std::vector<float> sparsity = {.01, .05, .1, .2, .35, .5, .75, 1.0};
    std::vector<double> create_times;
    std::vector<double> avg_append_times;
    std::vector<double> get_at_rank_times;
    std::vector<double> get_at_index_times;
    std::vector<double> num_elem_at_times;
    std::vector<double> num_elem_times;
    std::vector<double> save_times;
    std::vector<double> load_times;
    std::vector<double> space_savings;
    duration<double, std::milli> ms_double;
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    string temp;
    long int temp_pos;
    long int return_pos;
    string fname = "test.txt";
    long int empty_size;
    long int sparse_size;

    for(int i = 0; i < sizes.size(); i++)
    {
        for(int j = 0; j < sparsity.size(); j++)
        {
            ass1::sparse_array sp;

            t1 = high_resolution_clock::now();
            sp.create(sizes[i]);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            create_times.push_back(ms_double.count());

            long int elem_cnt = sizes[i] * sparsity[j];
            long int spacing = sizes[i]/ elem_cnt;
            long int curr_pos = 0;

            cout << "At size: " << sizes[i] << " Sparsity: "<< sparsity[j] << " Spacing " << spacing << endl;
            std::vector<double> append_times;
            for(int cnt = 0; cnt < elem_cnt; cnt++)
            {   
                t1 = high_resolution_clock::now();
                sp.append(gen_random_str(10), curr_pos);
                t2 = high_resolution_clock::now();
                ms_double = t2 - t1;
                append_times.push_back(ms_double.count());
                curr_pos += spacing;
            }
            // cout << "Average append time is " << average(append_times) << " Length of append is " << append_times.size() << endl;
            avg_append_times.push_back(average(append_times));
            // Testing get_at_rank
            t1 = high_resolution_clock::now();
            sp.get_at_rank(gen_random(0, elem_cnt), temp);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            get_at_rank_times.push_back(ms_double.count());

            // Testing get_at_index
            t1 = high_resolution_clock::now();
            temp_pos = gen_random(0, sizes[i] - 1);
            sp.get_at_index(temp_pos, temp);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            get_at_index_times.push_back(ms_double.count());

            // Testing num_elem_at
            t1 = high_resolution_clock::now();
            temp_pos = gen_random(0, sizes[i] - 1);
            return_pos = sp.num_elem_at(temp_pos);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            num_elem_at_times.push_back(ms_double.count());

            // Testing num_elem
            t1 = high_resolution_clock::now();
            return_pos = sp.num_elem();
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            num_elem_times.push_back(ms_double.count());

            // Testing save
            t1 = high_resolution_clock::now();
            sp.save(fname);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            save_times.push_back(ms_double.count());

            // Testing load
            t1 = high_resolution_clock::now();
            sp.load(fname);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            load_times.push_back(ms_double.count());

            // Testing space savings. We have ignored the actual data stored in the the indices since it will be same for both cases. Insted we would just be comparing 
            // the potential saving between initializing an empty string vectorr vs that of our sparse implementation
            empty_size = sizeof(string) * sizes[i];
            sparse_size =  (sizeof(string) * sizes[i] * sparsity[j]) + sp.overhead();
            space_savings.push_back((float) sparse_size / (float) empty_size);

        }
    }
    // Saving these numbers in a text file beause I have no idea how to do plots in C++. Will use python for plots
    ofstream myfile;
    myfile.open ("test_sparse.csv");
    myfile << "Size,Sparsity,Create,Append,get_at_rank,get_at_index,num_elem_at,num_elem,save,load, space\n";
    long int cnt = 0;
    for(int i = 0; i < sizes.size(); i++)
    {
        for(int j = 0; j < sparsity.size(); j++)
        {
            myfile << sizes[i] << "," << sparsity[j] << "," << create_times[cnt] << "," << avg_append_times[cnt] << "," << get_at_rank_times[cnt] << "," << get_at_index_times[cnt]
            << "," << num_elem_at_times[cnt] << "," << num_elem_times[cnt] << "," << save_times[cnt] << "," << load_times[cnt] << "," << space_savings[cnt] << "\n";
            cnt++;
        }
    }
    myfile.close();
}

void test_select_support()
{
    std::vector<long int> sizes = {1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000, 50000000, 100000000};
    std::vector<double> avg_query_times;
    std::vector<long int> overheads;
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double;
    long int max_rank;
    long int pos;
    for(int i = 0; i < sizes.size(); i++)
    {
        cout << "At size " << sizes[i] << endl;
        size_type n = sizes[i];
        bit_vector temp(n, 0);
        for(int j=0; j < temp.size(); j++)
        {
            temp[j] = rand()%2;
        }
        ass1::select_support sel_sup(&temp);
        // Measuring a fixed number (1000) of select operations
        max_rank = sel_sup.rank1(sizes[i] - 1);
        std::vector<double> query_times;
        for(int j =0; j < 1000; j++)
        {
            pos = gen_random(0, max_rank);
            t1 = high_resolution_clock::now();
            sel_sup.select1(pos);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            query_times.push_back(ms_double.count());
        }
        avg_query_times.push_back(average(query_times));
        // Getting the overhead of ds
        overheads.push_back(sel_sup.overhead());
    }

    // Saving these numbers in a text file beause I have no idea how to do plots in C++. Will use python for plots
    ofstream myfile;
    myfile.open ("test_select.csv");
    myfile << "Size,query,overhead\n";
    for(int i = 0; i < sizes.size(); i++)
    {
        myfile << sizes[i] << "," << avg_query_times[i] << "," << overheads[i] << "\n";
    }
    myfile.close();
}

void test_rank_support()
{
    std::vector<long int> sizes = {1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000, 50000000, 100000000};
    std::vector<double> avg_query_times;
    std::vector<long int> overheads;
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double;
    long int pos;
    for(int i = 0; i < sizes.size(); i++)
    {
        cout << "At size " << sizes[i] << endl;
        size_type n = sizes[i];
        bit_vector temp(n, 0);
        for(int j=0; j < temp.size(); j++)
        {
            temp[j] = rand()%2;
        }
        ass1::rank_support rank_sup(&temp);
        // Measuring a fixed number (1000) of select operations
        std::vector<double> query_times;
        for(int j =0; j < 1000; j++)
        {
            pos = gen_random(0, sizes[i] - 1);
            t1 = high_resolution_clock::now();
            rank_sup.rank1(pos);
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1;
            query_times.push_back(ms_double.count());
        }
        avg_query_times.push_back(average(query_times));
        // Getting the overhead of ds
        overheads.push_back(rank_sup.overhead());
    }

    // Saving these numbers in a text file beause I have no idea how to do plots in C++. Will use python for plots
    ofstream myfile;
    myfile.open ("test_rank.csv");
    myfile << "Size,query,overhead\n";
    for(int i = 0; i < sizes.size(); i++)
    {
        myfile << sizes[i] << "," << avg_query_times[i] << "," << overheads[i] << "\n";
    }
    myfile.close();
}

int main(){
    // test_rank_support();
    // test_select_support();
    test_sparse_array();
}