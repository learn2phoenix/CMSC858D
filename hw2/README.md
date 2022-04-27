**Requirements**

* Used sdsl library for standard int vector data structure. It can be installed from here: https://github.com/simongog/sdsl-lite
* Used KSeq library to parse FASTA file. It can be installed from here: https://github.com/cartoonist/kseqpp.git
* Used Cereal for serialization. It can be installed from here: https://github.com/USCiLab/cereal.git
* C++2O compatible compiler
* Python3.8 with pandas and matplotlib for plotting the graphs

**Resources Used**

* Random number generator in C++: https://stackoverflow.com/questions/7560114/random-number-c-in-some-range
* LCP search: https://www.techiedelight.com/find-longest-common-prefix-lcp-strings/


**Brief Description**

* There are 2 binaries `querysa` and `buildsa`
* The functionality of normal arguments is as defined in the requirements of the assignment
* In addition to above, corresponding to each of these there is a profile option which is disabled by default. Profiling `querysa` first requires `buildsa` to be run in profile mode. Both these binaries can be run in profile mode by `--profile`
* The profiling code dumps a csv file for each binary which is then used in a IPython notebook to generate graphs. The notebook is also avaiable in the repo.