# AltHap
By [Abolfazl Hashemi](https://www.linkedin.com/in/abolfazlh/), [Banghua Zhu](https://github.com/13aeon), and [Haris Vikalo](http://users.ece.utexas.edu/~hvikalo/).

### Table of Contents
0. [Introduction](#introduction)
0. [Citation](#citation)
0. [MATLAB implementation](#matlab-implementation)
0. [Input file format](#input-file-format)
0. [Ouput file format](#ouput-file-format)

### Introduction

AltHap is a haplotyping assembly method for diploid, polyploid, and polyallelic polyploid organisms. The method formulates the haplotype assembly problem as a sparse tensor decomposition proplem and exploits the special structures of the factors.
AltHap is implemented in both Python and MATLAB. 

**Note**

The current implementations assume that the input file contains a connected component. The preprocessing step to handle disconnected componenets will be added in near future.

### Citation

If you use AltHap in your research, please cite:

        @article{hashemi2017sparse,
            title={Sparse Tensor Decomposition For Haplotype Assembly Of Diploids And Polyploids},
            author={Hashemi, Abolfazl and Zhu, Banghua and Vikalo, Haris},
            journal={bioRxiv},
            pages={130930},
            year={2017},
            publisher={Cold Spring Harbor Labs Journals}
        }
        
        
### MATLAB implementation

1. Usage in terminal:

```matlab -nojvm -r "AltHap(ploidy,'input_file.txt','output_file.txt');exit"```

The first two inputs are mandatory. If ```output_file``` is not given, the results will be saved in ```AltHap-output.txt```.
Please refer to the source code for a detailed explanation of other options which are mainly related to the stopping criteria.
Please refer to the end of this file for a detailed explanation of input and output files format.

2. Run a simple example:

In this example, ```sim0.txt``` is the input fragment file, ```ploidy = 3```, and the results are saved in ```AltHap-output.txt```. 
There are 3 ways to run this example:

- In terminal write: ``` matlab -nojvm -r "AltHap(3,'sim0.txt');exit" ```

- In terminal write: ``` matlab -nojvm -nodesktop < ./simple_example.m ```

- Run the ```AltHapSimple.sh``` file in terminal. This file would then run the ```simple_example.m``` file.

### Input file format

The input file format is as follows

Number of reads <br />
Number of columns <br />
Number_of_contiguos_segments    Read_identifier    Position_of_first_SNP_segment    Continuous_bases_in_read    Position_of_next_SNP_segment	Continuous_bases_in_read ..... Quality_scores (in fastq format)

- Example for Biallelic:

5568 <br /> 
22801 <br />
2 chr22_SPA9_8733 2 0 5 0 == <br />
1 chr22_SPH2_1940 3 100 C== <br />

- Example for Polyallelic:

4500 <br />
1000 <br />
2	chr3_1	1	13103	37	1132	IIIIIIIII <br />
2	chr3_2	1	1310321	46	110302	IIIIIIIIIIIII <br />
2	chr3_3	1	3220001	44	3001	IIIIIIIIIII <br />
2	chr3_4	1	322000	40	30023	IIIIIIIIIII <br />



### Output file format 

The output file contains MEC score, CPU Time, and Recovered Haplotypes. Subsequently, the phased haplotype is printed in the following format:
first haplotype		second haplotype 	third haplotype 	...

- Example for Biallelic:

MEC: 353 <br />
CPU Time: 38.363051 <br />
Recovered Haplotype: <br />
0 1 1 <br /> 
1 0 1  <br />
0 0 1  <br />
1 1 1 <br /> 
0 1 1  <br />
1 0 0  <br />
1 0 1 <br /> 
1 1 1 <br /> 

- Example for Polyallelic:

MEC: 51 <br />
CPU Time: 28.912201 <br />
Recovered Haplotype: <br />
1 1 1  <br />
1 0 0 <br /> 
3 3 1  <br />
2 2 0  <br />
1 3 1  <br />
3 0 0  <br />
0 0 2  <br />
3 2 0  <br />

For higher ploidy, there are K phased haplotypes instead of 3.
