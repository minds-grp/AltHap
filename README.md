# AltHap
AltHap is a haplotyping assembly method for diploid, polyploid, and polyallelic polyploid organisms. The method formulates the haplotype assembly problem as a tensor factorization proplem and exploits the special structures of the factors.

AltHap is implemented in Python and MATLAB. 
Note: The current implementations assume that the input file contains a connected component. The preprocessing step to handle disconnected componenets will be added in near future.


------------------------------------------------------------------------------------------------
MATLAB implementation

1. Usage in terminal:

matlab -nojvm -r "AltHap(ploidy,'input_file.txt','output_file.txt');exit"

The first two inputs are mandatory. If output_file is not given, the results will be saved in 'AltHap-output.txt'.
Please refer to the source code for a detailed explanation of other options which are mainly related to the stopping criteria.
Please refer to the end of this file for a detailed explanation of input and output files format.



2. Run a simple example:

In this example, 'sim0.txt' is the input fragment file, ploidy = 3, and the results are saved in 'AltHap-output.txt'. 
There are 3 ways to run this example:

- In terminal write: matlab -nojvm -r "AltHap(3,'sim0.txt');exit"

- In terminal write: matlab -nojvm -nodesktop < ./simple_example.m

- Run the 'AltHapSimple.sh' file in terminal. This file would then run the 'simple_example.m' file.