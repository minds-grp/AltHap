function [Vt,MEC,cpu_time,err_hist,Iter_count] = AltHap(ploidy,input_filename,...
    output_filename,Sum_Proj,thr,maxit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstructing Haplotype sequences of a single Individual
% Using Alternating Hplotyping algorithm
% Objective function : Frobenius norm of error
%
% Input parameters
% ploidy	      : Number of haplotype sequences (e.g., 2 for humans)
% input_filename  : Name of SNP fragment matrix file
% output_filename : Name of output file
% Sum_Proj        : set 1 to normalize entries of Vt when updating
% thr             : Stopping threshold
% maxit           : Maximum number of iteration
% Output parameters
% Vt              : Estimate of Haplotype sequence matrix
% MEC             : Minimum Error Correction
% cpu_time        : CPU time
% err_hist        : Histoty of objective value
% Iter_count      : Maximum iteration reached by algorithm
%
% Written by Abolfazl Hashemi and Banghua Zhu, September 2016
% ECE department, UT Austin, Austin, TX, 78712, US
% Email: abolfazl@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
