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
%% Default inputs
if nargin == 2
    output_filename = 'AltHap-output.txt';
    thr = 1e-5;
    maxit = 200;
    Sum_Proj = 0;
elseif nargin == 3
    thr = 1e-5;
    maxit = 200;
    Sum_Proj = 0;
elseif nargin == 4
    thr = 1e-5;
    maxit = 200;
end
%% Reading the data
start = 3;
fid = fopen(input_filename);
file = textscan(fid,'%s','delimiter','\n');
content = char(file{1,1});
[read_num, ~] = size(content);
read_num = read_num-2;
line = textscan(content(2,:),'%s','delimiter',' ');
hap_len = str2double(line{1}{1});
R = zeros(read_num,hap_len);
for i = start:read_num+2
    line = textscan(content(i,:),'%s','delimiter',' ');
    line = line{1};
    for j = 1:str2double(line{1})
        R(i - start + 1,str2double(line{j*2+1})+1:str2double(line{j*2+1})+length(line{j*2+2})) = ...
            line{j*2+2} - '0' + 1;
    end
end
alleles_num = length(unique(R))-1;
%% Formulation for Polyploid - Polyallelic

