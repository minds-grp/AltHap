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
if ploidy > 2 && alleles_num > 2
    % Preprocessing
    BS = eye(ploidy);
    Index = double(R ~= 0);
    Index = repmat(Index,1,alleles_num);
    % Encode R to a binary tensor Rt and mode-1 unfolding
    [read_num,hap_len] = size(R);
    Rin = R;
    R = reshape(R,[],1);
    R(R == 0) = alleles_num+1;
    I = [eye(alleles_num);zeros(1,alleles_num)];
    R=reshape(reshape(I(R,:),[read_num,hap_len,alleles_num])...
        ,[read_num,alleles_num*hap_len]);
    [~,hapXallele] = size(R);  
    % Initialization
    [~, ss, Vt] = svds(R,ploidy);
    Vt = Vt*sqrt(ss);
    Vt = Vt';
    % Optimization
    err = inf; err_inV = inf; err_hap = inf;
    err_hist = zeros(1,maxit);
    Vt_last = 100*ones(ploidy,hapXallele);
    iter=0; tic
    while iter < maxit && err > thr && err_inV > thr && err_hap > thr
        iter = iter + 1;
        Uh = zeros(read_num,ploidy);
        for ii = 1:ploidy
            Uh(:,ii) = sum((bsxfun(@minus,R,Vt(ii,:)).*Index).^2,2);
        end
        [~,Uh]=min(Uh,[],2);
        Uh = BS(Uh, :);       
        G_Vt = -Uh'*((R-Uh*Vt).*Index);
        step_size = max(0.09,(0.9*norm(G_Vt,'fro')...
            /norm((Uh*G_Vt).*Index,'fro'))^2);
        Vt = Vt-step_size*G_Vt;
        Vt = max(min(Vt,1),10^(-16));
        if Sum_Proj == 1
            Vt = Vt./repmat(sum(reshape(Vt,[ploidy,hap_len,alleles_num])...
                ,3),[1,alleles_num]);
        end
        % Termination criteria
        err = norm((R-Uh*Vt).*Index,'fro'); err_hist(iter) = err;
        if iter > 1
            err_inV = abs(err_hist(iter) - err_hist(iter-1));
        end
        err_hap = norm(Vt-Vt_last,'fro')/sqrt(hapXallele/ploidy);
        Vt_last = Vt;
    end
    cpu_time = toc;
    Iter_count = iter;
    err_hist = err_hist(1:iter);   
    % Final projection and MEC calculation
    [~,Vt] = max(reshape(Vt,[ploidy,hap_len,alleles_num]),[],3);
    MEC = 0; R = Rin;
    for ii = 1:read_num
        MEC = MEC+min(sum(bsxfun(@ne,bsxfun(@times,Vt,double(R(ii,:)~= 0))...
            ,R(ii,:)),2));
    end
    Vt = Vt'-1;

