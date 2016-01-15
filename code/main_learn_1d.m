% This code implements convolutional tensor decomposition
% copyright Furong Huang, furongh@uci.edu
% Cite paper arXiv:1506.03509 
% This function estimates the filters based on conf.sample. 

clear;clc;
L = 1;
load(['../data/syntheticData_L',num2str(L),'.mat']);
conf.maxIter = 100;
conf.minIter = 1;
conf.tol = 1e-4;
conf.IniTrue = 0;
addpath('fn/');
Tensor = Construct_Tensor_from_Data(conf.sample, conf.N);
estimate = ALS(conf, Tensor);

save(['../data/syntheticData_L',num2str(L),'_estimate.mat'],'conf','estimate');