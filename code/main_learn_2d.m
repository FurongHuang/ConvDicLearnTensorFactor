% This code implements convolutional tensor decomposition
% copyright Furong Huang, furongh@uci.edu
% Cite paper arXiv:1506.03509 
% This function estimates the filters based on conf.sample. 

clear;clc;
L = 1;
load(['../data/syntheticData_2d_L',num2str(L),'.mat']);
conf.maxIter = 100;
conf.minIter = 1;
conf.tol = 1e-4;
conf.IniTrue = 1;
addpath('fn-2d/');
Tensor = Construct_Tensor_from_Data(conf.sample, conf.N);
% Tensor = circulant_2d(conf.f)*diag(ones(conf.n*conf.n*conf.L,1)*conf.lambda)*(matrix_katri_rao(circulant_2d(conf.f),circulant_2d(conf.f)))';

estimate = ALS_2d(conf, Tensor);

save(['../data/syntheticData_2d_L',num2str(L),'_estimate.mat'],'conf','estimate');