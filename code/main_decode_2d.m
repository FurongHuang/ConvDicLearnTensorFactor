% This code implements convolutional tensor decomposition
% copyright Furong Huang, furongh@uci.edu
% Cite paper arXiv:1506.03509
% This function estimates the filters based on conf.sample.

clear;clc;
addpath('fn/')
L = 1;
load(['../data/syntheticData_2d_L',num2str(L),'_estimate.mat']);
estimate.H = zeros(size(conf.sample,2),conf.n*conf.n);
for id_sample = 1 : size(conf.sample,2)
    fprintf('id_sample:%d\n',id_sample);
    filters = estimate.f;
    inv_concated_circulant_filters = cir_inv_2d(filters);
    thisH = inv_concated_circulant_filters*conf.sample(:,id_sample);
    estimate.H(id_sample,:)  = thisH';   
end

