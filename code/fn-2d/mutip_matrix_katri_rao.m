function  C = mutip_matrix_katri_rao(sample)
% sample * matrix_katri_rao(sample,sample)'.*(1/N);
n=size(sample,1);
N=size(sample,2);
C = zeros(n,n*n);
for idx_sample = 1 : size(sample,2)
    this_sample = sample(:,idx_sample);
    C = C + this_sample*katri_rao(this_sample,this_sample)';
end
C= C.*(1/N);