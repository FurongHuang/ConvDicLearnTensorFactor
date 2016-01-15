 function Tensor = Construct_Tensor_from_Data(sample,N)
% convert sample into per column per sample format
if size(sample,2) ~= N
    sample = sample';
end
M1 = mean(sample,2);
M2 = sample*sample'.*(1/N);
M3 = sample * matrix_katri_rao(sample,sample)'.*(1/N);

M2_otimes_M1 = zeros(size(sample,1),size(sample,1)^2);
for i = 1: size(sample,1)
    M2_otimes_M1(:,(i-1)*size(sample,1)+1:i*size(sample,1))= M2.*M1(i) + M1*M2(:,i)' + M2(:,i)*M1';
end
M1_otimes_M1_otimes_M1 = M1*katri_rao(M1,M1)';

Tensor  = M3-M2_otimes_M1 + 2.*M1_otimes_M1_otimes_M1;