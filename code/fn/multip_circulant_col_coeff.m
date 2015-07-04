function gamma = multip_circulant_col_coeff(alpha,beta)
% alpha is the first column of circulant A
% beta is the first column of circulant B
% gamma is the first column of result circulant
assert(length(alpha)==length(beta));
beta = beta(:);
alpha_coeff=alpha(:);
gamma = zeros(length(alpha),1);

for i =0: length(alpha)-1
    beta_coeff_part = beta(i+1:end);
    beta_coeff_part2 = beta(1:i); 
    gamma(i+1) = alpha_coeff'*[beta_coeff_part;beta_coeff_part2];
end