function Result = cir_inv(fa)
% compute pinv(Fa)
L = size(fa,2);

mat = [];
for i = 1: L
    vec = fa(:,i);
    mat = [mat,circulant(vec)];
end
Result = pinv(mat);
