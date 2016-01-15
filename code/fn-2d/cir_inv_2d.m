function Result = cir_inv_2d(fa)
% compute pinv(Fa)
L = size(fa,3);

mat = [];
for i = 1: L
    vec = fa(:,:,i);
    mat = [mat,circulant_2d(vec)];
end
Result = pinv(mat);
