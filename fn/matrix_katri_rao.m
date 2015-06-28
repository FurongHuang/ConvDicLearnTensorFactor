function C = matrix_katri_rao(A,B)
assert(size(A,2)==size(B,2));
for i =1: size(A,2)
    C(:,i) = katri_rao(A(:,i),B(:,i));
end