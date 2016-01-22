function C = matrix_katri_rao(A,B)
assert(size(A,2)==size(B,2));
C = zeros(size(A,1)*size(B,1),size(A,2));
for i =1: size(A,2)
    C(:,i) = katri_rao(A(:,i),B(:,i));
end