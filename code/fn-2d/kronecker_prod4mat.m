function [R] = kronecker_prod4mat(A,B)
[rowsA,colsA] = size(A);
[rowsB,colsB] = size(B);
R = zeros(rowsA*rowsB, colsA*colsB);
for id_B = 1 : rowsB
    for jd_B = 1 : colsB
        R(rowsA*(id_B-1)+1:rowsA*id_B, colsA*(jd_B-1)+1:colsA*jd_B) = A*B(id_B,jd_B);
    end
end