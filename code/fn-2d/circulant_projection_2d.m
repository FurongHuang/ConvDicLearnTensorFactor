function [f_new, lambda_new] = circulant_projection_2d(T,n,L)
% T is n^2 \times n^2 \times L
tol = 1e-6;
f_new = zeros(n,n,L);
[T,lambda_new] = normc(T);
for l = 1 : L
    currentT = T(:,(l-1)*n*n+1:l*n*n);    
    
    for q = 1 : n/2 
        currentColumn = zeros(n,n,n/2);
        for j_indx2 = 1 : n/2
            i_indx2 = j_indx2 + (q-1);
            currentColumn(:,:,j_indx2) = currentT((i_indx2-1)*n+1:i_indx2*n,(j_indx2-1)*n+1:j_indx2*n);
        end
        currentColumn_mean = mean(currentColumn,3);
        for p = 1 : n/2
            j_indx = 1 : n/2;
            i_indx = j_indx + (p-1);
            
            current = diag(real(currentColumn_mean(i_indx,j_indx)));
            f_new(p,q,l)= mean(current);
        end
    end
    f_new(n/2+1:n,n/2+1:n,l)=0;
end