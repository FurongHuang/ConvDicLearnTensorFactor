function [f_new, lambda_new] = circulant_projection(T,n,L)
tol = 1e-6;
f_new = zeros(n,L);
[T,lambda_new] = normc(T);
for l = 1 : L
    lambda_new(n*(l-1)+n/2+1:n*l)=0;
end
for l = 1 : L
    currentT = T(:,(l-1)*n+1:l*n);
    for p = 1 : n/2
        j_indx = 1: n/2;
        i_indx = j_indx + (p-1);
        % i_indx = circshift(j_indx',-(p-1))';
        current = diag(real(currentT(i_indx,j_indx)));
        f_new(p,l)= mean(current);
    end
    f_new(n/2+1:n,l)=0;
%    thisnorm = norm(f_new(:,l));
%     if thisnorm < tol
%         %f_new(:,l) = circulant_row2col(f_new(:,l));
%         f_new(:,l) = f_new(:,l);
%     else
%         %f_new(:,l) = circulant_row2col(f_new(:,l));
%         f_new(:,l) = f_new(:,l);
%     end      
end