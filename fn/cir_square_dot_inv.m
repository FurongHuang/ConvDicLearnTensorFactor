function Result = cir_square_dot_inv(fc,fb,n)
% compute pinv((Fc'*Fc).*(Fb'*Fb))
% tol = 1e-6;
assert(size(fc,1)==size(fb,1));
assert(size(fc,2)==size(fb,2));

n = size(fc,1);
L = size(fc,2);
Diag_Block_mat = sparse(zeros(n*L,n*L));
%% dummy matrix for U and Psi
dummy_vec = 0:n-1; 
dummy_mat = dummy_vec(:)*dummy_vec(:)';
Psi = exp(complex(0,dummy_mat.*(-2*pi/n)));
Vbasis = exp(complex(0,dummy_vec.*(-2*pi/n)));
%% compute block circulant row coefficinets, Construct Blocks of Diagonal Mat
%tic
for idx = 1 : L^2
    id_f1 = mod(idx,L)+(mod(idx,L)==0)*L;    
    id_f2 = ceil(idx/L);
    
    this_filter_2 = [fc(:,id_f2)];
    that_filter_2 = [fb(:,id_f2)];
    this_filter_1 = [fc(:,id_f1)];        
    that_filter_1 = [fb(:,id_f1)]; 
    
    current_circulant_coeff = circulant_row2col(multip_circulant_col_coeff(this_filter_1,this_filter_2)...
            .*multip_circulant_col_coeff(that_filter_1,that_filter_2));        
        current_block_diag = sparse(Psi * current_circulant_coeff(:));
        % assign to the Diag_Block_mat
        Diag_Block_mat((id_f1-1)*n+1:id_f1*n , (id_f2-1)*n+1:id_f2*n)=diag(current_block_diag);
    
end
%fprintf('Diag_Block_mat Construction takes time: %f seconds.\n',toc);
%% Inverse Blocks of Diagonal Mat
%tic
% Diag_Block_mat(Diag_Block_mat<tol) = 0;
Diag_Block_mat_inv = pinv(full(Diag_Block_mat));
% Diag_Block_mat_inv = diag_blk_inv(Diag_Block_mat,L);
% fprintf('Diag_Block_mat_inv takes time: %f seconds.\n',toc);

U = n^(-0.5).*Psi;
P = sparse(zeros(size(Diag_Block_mat)));
for id = 1 : L   
    P((id-1)*n+1:id*n, (id-1)*n+1:id*n) = U;
end
Result = real(P * Diag_Block_mat_inv * P');
