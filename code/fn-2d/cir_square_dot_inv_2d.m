function Result = cir_square_dot_inv_2d(fc,fb)
% compute pinv((Fc'*Fc).*(Fb'*Fb))
% for 2d filters, Fc \in R^{n^2 \times n^2L},Fb \in R^{n^2 \times n^2L}
% for 2d filters, fc \in R^{n\times n\times L}
tol = 1e-9;
assert(size(fc,1)==size(fb,1));
assert(size(fc,2)==size(fb,2));
assert(size(fb,1)==size(fb,2));
assert(size(fc,3)==size(fb,3));


n = size(fc,1);
L = size(fc,3);


%% compute block circulant row coefficinets, Construct Blocks of Diagonal Mat
Diag_Block_mat = sparse(zeros(n*n*L,n*n*L));
for id_f1 = 1 : L
    for id_f2 = 1 : L 
        tmp_1 = cconv2_reverse_revised(fc(:,:,id_f1),fc(:,:,id_f2));
        tmp_2 = cconv2_reverse_revised(fb(:,:,id_f1),fb(:,:,id_f2));
        tmp_1dot2 = tmp_1.*tmp_2;
        this_fft2 = fft2(tmp_1dot2);
        current_circulant_coeff_2d = sparse(this_fft2(:));
        current_block_diag = sparse(current_circulant_coeff_2d);
        Diag_Block_mat((id_f1-1)*n*n+1:id_f1*n*n , (id_f2-1)*n*n+1:id_f2*n*n)=diag(current_block_diag);        
    end
end
%% Inverse Blocks of Diagonal Mat
tic
% Diag_Block_mat(Diag_Block_mat<tol) = 0;
% Diag_Block_mat_inv2 = sparse(pinv(full(Diag_Block_mat)));
% fprintf('Diag_Block_mat_inv takes time: %f seconds.\n',toc);
tic
Diag_Block_mat_inv = diag_blk_inv(Diag_Block_mat,L,tol);
fprintf('Diag_Block_mat_inv_smart takes time: %f seconds.\n',toc);
% inv_diff = Diag_Block_mat_inv-Diag_Block_mat_inv2;
%% dummy matrix for U and Psi
dummy_vec = 0:n-1; dummy_mat = dummy_vec(:)*dummy_vec(:)'; Psi = exp(complex(0,dummy_mat.*(-2*pi/n)));
U = n^(-0.5).*Psi;
superU = kronecker_prod4mat(U,U);
P = sparse(zeros(size(Diag_Block_mat)));
for id = 1 : L
    P((id-1)*n*n+1:id*n*n, (id-1)*n*n+1:id*n*n) = superU;
end
Result = real(P * Diag_Block_mat_inv * P');
