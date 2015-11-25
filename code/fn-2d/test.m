n=3; L =1;
fb = zeros(n,n,L);
fc = zeros(n,n,L);
Fb = [];Fc=[];
for i = 1:L
    fb(:,:,i)=(magic(n)+i)./norm(magic(n)+i,'fro');
    fc(:,:,i)=(magic(n)+i)./norm(magic(n)+i,'fro');
    Fb=[Fb,circulant_2d(fb(:,:,i))];
    Fc=[Fc,circulant_2d(fb(:,:,i))];
end

%%
Result1 = (Fc'*Fc).*(Fb'*Fb);
%%
Diag_Block_mat = sparse(zeros(n*n*L,n*n*L));
for id_f1 = 1 : L
    for id_f2 = 1 : L 
        tmp_1 = cconv2_reverse(fc(:,:,id_f1),fc(:,:,id_f2));
        tmp_2 = cconv2_reverse(fb(:,:,id_f1),fb(:,:,id_f2));
        tmp_1dot2 = tmp_1.*tmp_2;
        this_fft2 = fft2(tmp_1dot2);
        current_circulant_coeff_2d = sparse(this_fft2(:));
        
        %this_filter_1 = fc(:,:,id_f1).*fb(:,:,id_f1); 
        %this_filter_1_fft2 = fft2(this_filter_1);
        %this_filter_1_fft2=this_filter_1_fft2(:);
        
        %this_filter_2 = fc(:,:,id_f2).*fb(:,:,id_f2); 
        %this_filter_2_fft2 = fft2(this_filter_2);
        %this_filter_2_fft2=this_filter_2_fft2(:);

        %current_circulant_coeff_2d = conj(this_filter_1_fft2).*this_filter_2_fft2;
        current_block_diag = sparse(current_circulant_coeff_2d);
        % assign to the Diag_Block_mat
        Diag_Block_mat((id_f1-1)*n*n+1:id_f1*n*n , (id_f2-1)*n*n+1:id_f2*n*n)=diag(current_block_diag);
        
    end
end

%% dummy matrix for U and Psi
dummy_vec = 0:n-1;
dummy_mat = dummy_vec(:)*dummy_vec(:)';
Psi = exp(complex(0,dummy_mat.*(-2*pi/n)));
U = n^(-0.5).*Psi;
superU = kronecker_prod4mat(U,U);
P = sparse(zeros(size(Diag_Block_mat)));
for id = 1 : L
    P((id-1)*n*n+1:id*n*n, (id-1)*n*n+1:id*n*n) = superU;
end
Result2 = full(real(P * Diag_Block_mat * P'));