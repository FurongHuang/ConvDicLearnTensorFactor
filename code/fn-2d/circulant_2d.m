function BlkCir = circulant_2d(f_mat)
assert(size(f_mat,1)==size(f_mat,2));
n= size(f_mat,1);
% FFT2
F = fft2(f_mat); F = F(:);
% U_super
U = EigenVector4Circulant(n);
U_super = kronecker_prod4mat(U,U);

BlkCir = U_super * diag(F) * U_super';
BlkCir = real(BlkCir);