function w = cconv2(f,h)
assert(size(f,1)==size(f,2));
assert(size(h,1)==size(h,2));
assert(size(h,1)==size(f,1));
n= size(f,1);

h_vec = h(:);
U = EigenVector4Circulant(n);
F = fft2(f); F = F(:);
BlkCir = kronecker_prod4mat(U,U) * diag(F) * kronecker_prod4mat(U,U)';
BlkCir = real(BlkCir);
w_vec = BlkCir * h_vec; w = reshape(w_vec,n,n);