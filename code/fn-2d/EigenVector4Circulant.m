function [U, U_hermconj] = EigenVector4Circulant(n)
dummy_vec = 0:n-1;
dummy_mat = dummy_vec(:)*dummy_vec(:)';
Psi = exp(complex(0,dummy_mat.*(-2*pi/n)));
U = sqrt(n).* (Psi)^(-1);
U_hermconj = (1/sqrt(n)).*Psi;