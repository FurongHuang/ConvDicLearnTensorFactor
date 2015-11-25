function [r] = kronecker_prod4vec(a,b)
a = a(:);
b = b(:);
r = a*b';