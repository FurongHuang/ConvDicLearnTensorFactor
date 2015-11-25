function [ R ] = blkinv2( A, m )
B = mat2cell(full(reshape(nonzeros(A), m, [])), m, ones(size(A,1)/m,1)*m)';
D = cellfun(@(X) sparse(inv(X)), B, 'UniformOutput', false);
R = blkdiag(D{:});
end