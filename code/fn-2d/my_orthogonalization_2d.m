function output=my_orthogonalization_2d(f)

assert(size(f,1)==size(f,2));
n = size(f,1);
L = size(f,3);
flattened_f = zeros(n*n,L);
for id = 1 : L
    flattened_f(:,id)=reshape(f(:,:,id),n*n,1);
    assert(norm(flattened_f(:,id))==1);
end
orthogonalized = qr(flattened_f);

output = zeros(n,n,L);
for id = 1 : L
    output(:,:,id)=reshape(orthogonalized(:,id),n,n);
end
