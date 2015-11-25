function Result = T_fc_odot_fb_2d(T, fc, fb, n)
% compute T(Fc\odot Fb)
% first compute Fc and Fb
% for 2d filters, Fc \in R^{n^2 \times n^2L},Fb \in R^{n^2 \times n^2L}
% for 2d filters, fc \in R^{n\times n\times L}
assert(size(fc,1)==size(fb,1));
assert(size(fc,2)==size(fb,2));
a = size(fc,1);
L = size(fc,2);
Result = zeros(n*n,n*n*L);
for i = 1 : L
    this_fc = zeros(n,n); this_fc(1:a,1:a)=fc(:,:,i);
    this_fb = zeros(n,n); this_fb(1:a,1:a)=fb(:,:,i);
    this_r = katri_rao4mat(circulant_2d(this_fc),circulant_2d(this_fb));
    Result(:,n*n*(i-1)+1:n*n*i) = T * this_r;
end
% this_r = this_fc*this_fb';
% Result(:,idx)=T*this_r(:);