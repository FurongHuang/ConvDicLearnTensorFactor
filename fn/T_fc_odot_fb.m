function Result = T_fc_odot_fb(T, fc, fb, n)
% compute T(Fc\odot Fb)
% first compute Fc and Fb
assert(size(fc,1)==size(fb,1));
assert(size(fc,2)==size(fb,2));
a = size(fc,1);
L = size(fc,2);
Result = zeros(n,n*L);
for idx = 1 : n*L
    ido =ceil(idx/n);
    idi = mod(idx-1,n);
    this_fc = zeros(n,1);
    this_fb = zeros(n,1);
    this_fc(1:a)=fc(:,ido);
    this_fb(1:a)=fb(:,ido);
    this_fc = circshift(this_fc,idi);
    this_fb = circshift(this_fb,idi);
    this_r = this_fc*this_fb';   
    Result(:,idx)=T*this_r(:);    
end