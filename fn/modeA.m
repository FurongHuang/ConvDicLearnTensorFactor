function [fa_new, lambda_new] = modeA(Tensor,fa_old,fb_old,fc_old,lambda_old,n,L)
firstTerm   = T_fc_odot_fb(Tensor, fc_old, fb_old, n);
secondTerm  = cir_square_dot_inv(fc_old,fb_old,n);
% thirdTerm   = diag_pinv(lambda_old);
T = firstTerm * secondTerm ;% * thirdTerm;
%[fa_new] = circulant_projection(T*thirdTerm,n,L);
[fa_new,lambda_new] = circulant_projection(T,n,L);
[fa_new,fa_new_colwise_norm] = normc(fa_new);
%lambda_new = zeros(n*L,1);
% for l = 1 : L
%     if fa_new_colwise_norm(l)<eps
%         lambda_new(n*(l-1)+1:n*(l-1)+n/2)=0;
%     else
%         lambda_new(n*(l-1)+1:n*(l-1)+n/2)=1/fa_new_colwise_norm(l);
%     end
%     
% end
%fa_new_colwise_norm = fa_new_colwise_norm(:)';
%fa_norm = reshape(repmat(fa_new_colwise_norm,n,1),n*L,1);

%for i = 1 : length(fa_norm)
%    if thirdTerm(i,i)<eps
%        lambda_new(i) = 0;
%    else
%    lambda_new(i) = fa_norm(i)./thirdTerm(i,i);
%    end
%end
%lambda_new = diag(cir_inv(fa_new,n)*T); % this is wrong, as cir_inv(fa_new) is a short fat matrix. It is true if it is a tall thin matrix. 
