function [fa_new, lambda_new] = modeA_2d(Tensor,fa_old,fb_old,fc_old,lambda_old,n,L)
firstTerm   = T_fc_odot_fb_2d(Tensor, fc_old, fb_old, n);
secondTerm  = cir_square_dot_inv_2d(fc_old,fb_old);
T = firstTerm * secondTerm; % * thirdTerm;
[fa_new,lambda_new] = circulant_projection_2d(T,n,L);
[fa_new,fa_new_colwise_norm] = normc_2d(fa_new);
% TODO: modifi cir_squre_dot_inv_2d
% TODO: 2d_circulant project
