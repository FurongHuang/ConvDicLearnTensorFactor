function estimate = ALS(conf, Tensor)
tic
OriginalTenosrNrom = sqrt(sum(sum(abs(Tensor).^2)));
%% Random Initialization
if conf.IniTrue ==0
    fa_old = orth(randn(conf.n,conf.L));
    fb_old = [orth(randn(conf.a,conf.L));zeros(conf.a,conf.L)];
    fc_old = [orth(randn(conf.a,conf.L));zeros(conf.a,conf.L)];
else
    %% True Initialization
    fa_old = [conf.f];
    fb_old = [conf.f];
    fc_old = [conf.f];
end
lambda_old = ones(conf.n*conf.L,1);
lambda_new = lambda_old;
TensorResidual = Tensor;
%% Iterative Steps
for deflat_id = 1 : conf.L
    fb_old_i = fb_old(:,deflat_id);
    fc_old_i = fc_old(:,deflat_id);
    fa_old_i = fa_old(:,deflat_id);
    lambda_old_i = lambda_old((1-1)*conf.n+1:1*conf.n);
    % f_old_i = [fa_old_i,fb_old_i,fc_old_i];    
    for iter = 1 : conf.maxIter       
        thisTensor = zeros(size(Tensor));
        for i = 0 : conf.n-1
            thisTensor = thisTensor + lambda_old_i(i+1)*circshift(fa_old_i,i)*matrix_katri_rao(circshift(fc_old_i,i),circshift(fb_old_i,i))';            
        end       
        % mode a       
        [fa_new_i, lambda_new_i] = modeA(TensorResidual,fa_old_i,fb_old_i,fc_old_i,lambda_old_i,conf.n,1);
        % mode b
        [fb_new_i, lambda_new_i] = modeA(TensorResidual,fb_old_i,fa_new_i,fc_old_i,lambda_new_i,conf.n,1);     
        % mode c
        [fc_new_i, lambda_new_i] = modeA(TensorResidual,fc_old_i,fa_new_i,fb_new_i,lambda_new_i,conf.n,1);              
        f_new_i = [fa_new_i,fb_new_i,fc_new_i];
        
        % f_old_i = f_new_i;
        fa_old_i = fa_new_i;
        fb_old_i = fb_new_i;
        fc_old_i = fc_new_i;
        lambda_old_i = lambda_new_i;
    end
    TensorResidual = TensorResidual - thisTensor;
    fa_new(:,deflat_id) = fa_new_i;
    fb_new(:,deflat_id) = fb_new_i;
    fc_new(:,deflat_id) = fc_new_i;
    lambda_new((deflat_id-1)*conf.n+1:deflat_id*conf.n) = lambda_new_i;
end
estimate.f = fa_new;
estimate.lambda = lambda_new;
estimate.runningtime = toc;