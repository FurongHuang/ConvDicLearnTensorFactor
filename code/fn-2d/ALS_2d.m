function estimate = ALS_2d(conf, Tensor)
tic
OriginalTenosrNrom = sqrt(sum(sum(abs(Tensor).^2)));
%% Random Initialization
if conf.IniTrue ==0
    fa_old = zeros(conf.n,conf.n,conf.L);
    fb_old = zeros(conf.n,conf.n,conf.L);
    fc_old = zeros(conf.n,conf.n,conf.L);
    for ind_L = 1 : conf.L
        fa_old(1:conf.a,1:conf.a,ind_L) = orth(randn(conf.a,conf.a));
        fb_old(1:conf.a,1:conf.a,ind_L) = orth(randn(conf.a,conf.a));
        fc_old(1:conf.a,1:conf.a,ind_L) = orth(randn(conf.a,conf.a));
    end
    lambda_old = ones(conf.n*conf.n*conf.L,1);
    lambda_new = lambda_old;
else
    %% True Initialization
    fa_old = [conf.f];
    fb_old = [conf.f];
    fc_old = [conf.f];
    lambda_old = conf.lambda*ones(conf.n*conf.n*conf.L,1);
    lambda_new = lambda_old;
end

TensorResidual = Tensor;
%% Iterative Steps
for deflat_id = 1 : conf.L
    fb_old_i = fb_old(:,:,deflat_id);
    fc_old_i = fc_old(:,:,deflat_id);
    fa_old_i = fa_old(:,:,deflat_id);
    lambda_old_i = lambda_old((1-1)*conf.n*conf.n+1:1*conf.n*conf.n);   
    for iter = 1 : conf.maxIter       
        % thisTensor = zeros(size(Tensor));
        thisTensor = circulant_2d(fa_old_i)*diag(lambda_old_i)*(matrix_katri_rao(circulant_2d(fc_old_i),circulant_2d(fb_old_i)))';
        assert(size(thisTensor,1)==size(Tensor,1));
        assert(size(thisTensor,2)==size(Tensor,2));
        assert(size(thisTensor,3)==size(Tensor,3));
        
        % mode a       
        [fa_new_i, lambda_new_i] = modeA_2d(TensorResidual,fa_old_i,fb_old_i,fc_old_i,lambda_old_i,conf.n,1);
        % mode b
        [fb_new_i, lambda_new_i] = modeA_2d(TensorResidual,fb_old_i,fa_new_i,fc_old_i,lambda_new_i,conf.n,1);     
        % mode c
        [fc_new_i, lambda_new_i] = modeA_2d(TensorResidual,fc_old_i,fa_new_i,fb_new_i,lambda_new_i,conf.n,1);              

        fa_old_i = fa_new_i;
        fb_old_i = fb_new_i;
        fc_old_i = fc_new_i;
        lambda_old_i = lambda_new_i;
    end
    TensorResidual = TensorResidual - thisTensor;
    fa_new(:,:,deflat_id) = fa_new_i;
    fb_new(:,:,deflat_id) = fb_new_i;
    fc_new(:,:,deflat_id) = fc_new_i;
    lambda_new((deflat_id-1)*conf.n*conf.n+1:deflat_id*conf.n*conf.n) = lambda_new_i;
end
estimate.f = fa_new;
estimate.lambda = lambda_new;
estimate.runningtime = toc;