function estimate = ALS_2d_noDeflation(conf, Tensor)
tic
OriginalTenosrNorm = sqrt(sum(sum(abs(Tensor).^2)));
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
    tmp_lambda = zeros(conf.n,conf.n,conf.L);
    tmp_lambda(1:conf.a,1:conf.a,:)=1;
    lambda_old = conf.lambda * tmp_lambda(:);
    % lambda_old = conf.lambda*ones(conf.n*conf.n*conf.L,1);
    lambda_new = lambda_old;
end

INTERVAL = 10;
%% Iterative Steps
for iter = 1 : conf.maxIter
    % mode a
    [fa_new, lambda_new] = modeA_2d(Tensor,fa_old,fb_old,fc_old,lambda_old,conf.n,conf.L);
    % mode b
    [fb_new, lambda_new] = modeA_2d(Tensor,fb_old,fa_new,fc_old,lambda_new,conf.n,conf.L);
    % mode c
    [fc_new, lambda_new] = modeA_2d(Tensor,fc_old,fa_new,fb_new,lambda_new,conf.n,conf.L);
    
    fa_old = fa_new;
    fb_old = fb_new;
    fc_old = fc_new;
    lambda_old = lambda_new;
    % orthogonalize every few iteration
    if mod(iter,INTERVAL)==0
        fa_new = my_orthogonalization_2d(fa_new);
        fb_new = my_orthogonalization_2d(fb_new);
        fc_new = my_orthogonalization_2d(fc_new);
    end
end

% compute tensor residual
thisTensor = zeros(size(Tensor));
for id_filter = 1  : conf.L
    tmp_fa = fa_new(:,:,id_filter);
    tmp_fb = fb_new(:,:,id_filter);
    tmp_fc = fc_new(:,:,id_filter);
    tmp_lambda = lambda_new((id_filter-1)*conf.n*conf.n+1:(id_filter)*conf.n*conf.n);
    thisTensor = thisTensor + circulant_2d(tmp_fa)*diag(tmp_lambda)*(matrix_katri_rao(circulant_2d(tmp_fc),circulant_2d(tmp_fb)))';
    
end
TensorResidual = Tensor - thisTensor;

estimate.TensorResidual = TensorResidual;
estimate.TensorResidualNorm = sqrt(sum(sum(abs(TensorResidual).^2)));
estimate.OriginalTenosrNorm = OriginalTenosrNorm;
estimate.f = fa_new;
estimate.lambda = lambda_new;
estimate.runningtime = toc;