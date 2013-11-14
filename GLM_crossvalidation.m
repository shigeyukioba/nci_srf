function [b,o] = GLM_crossvalidation( y, X, W, param )
% [b,o] = GLM_crossvalidation( y, X, W, param )
% NewtonRaphson method to estimate GLM with crossvalidation
%  y( 1:T, 1 ) in (0, 1)
%  X( 1:T, 1:(NN+K+1) )
%  W( 1, 1:(NN+K+1) ) parameters
%  P( y( t ) = 1 | W ) = p( X(t,:) * W )
% param should include
%   param.eta 
% optionally
%   param.MaxEpoch and others.

FoldNum = 5;
idx0 = find(y==0); r0 = randperm( length(idx0) );
idx1 = find(y==1); r1 = randperm( length(idx1) );
dum = y*0; dum( idx0 ) = r0; dum( idx1 ) = r1;
fidx = mod( dum, FoldNum ) + 1;
if ~isfield(param,'MaxEpoch')
   param.MaxEpoch = 100;
end
eta0 = param.eta;

for eta1 = 1:length(eta0)
    eta = eta0(eta1);
  for f = 1:FoldNum
    idx = ( fidx ~= f );
    [b,o] = GLM_NewtonRaphson( y(idx), X(idx,:), W, eta, param );
    idx = (fidx == f );
    Bx = X(idx,:)*b';
    p = 1./(1+exp(-Bx));
    loglikelihood(f) = sum( log( p( y(idx)==1 ) + eps ) ) ...
        + sum( log( 1-p( y(idx)==0 ) )+eps );
  end
  ml = mean(loglikelihood);
  sd = std(loglikelihood);
  disp( sprintf( 'eta=%g: L=%g (sd %g)', eta, ml, sd ) )
  recml(eta1) = ml;
end

[dum, besteta1 ] = max( recml );
eta = eta0( besteta1 );
b = GLM_NewtonRaphson( y, X, W, eta, param );
o.L = recml;
o.Lmax = dum;
o.eta0 = eta0;
o.eta = eta;
