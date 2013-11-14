function [qvalue, pi0_est, lfdr] = EBayes_qvalue( x0, x, param )
% Calculates Empirical Bayesian qvalue
% Input: x0(b,k), x(i,k)  test statistics from null and all samples
%  b  and  i  are indices of (simulated) null and observed null samples.
%  k  is an index of covariates

if nargin == 2
    param.doplot = false;
    param.polynomial = 2;
    param.eta = 1;
    param.ispositive = true;
    param.alg = 'poslog';
end
NK = length(param.polynomial);

for K1 = 1:NK
    n0 = size(x0,1); n = size(x,1);
    y = [ones(n0,1); zeros(n,1)];
    xd = [x0;x];
    switch param.alg
        case 'poslog' 
            pm.K = 30; pm.step = 1;
            X = Logistic_Base( xd, pm );
            X = [ones(length(xd),1), X'];
        case 'polynomial'
      K = param.polynomial(K1);
      disp(sprintf('K=%d',K))
      if K>1
        xd = xd/mean(xd);
        for k=2:K
          xd = [xd, xd.^k]; 
        end
      end
      X = [ones(n0+n,1), xd];
    end
    W = zeros( 1, size(X,2) );
   [b,o] = GLM_crossvalidation(  y, X, W, param ); 
   %b = GLM_NewtonRaphson( y, X, W, param.eta ); 
   L(K1) = o.Lmax;
   brec{K1} =b;
end
[dum, K1] = max( L );
K = param.polynomial(K1)
b = brec{K1};
L
Bx = X(:,1:length(b))*b';
r = 1./(1+exp(-Bx));
% Note: r( x ) approximates
%  n0*p0(x) / ( n0*p0(x) + n*p(x) )
%    = n0*p0(x) / ( (n*pi0+n0)*p0(x)+n*(1-pi0)*p1(x) ) 
r0 = r( 1:n0 );
r1 = r( n0 + (1:n) );

lambda = 0.1;
sr0 = sort( r0 );
r0bar = mean(r0(ceil(n0*lambda):end));
pi0_est = n0 * (1-r0bar) / (n*r0bar);
lfdr = n * pi0_est / n0 * ( r1 ./ (1-r1) );

pi0_est = 1.0;

[sldfr, idx] = sort(lfdr);
sq = cumsum( sldfr ) ./ [1:n]';
qvalue = x*0;
qvalue(idx) = sq;

if param.doplot
    plot( x, qvalue, '.')
    xlabel('x'), ylabel('qvalue')
end