function [qvalue, pi0_est, lfdr] = EBayes_qvalue_rowwise( x0, x, param )
% Calculates Empirical Bayesian qvalue
% Input: x0(i,j), x(i,j)  test statistics from null and all samples
%  i  and  j are row and column indices.
%  Common row  i  has common parameter in test statistic.

if nargin == 2
    param.doplot = false;
    param.polynomial = 2;
    param.rowweight = 1000;
    param.eta = 1;
end
NK = length(param.polynomial);

for K1 = 1:NK
[n0r,n0c] = size(x0); n0 = n0r*n0c;
[nr,nc] = size(x); n = nr*nc;
y = [ones(n0,1); zeros(n,1)];
xd = [x0(:);x(:)];
switch param.alg
        case 'poslog' 
            pm.K = 30; pm.step = 1;
            X = Logistic_Base( xd, pm );
            X = [ones(length(xd),1), X'];
        case 'polynomial'
          K = param.polynomial(K1);
          disp(sprintf('K=%d',K))
          if param.polynomial>1
            xq = xd/mean(xd);
            for k=2:param.polynomial
                xd = [xd, xd.^k]; 
            end
          end
end
% row-wise 
xdr = zeros( length(y), nr );
for i=1:nr
    xdr(i+nr*(0:(n0c+nc-1)),i)=param.rowweight;
end
%figure
%imagesc(xdr)
X = [xdr,xd];
W = zeros( 1, size(X,2) );

[b,o] = GLM_crossvalidation( y, X, W, param ); 
L(K1) = o.Lmax;
brec{K1}=b;
end
[dum,K1] = max(L);
K = param.polynomial(K1);
b = brec{K1};
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