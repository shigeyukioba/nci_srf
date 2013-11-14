function qvalue = chi2fit_bayesian( x0, x )
% Bayesian estimation of degree of freedom of chi2distribution
% with calculation of empirical Bayesian p-values
% Input: x0, x
%  statistics from null and all samples
n0 = length(x0); n = length(x);
range = [0.1, 50];
nbins = 100;
logbincenters = log(range(1)) + ...
    (log(range(2))-log(range(1)))*[0:nbins-1]/nbins;
bin = exp( logbincenters );
binprior = conv( bin, [0.5, 0, -0.5]);
binprior = binprior(1:end-2); binprior = binprior/sum(binprior);

prchi2bin = zeros(nbins, nbins);
prchi2_0 = zeros(nbins, n0 );
prchi2_1 = zeros(nbins, n );
cdfchi2 = zeros(nbins, n );
for i = 1:nbins
    prchi2_0(i,:) = chi2pdf( x0, bin(i) );
    prchi2_1(i,:) = chi2pdf( x, bin(i) );
    prchi2bin(i,:) = chi2pdf( bin, bin(i) );
    cdfchi2(i,:) = 1 - chi2cdf( x, bin(i) );
end
loglikelihood = sum( log( eps+prchi2_0 ), 2 )';
likelihood = exp( loglikelihood );
posterior = likelihood.*binprior;
posterior = posterior / sum( posterior );
[dum, i] = max( likelihood ); % Map index

pvalueBayes = posterior * cdfchi2;
pvalueMAP = cdfchi2(i,:);
w=0.025; pvbin = w/2:w:1;
hBayes = hist( pvalueBayes, pvbin );
hMAP = hist( pvalueMAP, pvbin );

y = [ones(n0,1); zeros(n,1)];
%X = [ones(length(y),1),[x0;x]];
%X = [ones(length(y),1),[x0;x], [prchi2_0, prchi2_1]'];
%X = [ones(length(y),1), [prchi2_0, prchi2_1]'];
xd = [x0;x];
X = [xd*0+1, xd, xd.^2/10, xd.^3/100,xd.^4/1000];
W = zeros( 1, size(X,2) );

eta = 1e-6;
b = GLM_NewtonRaphson( y, X, W, eta ); 
Bx = X*b';
r = 1./(1+exp(-Bx));
p0 = r( 1:n0 );
p = r / mean(p0); % local fdr
p0 = p0/mean(p0);
p1 = p( n0+(1:n) );

lfdr = p( n0+(1:n) );
[sldfr, idx] = sort(lfdr);
sq = cumsum( sldfr ) ./ [1:n]';
qvalue = x*0;
qvalue(idx) = sq;

return
plot( x, qvalue, '.')
xlabel('x'), ylabel('qvalue')
return
%%
plot( x, p1, 'r.', ...
    x0, p0,'b.')
return
%%
subplot(2,2,1)
likelihood = posterior;
plot( bin, likelihood )
xlabel('DoF of chi2dist')
ylabel('Likelihood')
hold on
plot( bin(i), dum, 'o' )
text( bin(i), dum, sprintf('DoF=%g',bin(i)))
DoF = bin(i);

subplot(2,2,3)
%plot( bin, prchi2bin, 'b-' )
hold on
plot( bin, posterior * prchi2bin, 'r-', 'LineWidth', 2 )

subplot(2,2,2)
plot( pvbin, hBayes, 'r-', ...
    pvbin, hMAP, 'b-' )
xlabel('P-value')
