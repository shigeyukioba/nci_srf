function DoF = chi2fit( x, doplot )
% Estimate degree of freedom of chi2distribution
candidate = [0.02:0.02:4 5:0.5:10 11:30];
for i = 1:length(candidate)
    loglikelihood(i) = sum( log( eps+chi2pdf( x, candidate(i) ) ) );
end
[dum, i] = max( loglikelihood );
DoF = candidate(i)
mean(x)
var(x)
DoF = mean(x)

if nargin == 1
    doplot = false;
end

if doplot
plot( candidate, loglikelihood )
%%
xlabel('DoF of chi2dist')
ylabel('Log Likelihood')
hold on
plot( candidate(i), dum, 'o' )
text( candidate(i), dum, sprintf('DoF=%g',candidate(i)))
end
