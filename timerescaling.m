function [o, th] = timerescaling( N, p, doplot )
% function  o = timerescaling( N, p )
% calculates statistics on time rescaling analysis
% Input
%  N( 1:T, 1 ) in {0,1} spike train data
%  p( 1:T, 1 ) in [0,1] spike probability
%

% calculate cumulative distribution function
[C,T] = size(N);
if ~(C==size(p,1) & T==size(p,2))
    error('Input size mismatch')
end
if nargin<3
doplot = false;
end

for c=1:C
    idx = find( N(c,:) );
    ns = length(idx); % #spikes
    i = 1;
    cdf = zeros( ns, 1 );
    for j=1:length(idx)
        id0 = i:idx(j)-1;
        cdf( j ) = 1 - exp( - sum( p( c, id0 ) ) );
        i = idx(j);
    end
    scdf = sort(cdf);
    F = ((1:ns)-1/2 )/ns;
    Fu = min( F + 1.36 / sqrt(ns), 1 );
    Fl = max( F - 1.36 / sqrt(ns), 0 );
    Fn = scdf';
    d = F - Fn;
    o(c) = max( abs( d ) );
    th(c) = 1.36 / sqrt(ns);
    if doplot
        subplot(ceil(sqrt(C)),ceil(sqrt(C)),c)
        plot( F, F, 'k-', ...
            F, Fu, 'k--', ...
            F, Fl, 'k--')
        hold on
        plot(F, Fn, 'r-','LineWidth',2)
        title(sprintf('%d:D=%4.3g',c,o(c)))
    end
end
