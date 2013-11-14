function cc = crosscorrelogram( N1, N2, M )
cc = zeros( 1, M*2 + 1 );
T = length(N1);
idx = find( N1 );

for i=idx
    %i
    n1 = max( 1, i-M );
    n2 = min( T, i+M );
    cc( M+1+ ((n1-i):(n2-i)) ) ...
        = cc( M+1+ ((n1-i):(n2-i)) ) + N2( n1:n2 );
end
