function X = spikePCA( N, w )
% Calculate spike rate PCA

[C,T] = size(N);
N2 = conv2( full(N), w, 'same' );
[U,S,X] = svds(N2, C);

