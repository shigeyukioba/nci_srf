function [X,Z]=spikedetectNESTvoltage( x )
neuronIndex = unique( x(:,1) );
timeIndex = unique( x(:,2) );
NN = length( neuronIndex );
T = length(timeIndex);
X = zeros( NN, T );
for j = 1:NN
    i = neuronIndex(j);
    idx = find( x(:,1) == i );
    X( i, 1:length(idx) ) = x(idx,3)'-mean(x(idx,3));
end

Y = convn( X, [-1,3,-2], 'same' );
Y1 = convn( X, [-1,1,0], 'same' );
Y2 = convn( X, [0,1,-1], 'same' );
Z = (Y>30&Y1>0&Y2>0);
