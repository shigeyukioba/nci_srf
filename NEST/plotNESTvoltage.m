function Z = plotNESTvoltage( x )

neuronIndex = unique( x(:,1) )
w = 100;
for j = 1:length(neuronIndex)
    i = neuronIndex(j);
    idx = find( x(:,1) == i );
    plot( x(idx,2), x(idx,3) + j*w, '-' )
    hold on
end

[dum, Z] = spikedetectNESTvoltage( x );
for j=1:length(neuronIndex)
    i = neuronIndex(j);
    idx = find( Z(i,:)==1 );
    plot( idx*2, idx*0 + j*w, 'ro' )
end
