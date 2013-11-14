function FDPplot( q, y, stl )

[dum,idx] = sort(q); 
fdp = cumsum( 1-y(idx) ) ./ [1:length(y)]';
sq = q(idx);

plot( [0,1],[0,1],'k-')
hold on
plot( sq, fdp, stl, 'LineWidth',2)
axis([0,1,0,1])
xlabel('qvalue')
ylabel('fdp')
