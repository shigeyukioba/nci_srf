function draw_distance( c, mx, my, stat )
%
stl = {'b.','r.'};
C = size(c,1);
for i=1:C
    for j=1:C
        if c(i,j)~=0;
            d = sqrt( (mx(i) - mx(j))^2 + (my(i)-my(j))^2 ) ; 
            plot( log10(d), log10( stat(i,j) ), stl{c(i,j)},...
                'MarkerSize',15)
            hold on
        end
    end
end
xlabel( 'Distance[\mu m]')
set(gca, 'XTick', [0,1,2], ...
    'XTickLabel', {'10^0','10^1','10^2'}, ...
    'YTick', [0,1,2], ...
    'YTickLabel', {'10^1','10^2','10^3'})

