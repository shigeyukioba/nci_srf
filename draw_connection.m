function draw_connection( c, mx, my )
%
plot( mx, my, 'bo' )
hold on
stl = {'b-','r-'};
C = size(c,1);
for i=1:C
    for j=1:C
        if c(i,j)~=0;
            arrow( [mx(i);my(i)], [mx(j);my(j)], stl{c(i,j)} )
        end
    end
end


function arrow( x0, x1, stl )

d = x1 - x0;
d = d / sqrt( sum( d.^2 ) );
theta = 15/360 * 2*pi;
d1 = 3*[cos(theta) sin(theta);-sin(theta) cos(theta)]*d;
theta = -theta;
d2 = 3*[cos(theta) sin(theta);-sin(theta) cos(theta)]*d;

x = [x0 x1 x1-d1 x1 x1-d2];
plot( x(1,:), x(2,:), stl, 'LineWidth', 2 )
