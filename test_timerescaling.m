% test timerescaling analysis

T = 100000;
t = 1:T;
pt = 0.01 * (1 + sin( 2*pi*t/200)) / 2;
Nt = rand(1,T) < pt;


figure
ooo = timerescaling( Nt, pt );
title( sprintf('KSstat=%f',ooo.KSstat))

