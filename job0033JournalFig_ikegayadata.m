%
if true
load ../gaya2011/rec072_2.mat
save ikegaya_data_u.mat u
end
load ikegaya_data
%%
figure
[mx, my] = roi_plot( u, 1:60 );
set(gca,'XTick',[],'YTick',[])
set(gcf,'Position',[100,100,1200,500])
saveAsPNGHQ('fig/fig0033ROI.png')
%%

%load rec/R0024_003_C60_C010EK2.mat
%load rec/R0024_002_C60_C010EK0.mat
%load rec/R0024_005_C60_C010EK0.mat
load rec/R0024_006_C60_C010EK0.mat
doLR=false;

C = 60; C0 = 10;

param.doplot = false;
param.polynomial = [1,2,3];
param.eta = [1, 3, 10, 30];
param.rowweight = 1000;
param.MaxEpoch = 10000;

%param.eta = [1, 3, 10, 30, 100, 300, 1000]*10;
%param.alg = 'poslog'; param.polynomial = 0;

if doLR
LR0 = ot.LR(1:C, C+(1:C0));
LR = ot.LR(1:C, 1:C);
[qLR,pi0_est] = ...
    EBayes_qvalue_rowwise( LR0, LR, param);
figure
subplot(2,2,1); plot( q(:), LR(:), '.')
subplot(2,2,2); imagesc(q)
end

%%
st = calcSRFstat( R_est, 'surrogate', C+(1:C0) );
names = {'AS','peak','surface'};
%names = fieldnames( st );
nn = length(names);
figure;
for i=1:nn
   subplot( 1,nn,i )
   x = st.(names{i});
   x1 = x( 1:C, 1:C );
   x0 = x( 1:C, C + (1:C0) );
   [q,pi0_est] = ...
    EBayes_qvalue_rowwise( x0, x1, param);
   qr{i} = q;
   plot( log10(qLR(:)), log10(q(:)), '.')
   xlabel('GC Log_10 q-value')
   ylabel( [names{i}, ' Log_10 q-value'])
end
%%
figure
for i=1:nn+1
   subplot(2,2,i)
   if i==1
       if doLR
      qq = qLR; nm = 'LR';
       else
           qq = zeros(C,C);
       end
   else
      qq = qr{i-1}; nm = names{i-1};
   end
   qq1 = qq*0;
   qq1(qq<0.5)=1;qq1(qq<0.3)=2;qq1(qq<0.2)=3;
   qq1(qq<0.1)=4;qq1(qq<0.01)=5;
   imagesc( qq1)
   title([nm, ' q-value'])
end
%%
figure(1);
qq = qr{1};
for i=1:nn
   subplot( 1,nn,i )
   x = st.(names{i});
   x1 = x( 1:C, 1:C );
   plot( qq(:), x1(:), '.')
   xlabel('Likelihood Ratio')
   %xlim([min(LR(:)),max(LR(:))])
   ylim([min(x1(:)),max(x1(:))])
   ylabel( names{i})
end

%%
%% qq = qLR;
qq = qr{1}; 
figure
qq = qq -diag(diag(qq)) + eye(length(qq));
draw_SRFexamples( R_est, qq )
saveAsPNGHQ('fig/fig0033SRFexamples_2.png')

%%
flag = qq<0.1;
c = zeros(C,C); c( flag ) = 2;
a = sign( st.max(1:C,1:C)+st.min(1:C,1:C));
c( flag & a==-1) = 1;
figure
[mx, my] = roi_plot( u, 1:60 );
set(gca,'XTick',[],'YTick',[])
set(gcf,'Position',[100,100,1200,500])
draw_connection( c, mx, my )
title(sprintf('%d significant links (q<0.1)',sum(c(:)~=0)))
saveAsPNGHQ('fig/fig0033connection.png')
%%


%%
figure;
draw_distance( c, mx, my, st.delay )
ylabel('Latency [msec]')
saveAsPNGHQ('fig/fig0033distance.png')
