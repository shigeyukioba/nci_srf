switch 3
    case 3
     recfilename = 'rec/R0026d_001_C15_C010_EK2.mat'; C=15; C0=10;
     load( 'rec/simu0025_N15d.mat' ) %'R_est','o','ot'
     T=50000; Cidx = 1:15;
    case 4
     recfilename = 'rec/R0026d_002_C15_C010_EK0.mat'; C=15; C0=10;
     load( 'rec/simu0025_N15d.mat' ) %'R_est','o','ot'
     T=50000; Cidx = 1:15;
end
Y = N0( Cidx, 1:T);
load(recfilename)

%%
LR = ot.LR(1:C, 1:C);
LR0 = ot.LR(1:C, C+(1:C0));

%%
param.doplot = false;
param.polynomial = 2;
param.rowweight = 1000;
param.eta = 1;
%[q,pi0_est] = EBayes_qvalue( LR0(:), LR(:), param);
[q,pi0_est] = EBayes_qvalue_rowwise( LR0, LR, param);
param

st = calcSRFstat( R_est, 'surrogate', C+(1:C0) );
%%
names = fieldnames( st );
nn = length(names);
figure;
for i=1:nn
   subplot( 1,nn,i )
   x = st.(names{i});
   x1 = x( 1:C, 1:C );
   plot( q(:), x1(:), '.')
   xlabel('GC q-value')
   ylabel( names{i})
end
figure(1);
for i=1:nn
   subplot( 1,nn,i )
   x = st.(names{i});
   x1 = x( 1:C, 1:C );
   plot( LR(:), x1(:), '.')
   xlabel('Likelihood Ratio')
   xlim([min(LR(:)),max(LR(:))])
   ylim([min(x1(:)),max(x1(:))])
   ylabel( names{i})
end
