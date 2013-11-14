%dataname = ''; % <= N9
dataname = 'N15'; 
%dataname = 'HH004'
version0 = 'b';
%version0 = 'a';

%%
algcase0 = {'Kim1','Kim2','Smooth','peak','AS'}
%algcase0 = {'Smooth','AS'}
stl      = {'b--','b-','r-', 'g-','m-'};
NT0 = [1000,2000,5000,10000,20000,50000,100000,200000]
%NT0 = [2000,5000,10000,20000]
addpath ~/NAISTtoolbox/lightspeed/

%%
figure
for NT1 = 1:length(NT0)
  for alg1 = 1:length(algcase0)
    NT = NT0(NT1) 
    algcase = algcase0{alg1}

    param.doplot = false;
    param.polynomial = [2:4];
    param.eta = [1, 3, 10, 30, 100, 300, 1000];

    switch algcase
       case 'Kim1' %chi2 fit
          recfilename = sprintf('rec/rec0015_%s%sNT%d%s.mat',...
              dataname,'Kim',NT,version0)
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'R_est', 'o'
          % recL1( trial, link ): trial = 1:10, link = 1:225
          %    where 225 = 15 neurons *15 neuron
          [LR,LR0,y] = LRrearrange0030( recLR, recRL);
          [K,M] = size(R_est(1).Bases);NU = K;
          NU = chi2fit( LR0 );
          pvalue = 1 - chi2cdf( LR, NU );
          v = ranking(pvalue);
          pi0 = 1;
          q = min( pi0 * length(pvalue) * pvalue ./ v, 1.0 );
          %[q, pi0_est(NT1,alg1)] = qvalue( pvalue );
       case 'Kim2' % Empirical Bayes
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'R_est', 'o'
          [LR,LR0,y] = LRrearrange0030( recLR, recRL );
         [q,pi0_est(NT1,alg1)] = EBayes_qvalue( LR0((1:100)), LR, ...
             param);
       case 'Smooth'
          recfilename = sprintf('rec/rec0015_%s%sNT%d%s.mat',...
              dataname,'Smooth',NT,version0)
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'R_est', 'o'
          [LR,LR0,y] = LRrearrange0030( recLR, recRL );
         [q,pi0_est(NT1,alg1)] = EBayes_qvalue( LR0((1:100)), LR, ...
             param);
        case 'peak'
          stat = calcSRFstat( R_est );
          [s, s0, y] = LRrearrange0030( stat.peak(:), recRL(1,:)' );
         [q,pi0_est(NT1,alg1)] = EBayes_qvalue( s0((1:100)), s, ...
             param);
        case 'AS'
          stat = calcSRFstat( R_est );
          [s, s0, y] = LRrearrange0030( stat.AS(:), recRL(1,:)' );
         [q,pi0_est(NT1,alg1)] = EBayes_qvalue( s0((1:100)), s, ...
             param);
    end
    % Draw fdp
    subplot(2,6,NT1+6)
    FDPplot( q, y, stl{alg1} )
    % Draw ROC
    subplot(2,6,NT1)
    or(NT1, alg1) = drawROCcurve( LR, y, stl{alg1} );
    auc( alg1 ) = or( NT1, alg1 ).ROCsurface;
    axis([0,1.01,0,1.01])
  end
  title( {sprintf('T=%d',NT), sprintf('AUC=(%5.3g, %5.3g)', auc)})
end
set( gcf, 'Position', [130,451,1512,470])
