for iii=3
switch iii
    case 1
setting1 = 1; setting2 = 1; 
figfilename='fig/fig0032A_ROCK_N15d.png'
    case 2  
setting1 = 2; setting2 = 1;
figfilename='fig/fig0032B_ROCK_HH004.png'
    case 3
setting1 = 1; setting2 = 2;  % Figure A
figfilename='fig/fig0032C_ROCS_N15d.png'
    case 4
setting1 = 2; setting2 = 2;  % Figure A
figfilename='fig/fig0032D_ROCS_HH004.png'
end


switch setting1
    case 1
        dataname = 'N15d'; version0='b';
    case 2
        dataname = 'HH004';version0 = 'b';
end
C=15; C0 = 10;

switch setting2
    case 1
        algcase0 = {'Kim0','Kim1','Kim2','Smooth'}
        rocidx = [1,4]
        stl      = {'b:','b--','b-','r-'};
    case 2
        algcase0 = {'Smooth','peak','AS','surface'}
        stl      = {'r-','g-','m-','c-'};
        rocidx = [1,2,3,4];
end
%NT0 = [1000,2000,5000,10000,20000,50000]
NT0 = [1000,5000,10000,50000]
NNT0 = length(NT0);
addpath ~/NAISTtoolbox/lightspeed/

param.doplot = false;
param.polynomial = [3:4];
param.eta = [1, 3, 10, 30, 100, 300, 1000]/1000;
param.alg = 'poslog'; param.polynomial = 0;

%%
figure
for NT1 = 1:length(NT0)
  for alg1 = 1:length(algcase0)
    NT = NT0(NT1) 
    algcase = algcase0{alg1}

    switch algcase
       case 'Kim0' %Original K
          recfilename = sprintf('rec/rec0031_%s%sNT%d%sEK0.mat',...
              dataname,'Kim',NT,version0)
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'rec_R_est', 'o'
          [LR,LR0,y] = LRrearrange0032( recLR, recRL );
          [K,M] = size(rec_R_est{1}(1).Bases);NU = K;
          pvalue = 1 - chi2cdf( LR, NU );
          v = ranking(pvalue);
          pi0 = 1;
          q = min( pi0 * length(pvalue) * pvalue ./ v, 1.0 );
          %[q, pi0_est(NT1,alg1)] = qvalue( pvalue );
       case 'Kim1' %chi2 fit
          recfilename = sprintf('rec/rec0031_%s%sNT%d%sEK0.mat',...
              dataname,'Kim',NT,version0)
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'rec_R_est', 'o'
          [LR,LR0,y] = LRrearrange0032( recLR, recRL );
          [K,M] = size(rec_R_est{1}(1).Bases);
          NU = chi2fit( LR0 );
          pvalue = 1 - chi2cdf( LR, NU );
          v = ranking(pvalue);
          pi0 = 1;
          q = min( pi0 * length(pvalue) * pvalue ./ v, 1.0 );
          %[q, pi0_est(NT1,alg1)] = qvalue( pvalue );
       case 'Kim2' % Empirical Bayes
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'R_est', 'o'
          [LR,LR0,y] = LRrearrange0032( recLR, recRL );
         [q,pi0_est(NT1,alg1)] = EBayes_qvalue( LR0, LR, ...
             param);
       case 'Smooth'
          recfilename = sprintf('rec/rec0031_%s%sNT%d%sEK0.mat',...
              dataname,'Smooth',NT,version0)
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'R_est', 'o'
          [LR,LR0,y] = LRrearrange0032( recLR, recRL );
         [q,pi0_est(NT1,alg1)] = EBayes_qvalue( LR0, LR, ...
             param);
       case 'SmoothEK2'
          recfilename = sprintf('rec/rec0031_%s%sNT%d%sEK2.mat',...
              dataname,'Smooth',NT,version0)
          load(recfilename) %'recL1', 'recL0', 'recLR', 'recRL', 'R_est', 'o'
          [LR,LR0,y] = LRrearrange0032( recLR, recRL );
         [q,pi0_est(NT1,alg1)] = EBayes_qvalue( LR0, LR, ...
             param);
        case 'peak'
          so = [];
          for trial = 1:length(rec_R_est)
            stat = calcSRFstat( rec_R_est{trial} );
            s = stat.peak(1:C,1:(C+C0));
            so = [so; s(:)']; 
          end
          [s, s0, y] = LRrearrange0032( so, recRL );
          [q,pi0_est(NT1,alg1)] = EBayes_qvalue( s0, s, param);
          LR = s;
        case 'AS'
          so = [];
          for trial = 1:length(rec_R_est)
            stat = calcSRFstat( rec_R_est{trial} );
            s = stat.AS(1:C,1:(C+C0));
            so = [so; s(:)']; 
          end
          [s, s0, y] = LRrearrange0032( so, recRL );
          [q,pi0_est(NT1,alg1)] = EBayes_qvalue( s0, s, param);
          LR = s;
        case 'surface'
          so = [];
          for trial = 1:length(rec_R_est)
            stat = calcSRFstat( rec_R_est{trial} );
            s = stat.surface(1:C,1:(C+C0));
            so = [so; s(:)']; 
          end
          [s, s0, y] = LRrearrange0032( so, recRL );
          [q,pi0_est(NT1,alg1)] = EBayes_qvalue( s0, s, param);
          LR = s;
    end
    % Draw fdp
    subplot(2,NNT0,NT1+NNT0)
    FDPplot( q, y, stl{alg1} )
    % Draw ROC
    subplot(2,NNT0,NT1)
    or(NT1, alg1) = drawROCcurve( LR, y, stl{alg1} );
    auc( alg1 ) = or( NT1, alg1 ).ROCsurface;
    axis([0,1.01,0,1.01])
  end
  ss = '%4.2g';
  for i00 = 2:length(rocidx)
      ss = [ss, ',%4.2g']
  end
  %title( {sprintf('T=%d',NT), sprintf(['AUC=(' ss ')'], auc(rocidx))})
  title( {sprintf('T=%d,AUC=',NT), sprintf(['(' ss ')'], auc(rocidx))})
end
set( gcf, 'Position', [130,451,900,470])
drawnow
%%
saveAsPNGHQ( figfilename )
%set(gcf,'PaperPositionMode','auto')
%print('-dpng', '-r150',figfilename)
end %for iii
