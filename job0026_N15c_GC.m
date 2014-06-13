% job0020_GayaSRF


%load rec/simu0025_N15c.mat
load rec/simu0025_N15d.mat


if false
figure;
bar( sum(N0'))
xlabel('Neuron index','FontSize',14)
ylabel('#spikes in 6e5 frames','FontSize',14)
xlim([0,16])
end

%%
C0 = 10, iC0 = 1:C0;
T = 50000;
N = N0(:,1:T);
[C, T] = size(N);

N0 = N(iC0, [(301:T), (1:300)] );
size(N0)

caseno = 6


switch caseno
    case 1
      eta = 3;
      EK = 2;
      mw = 200
      Cidx = 1:15;
      fn = sprintf('rec/R0026d_%03d_C%d_C0%d_EK%d.mat',caseno,C,C0,EK);
    case 2
      eta = 3;
      EK = 0;
      mw = 200
      Cidx = 1:15;
      fn = sprintf('rec/R0026d_%03d_C%d_C0%d_EK%d.mat',caseno,C,C0,EK);
    case 3
      eta = 3;
      EK = 0;
      mw = 200
      C = 8; Cidx = [1:5, 6,7, 11];
      fn = sprintf('rec/R0026d_%03d_C%d_C0%d_EK%d.mat',caseno,C,C0,EK);
    case 4
      eta = 3;
      EK = 2;
      mw = 200
      C = 8; Cidx = [1:5, 6,7, 11];
      fn = sprintf('rec/R0026d_%03d_C%d_C0%d_EK%d.mat',caseno,C,C0,EK);
    case 5
      eta = 3;
      EK = 0;
      mw = 200
      C = 3; Cidx = [1, 6, 11];
      fn = sprintf('rec/R0026d_%03d_C%d_C0%d_EK%d.mat',caseno,C,C0,EK);
    case 6
      eta = 3;
      EK = 2;
      mw = 200
      C = 3; Cidx = [1, 6, 11];
      fn = sprintf('rec/R0026d_%03d_C%d_C0%d_EK%d.mat',caseno,C,C0,EK);
end
E = spikePCA(N(Cidx,:), ones(1,mw)/mw );
E = E(:,1:EK)';
[R_est,o,ot] = estimateSRF( [N(Cidx,:);N0],'iN',1:C, ...
      'MaxEpoch',100,'M',100,...
      'eta',eta,'Bases','Gamma','Link','Logistic',...
      'GCLR',true,'external',E);

fn
save( fn, 'R_est','o','ot' )

%%
figure
plotSRF(R_est,'out',1:10,'in',51:60,'plusminus','BG','ylim',[-1,1]*3)

