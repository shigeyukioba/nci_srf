% job0020_GayaSRF

load ikegaya_data.mat

figure;
bar( sum(Y'))
xlabel('Neuron index','FontSize',14)
ylabel('#spikes in 1min','FontSize',14)
xlim([0,61])

%%
C = 60
C0 = 10, iC0 = (1:C0)*2;
N = sparse(Y(1:C,:));
[C, T] = size(N);
N0 = N(iC0, [(301:T), (1:300)] );
size(N0)

caseno = 6

M=100;MaxEpoch=100;
switch caseno
    case 1
      eta = 0.3;
      EK = 2;
      mw = 200
fn = sprintf('rec/R0024_%03d_C%d_C0%d.mat',caseno,C,C0)
    case 2
      eta = 0.3;
      EK = 0;
      mw = 200
fn = sprintf('rec/R0024_%03d_C%d_C0%dEK0.mat',caseno,C,C0)
    case 3
      eta = 1.0;
      EK = 2;
      mw = 200
fn = sprintf('rec/R0024_%03d_C%d_C0%dEK2.mat',caseno,C,C0)
    case 4
      eta = 1.0;
      EK = 2;
      mw = 200
fn = sprintf('rec/R0024_%03d_C%d_C0%dEK0.mat',caseno,C,C0)
    case 5
      eta = 1.0;
      EK = 0;
      mw = 200
fn = sprintf('rec/R0024_%03d_C%d_C0%dEK0.mat',caseno,C,C0)
    case 6
      eta = 1.0;
      EK = 0;
      mw = 200
      fn = sprintf('rec/R0024_%03d_C%d_C0%dEK0.mat',caseno,C,C0)
      %%
      K = 25;
      M = 300;
      MaxEpoch = 200;
      step = 0.4;
      Bases = zeros(K,M);
      for k=1:K
         t = 1:M;
         m = step*k^2; % mean = a/b
         v = step*k^2; % variance a/b^2
         b = m/v; a = m/b;
         tmp = (a-1)*log(t)-b*t;
         Bases(k,:) = tmp - max(tmp);
      end
      Bases = exp(Bases);
      Bases = Bases ./ repmat( max( Bases,[],2 ), 1, M );
      plot( Bases' )          
      %%
end

E = spikePCA(N, ones(1,mw)/mw );
E = E(:,1:EK)';

%[R_est,o,ot ] = estimateSRF( [N;N0],'iN',1:C, ...
%      'MaxEpoch',MaxEpoch,'M',M,...
%      'eta',eta,'Bases',Bases,'Link','Logistic',...
%      'GCLR',true,'external',E);

[R_est,o] = estimateSRF( [N;N0],'iN',1:C, ...
      'MaxEpoch',MaxEpoch,'M',M,...
      'eta',eta,'Bases',Bases,'Link','Logistic',...
      'GCLR',false,'external',E);

fn
save( fn, 'R_est','o','ot' )

%%
figure
plotSRF(R_est,'out',1:10,'in',51:60,'plusminus','BG','ylim',[-1,1]*3)

