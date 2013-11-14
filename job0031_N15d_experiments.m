switch 1.2
    case 1
      load( 'rec/simu0025_N15d.mat' ) %'N0','R','RL'
      dataname = 'N15d'
      NT0 = [1000,2000,5000,10000,20000,50000]
      maxepoch0 = [10,10,10,10,10,10]
      etaK = 10;  etaS = 1; version0 = 'b';
      C = 15; % Num neurons
      C0 = 10; % Num surrogates
      mw = 200;
      EK = 0; EKT='EK0';
    case 1.2
      load( 'rec/simu0025_N15d.mat' ) %'N0','R','RL'
      dataname = 'N15d'
      NT0 = [1000,2000,5000,10000,20000,50000]
      maxepoch0 = [10,10,10,10,10,10]
      etaK = 10;  etaS = 1; version0 = 'b';
      C = 15; % Num neurons
      C0 = 10; % Num surrogates
      mw = 200;
      EK = 2; EKT='EK2';
    case 2
      dataname = 'HH004'
%      etaK = 1; etaS = 2; version0 = 'a';
      etaK = 10; etaS = 1; version0 = 'b';
      load( 'rec/spikesHH004.mat' ) % T = 150000
      RL = eye(15);
      RL([1:5]+5, [1:5]) = ...
          [2,0,0,0,0;...
          2,2,2,2,2;...
          0,0,0,0,0;...
          0,0,2,2,0;...
          0,2,0,0,0];
      RL([1:5],[1:5]+10) = ...
          [0,0,3,3,0;...
          3,3,0,3,3;...
          3,0,3,0,0;...
          0,3,0,0,0;
          0,0,0,0,0];
      RL([1:5]+10,[1:5]+5) = ...
          [2,2,0,0,0;...
          2,0,0,2,2;...
          0,0,2,0,0;...
          0,2,2,0,2;
          0,0,0,2,0];
      NT0 = [1000,2000,5000,10000,20000,50000]
      maxepoch0 = [10,10,10,10,7,2]
            C = 15; % Num neurons
      C0 = 10; % Num surrogates
      mw = 200;
      EK = 0;EKT='EK0';
end      
% Experimental setting

% common
algcase0 = {'Kim','Smooth'}
for algcase1 = 1:2
    algcase = algcase0{algcase1}
%for q=1:length(NT0)
for q=5:6
NT = NT0(q) 
maxepoch = maxepoch0(q)
recfilename = sprintf('rec/rec0031_%s%sNT%d%s%s.mat',dataname,algcase,NT,version0,EKT)

if true
  recL1 = [];recL0 = [];recLR = [];recRL = [];
  clear rec_R_est
end

%% Main body of the experiment
tic
for epoch = 1:maxepoch
 switch algcase
     case 'Kim'
       eta = etaK;
       M = 20;
       base = 'NoBase';
     case 'Smooth'
       eta = etaS;
       M = 20;
       base = 'NGamma';
 end
 t = (1:NT) + (NT*(epoch-1));
 t_100 = [100:NT,1:99] + (NT*(epoch-1));
 N = [N0( :, t ); N0( 1:C0, t_100) ];
 
 E = spikePCA(double(N(1:C,:)), ones(1,mw)/mw );
 E = E(:,1:EK)';

 [R_est,o,ot] = estimateSRF( N,...
    'iN', 1:C, ...
    'MaxEpoch',100,...
    'eta',eta, 'M', M, ...
    'Bases',base,'Link','Logistic','GCLR',true,'external',E); 
 recL1 = [recL1; ot.L1(:)'];
 recL0 = [recL0; ot.L0(:)'];
 recLR = [recLR; ot.LR(:)'];
 recRL = [recRL; RL(:)'];
 rec_R_est{epoch} = R_est;
 save( recfilename, 'recL1', 'recL0', 'recLR', 'recRL', 'rec_R_est', 'o'  )
end
end
end
toc
