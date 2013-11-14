function [R,o,ot] = estimateSRF( N, varargin )
% [R] = estimateSRF( N, iN, varargin )
% Estimates spike response function (SRF)
% Input
%  N: Sparse matrix of spikes
%     N(i,t)=1 denotes spike of the i-th neuron at time  t.
%
%  varargin: Options can be cast by 'OptionName', value
%   'iN', iN: specifies a set of neuron indices whose SRFs are estimated.
%      When  iN = [1,3], R(1) and R(3) are calculated and
%      R(2), R(4), ..., R(NN) are set at blank values.
%   'Bases', B: specifies a factor loading matrix  B(1:K,1:M).
%     If  B='Gamma', Gamma distribution density function is used as bases.
%     If  B='NoBase', B(1:M, 1:M) = eye(M) is used as bases.
%   'Link', 'Logistic': specifies a link function,
%      which should be either 'Poisson' or 'Logistic'. If omitted 'Poisson'
%      is used as a default.
%   'Rinit', R: specifies initial values of the SRFs.
%   'eta', eta: specifies a (set of) regularization hyperparameter(s)
%      When eta is a vector,
%      corresponding SRFs  R{1:length(eta)} are output as a cell array.
%   'Validation', Nv: specifies independent validation data.
%      If it is specified, test-likelihood is output as  o(i).E
%   'OmitSRF': if specified, R(:).SRF is not calculated.
%       This setting reduces memory consumption.
%   'external', E: if specified, external input E is added.
% Output
%  R: Structure variable of SRFs
%    R(i).SRF{c}(s) denotes spike response
%      of the i-th neuron
%      to the spike at time  s * dt
%      of the c-th neuron.
%    R(i).A0 denotes constant term.
%    R(i).SRF{c} = 0  represents
%    R(i).SRF{c}(s) = 0 for all  s.
%     If multiple values of  eta  is specified,
%     R: Cell variable of multiple SRFs
%    R{ieta}  is the ieta-th structure whose value is as explained above.
%  o: Structure variable of some statistics
%    o.LogLikelihood(i)  is likelihood of  R(i)  calculated based on  N
%    when an option 'dovalidate' is specified with 'Rinit',
%    it is used as a cross validation score.
%  ot: Structure variable of causality test statistics
%

[NN, T] = size( N );
%% Set default values
iN = 1:NN; 
M=10;K=10;
eta=1;

param.MaxEpoch = 10;
param.convergenceThreshold = 1e-6;
param.divergenceThreshold = 0.1;
param.withbases = false;
param.doplot = false;
param.linkfunction = 'Poisson';
param.dovalidate = false;
param.dofixedhessian = true;
%param.dofixedhessian = false;
param.omitsrf = false;
param.GCLR = false; % do calculate pvalue?
param.external = [];
param.numExtInput = 0;
Bases = [];

%% Specific settings
for i=1:nargin-1
    if i<nargin && ischar(varargin{i})
      switch varargin{i}
          case 'iN'
              iN = varargin{i+1};
          case 'M'
              M = varargin{i+1};
          case 'K'
              K = varargin{i+1};
          case 'eta'
              eta = varargin{i+1};
          case 'Rinit'
              R = varargin{i+1};
          case 'MaxEpoch'
              param.MaxEpoch = varargin{i+1};
          case 'Bases'
              Bases = varargin{i+1};
              if length(Bases)>1
                 param.withbases = true;
              else
                 param.withbases = false;
              end
          case 'Link'
              param.linkfunction = varargin{i+1};
          case 'doplot'
              param.doplot = true;
          case 'dopvalue'
              param.GCLR = true;
          case 'Validation'
              param.dovalidate = true  ;
              Nv = varargin{i+1};
          case 'OmitSRF'
              param.omitsrf = true;
          case 'GCLR'
              param.GCLR = varargin{i+1};
          case 'external'
              param.external = varargin{i+1};
              param.numExtInput = size( param.external, 1 );
          case 'param'
              param = varargin{i+1};
      end
    end
end

%% If length(eta) is not unity, start multiple hyperparameter mode
if length(eta)>1
    numEtaCandidate = length(eta);
    if exist('R')
       Rinit = R;
    end
    R = cell(numEtaCandidate,1);
    o.LogLikelihood = zeros( numEtaCandidate, length(iN) );
    o.TestLogLikelihood = zeros( numEtaCandidate, length(iN) );
    for ieta = 1:numEtaCandidate
        fprintf('Calling estimateSRF with eta=%f (%d/%d)\n',...
            eta(ieta),ieta,numEtaCandidate);
        [Rout,o1,ot1] = estimateSRF(N,'eta',eta(ieta),...
                'K',K,'M',M,'iN',iN,'Bases',Bases,...
                'param', param, ...
                'Validation',Nv);
        Rinit = Rout;
        R{ieta} = Rout;
        for i=1:length(iN)
          o.LogLikelihood(ieta,i) = o1(i).LogLikelihood;
          o.TestLogLikelihood(ieta,i) = o1(i).E;
          o.eta(ieta) = eta(ieta);
        end
        ot(ieta)=ot1;
    end
    return
end

%% Initialize SRF
if ~exist('R')
for i=1:NN
    Rc.SRF = cell(1,NN);
    p = sum(N(i,:))/T;
    switch param.linkfunction
        case 'Logistic'
      Rc.A0 = log( p/(1-p) );
        case 'Poisson'
      Rc.A0 = log( -log(1-p) );
    end
    if ~param.omitsrf
    for c=1:NN
        Rc.SRF{c} = zeros(1,M);
    end
    end
    Rc.AE = zeros(1,param.numExtInput);
    if param.withbases
        if ischar(Bases)
        switch Bases
            case 'Gamma'% Gamma whose peak is normalized to be one
                Bases = zeros(K,M);
                for k=1:K
                    t = 1:M;
                    m = 0.5*k^2; % mean = a/b
                    v = 0.5*k^2; % variance a/b^2
                    b = m/v; a = m/b;
                    Bases(k,:) = t.^(a-1) .* exp( - b*t );
                end
                Bases = Bases ./ repmat( max( Bases,[],2 ), 1, M );
            case 'NGamma' % Gamma whose area is normalized to be one
                Bases = zeros(K,M);
                for k=1:K
                    t = 1:M;
                    m = 0.5*k^2; % mean = a/b
                    v = 0.5*k^2; % variance a/b^2
                    b = m/v; a = m/b;
                    Bases(k,:) = t.^(a-1) .* exp( - b*t )...
                                        * b^a/gamma(a);
                end
            case 'NoBase'
                K = M;
                Bases = eye(K);             
        end
        end
        Rc.Bases = Bases;
        [K,M] = size(Bases);
        for c=1:NN
            Rc.A{c} = zeros(1,K);
        end
    end
    R(i)=Rc;
end
else
    if isfield(R(1), 'Bases')
        Bases = R(1).Bases;
        param.withbases = true;
    end
end
if ~param.withbases
    K = M;
end

fprintf('Starting %s regression\n',param.linkfunction)

X = applyBases( N, Bases );
if param.dovalidate
  Xv = applyBases( Nv, Bases );
end

X = [X, param.external'];

%% Parameter estimation for each output neuron  i
for i0=1:length(iN)
    i = iN(i0);
    fprintf( 'Estimating SRF of Neuron:%02d (%d/%d)\n',...
        i, i0, length(iN) )
    W = R2W( R, i, param, 1:NN );
    %
    [W,o(i0)] = GLM_NewtonRaphson(N(i,:)', X, W, eta, param);
    %
    if param.dovalidate
        param2 = param;
        param2.MaxEpoch = 1;
        [dum,ov] = GLM_NewtonRaphson(Nv(i,:)', Xv, W, eta, param2);
        o(i0).E = ov.LogLikelihood;
    end
    %
    R = W2R( R, i, W, param, 1:NN );
    %
end % for each Neuron

if param.GCLR
    ot = GCLR( N, R, X, eta, param, iN )
end

function ot = GCLR( N, R, X, eta, param, iN )
NN = length(R)
if param.withbases
  K = length(R(1).A{1});
else
  K = length(R(1).SRF{1});
end

%% Parameter estimation for each output neuron  i
for i0=1:length(iN)
    i = iN(i0);
    W1 = R2W( R, i, param, 1:NN );
    [dum,o] = GLM_NewtonRaphson(N(i,:)', X, W1, eta, param);
    L1 = o.LogLikelihood;
    ot.L1(i,1) = L1;
    for c = 1:NN
        fprintf('Evaluating GC at (%d,%d)\n',i,c)
        nidx = [1:(c-1),(c+1):NN];
        xidx = 1;
        for cc = nidx
            xidx = [xidx, (cc-1)*K+[1:K]];
        end
        xidx = [xidx, cc*K+[1:param.numExtInput]];
        W0 = R2W( R, i, param, nidx );       
        [dum,o] = GLM_NewtonRaphson(N(i,:)', X(:,xidx), W0, eta, param);
        L0 = o.LogLikelihood;  
        ot.L0(i,c) = L0;
    end
end % for each Neuron
LR = -2 * ( ot.L0 - repmat(ot.L1,1,NN) );
DOF = K;
ot.pvalue = 1 - chi2cdf( LR, K );
ot.LR = LR;


function W = R2W( R, i, param, nidx )
NN = length(nidx);
if param.withbases
  K = length(R(i).A{1});
else
  K = length(R(i).SRF{1});
end
NE = param.numExtInput;
W = zeros(1, NN*K+1+NE);
for c=1:NN
  if param.withbases
      W((c-1)*K+(1:K)+1)=R(i).A{nidx(c)}(1:K);
  else
      W((c-1)*K+(1:K)+1)=R(i).SRF{nidx(c)}(1:K);
  end
end
W(NN*K+1+(1:NE)) = R(i).AE;
W(1) = R(i).A0;

function R = W2R( R, i, W, param, nidx )
NN = length(nidx);
NE = param.numExtInput;
if param.withbases
  K = length(R(i).A{1});
else
  K = length(R(i).SRF{1});
end
R(i).A0 = W(1);
for c=1:NN
   for k=1:K
      if param.withbases
         R(i).A{nidx(c)}(k) = W(1+(c-1)*K+k);
         if ~param.omitsrf
           R(i).SRF{nidx(c)} = R(i).A{nidx(c)}*R(i).Bases;
         end
      else
         R(i).SRF{nidx(c)}(k) = W(1+(c-1)*K+k);
      end
   end
end
R(i).AE = W(NN*K+1+(1:NE));



function X = applyBases( N, Bases )
[NN,T] = size(N);
[K,M] = size(Bases);
X = zeros(T, NN+K+1);
X(:,1) = ones(T,1);
for c=1:NN
  for k=1:K
    dum = conv( double(full(N(c,:))), Bases(k,1:end) );
    X(:,(c-1)*K+k+1) = [0,dum(1:T-1)];
  end
end


function [betaR, dbetaR] = calcbeta( N, z, A0, Rh )
% A function that is no longer used; left as a reference
[NN,T] = size(N);
[NN,K] = size(Rh);

    betaR = zeros(1,1+K*NN);
    dbetaR = betaR;
    betaR(1) = A0;
    dbetaR(1) = sum( z );
    for c=1:NN
        idx = find(N(c,:));
        for k=1:K
           idxs = idx((idx+k)<=T)+k;
           betaR(1+(c-1)*K+k) = Rh(c,k);
           dbetaR(1+(c-1)*K+k) = sum( z( idxs ) );
        end
    end
    
function [betaR, dbetaR] = calcbetaB( N, z, A0, Rh )
% A function that is no longer used; left as a reference
[NN,T] = size(N);
[NN,K] = size(Rh);

betaR = zeros(1,1+K*NN);
dbetaR = betaR;
betaR(1) = A0;
dbetaR(1) = sum( z );

global BN

for c=1:NN
    for k=1:K
       betaR(1+(c-1)*K+k) = Rh(c,k);
       dbetaR( 1+(c-1)*K+k ) = BN{c,k}*z;
    end
end


function Bx = calc_Bx( N, A0, Rh, Bases )
% A function that is no longer used; left as a reference
if nargin == 4 && length(Bases)>1
  [K,M] = size(Bases);
  withbase = true;
else
  [NN,M]=size(Rh);
  withbase = false;
  K=M;
end
[NN, T] = size(N);
if size( Rh, 1 )==NN && size( Rh, 2 )==K
    %OK
else
    disp( 'Something wrong' )
end
Bx = zeros(T,1);
for t=2:T
   idx = 1:(min( t-1, M ));
   Nh = N(:,t-idx);
   if withbase
       Rh0 = Rh * Bases(:,1:length(idx));
   else
       Rh0 = Rh(:,1:length(idx));
   end
   Bx(t) = A0 + Nh(:)'*Rh0(:);
end % t


function H = calcH_B( N, Bases, w )
% A function that is no longer used; left as a reference
global BN
[NN,T] = size(N);
[K,M] = size(Bases);
H = zeros( K*NN, K*NN );

for c=1:NN
  for cc=c:NN
  for k=1:K
  for kk=1:K
        H((c-1)*K+k,(cc-1)*K+kk) = BN{c,k}.*BN{cc,kk}*w;
  end
  end
  end
end
H00 = sum( w );
H0C = zeros( 1, K*NN );

for c=1:NN
   for k=1:K
     H0C(1,(c-1)*K+k) = BN{c,k}*w;
   end
end
H = [H00,H0C;...
     H0C',H];


function H = calcH_fast( N, M, w )
% A function that is no longer used; left as a reference
[NN,T] = size(N);
H = zeros( M*NN, M*NN );
for c=1:NN
   fprintf('-- calcH_fast. %d/%d',c,NN)
   for cc=1:NN
       idx_c = find( N(c,:) );
       for q=idx_c
           idxS0 = max(1,q-M+1);
           idxS1 = min(T,q+M-1);
           idx_c2 = find( N(cc,idxS0:idxS1) );
           idx_c2 = idx_c2+idxS0-1;
           for q2=idx_c2
               smin = max([1,q2-q+1]);
               smax = min([M,T-q,q2-q+M]);
               for s=smin:smax
                  ss = q + s - q2;
               H((c-1)*M+s,(cc-1)*M+ss)=...
                   H((c-1)*M+s,(cc-1)*M+ss) + ...
                   w(q+s);
               end
           end
       end
   end
end
H00 = sum( w );
H0C = zeros( 1, M*NN );
for c=1:NN
    for s=1:M
       idx_t=1:T;
       idx_t= find( idx_t-s>0 );
       if ~isempty(idx_t)
          idx_t = idx_t(  N(c,idx_t-s)  );
          H0C((c-1)*M+s)=sum(w(idx_t));
       end
    end
end
H = [H00,H0C;...
     H0C',H];

 
 
function H = calcH( N, M, w )
% A function that is no longer used; left as a reference
[NN,T] = size(N);
H = zeros( M*NN, M*NN );
for c=1:NN
   disp(sprintf('%d/%d',c,NN))
   for cc=c:NN
       for s=1:M
           for ss=1:M
               idx_t=1:T;
               idx_t= find( idx_t-max(s,ss)>0 );
               if length(idx_t)>0
               idx_t = idx_t( find( N(c,idx_t-s) ) );
               if length(idx_t)>0
               idx_t = idx_t( find( N(cc,idx_t-ss) ) );
               if length(idx_t)>0
               H((c-1)*M+s,(cc-1)*M+ss)=...
                 sum(w(idx_t));
               H((cc-1)*M+ss,(c-1)*M+s)=...
                 sum(w(idx_t));
               end
               end
               end
           end
       end
   end
end
H00 = sum( w );
H0C = zeros( 1, M*NN );
for c=1:NN
    for s=1:M
       idx_t=1:T;
       idx_t= find( idx_t-s>0 );
       if length(idx_t)>0
          idx_t = idx_t( find( N(c,idx_t-s) ) );
          H0C((c-1)*M+s)=sum(w(idx_t));
       end
    end
end
H = [H00,H0C;...
     H0C',H];
