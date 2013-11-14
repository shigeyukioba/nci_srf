function X = Gamma_Base( x, param )
% X = Gamma_Base( x, param )
% Apply Gamma distribution base function
% Input x( 1, N )
% Output X( K, N )

if nargin == 1
    param.K = 50;
    param.step = 1;
    param.doNormalize = false;
end
K = param.K;
M = length(x);
X = zeros(K,M);
for k=1:K
  m = param.step*k^2; % mean = a/b
  v = param.step*k^2; % variance a/b^2
  b = m/v; a = m/b;
  tmp = (a-1)*log(x)-b*x;
  X(k,:) = tmp - max(tmp);
end
X = exp(X);
X = X ./ repmat( max( X,[],2 ), 1, M );
if param.doNormalize
  X = X ./ repmat( sum(X,2), 1,M );
end
 