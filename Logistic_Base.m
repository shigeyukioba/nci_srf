function X = Logistic_Base( x, param )
% X = Gamma_Base( x, param )
% Apply Gamma distribution base function
% Input x( 1, N )
% Output X( K, N )

if nargin == 1
    param.K = 50;
    param.step = 1;
end
K = param.K;
M = length(x);
X = zeros(K,M);
k = 1:K
m = (param.step*k).^2; % mean = a/b
v = (param.step*k)/3; % variance a/b^2
for k=1:K
  X(k,:) = 1 ./ (1 +exp( -(x-m(k))/v(k)) );
end
X = X ./ repmat( max( X,[],2 ), 1, M );
 