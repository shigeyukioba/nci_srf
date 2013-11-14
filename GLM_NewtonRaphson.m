function [betaR,o] = GLM_NewtonRaphson(y, X, W, eta, param)
% NewtonRaphson method to estimate GLM
%  y( 1:T, 1 ) in (0, 1)
%  X( 1:T, 1:(NN+K+1) )  Input, X(:,1) must be a constant term.
%  W( 1, 1:(NN+K+1) ) parameters.
%  P( y( t ) = 1 | W ) = p( X(t,:) * W )
T = length(y);
NP = length(W);
idx1 = ( y==1 );
idx0 = ( y==0 );
past_E = inf;

paramDefault = struct(...
    'MaxEpoch',20,...
    'linkfunction', 'Logistic', ...
    'divergenceThreshold', 1e-6, ...
    'convergenceThreshold', 1e-6, ...
    'issilent', true,...
    'dofixedhessian', true, ...
    'ispositive', false);
if nargin <5
    param = paramDefault;
else
    names = fieldnames( paramDefault );
    for i=1:length(names)
        if ~isfield( param, names{i} )
            param.(names{i}) = paramDefault.(names{i});
        end
    end
end

W0 = W(1);W(1) = 0;
oldW = W;
counter = 0;
for epoch = 1:param.MaxEpoch
    %% Calc. p.
    Bx = X*W'+W0;
    switch param.linkfunction
        case 'Logistic'
            p = 1./(1+exp(-Bx));
        case 'Poisson';
            p = 1-exp(-exp(Bx));
    end
    if false
    figure(1)
    subplot(2,1,1)
    hold off
    %plot( sort( log10(abs(W)) ) );
    plot( sort( W ) );
    hold on
    %plot( 1, log10(abs(W(1))), 'ro')
    plot( 1, W(1), 'ro')
    %ylabel('log10(abs(W))')
    ylabel('W')
    subplot(2,1,2)
    hold off
    plot( log10( sort( p(idx1))), 'r-')
    hold on
    plot( log10( sort( p(idx0))), 'b-')
    ylabel('log10(p))')
    drawnow
    end
    
    %% Evaluate Reg-term
    Regterm = sum(W(2:end).^2)/2; % Omit W(1)
    %Regterm = sum(W(1:end).^2)/2;
    %% Evaluate Likelihood
    %LogLikelihood = sum(log(p(idx1)+eps))+sum(log(1-p(idx0)+eps));
    LogLikelihood = sum(log(p(idx1)))+sum(log(1-p(idx0)));
    Li = LogLikelihood; % The larger is the better fit
    E = -Li+eta*Regterm; % The smaller is the better
    o.LogLikelihood = Li;
    o.eta = eta;
    o.E = E;
    if past_E<E
        isupward=true;
        txtupward = sprintf('UPWARD%d',counter);
    else
        isupward=false;
        txtupward = '';
    end
    if ~param.issilent
    fprintf( 'Epoch=%2d, Log-likelihood=%f, R=%f, E=%f %s\n',...
              epoch, Li, Regterm, E, txtupward)
    end
    if epoch>2 && (E-past_E)/past_E > param.divergenceThreshold
        disp(sprintf('=== Diverged at epoch %d ===',epoch));
        o.E = nan;
        break
    end
    if abs( past_E-E )/E < param.convergenceThreshold
        disp(sprintf('Converged at epoch=%d',epoch))
        break
    else
        past_E = E;
    end
        counter = 0;
        oldW = W;
    %% Update SRF
    switch param.linkfunction
        case 'Logistic'
          a = p .* (1-p);
          z = y - p;
        case 'Poisson'
          a = y .* (1-p)./p.^2;
          z = (y - p)./p;
    end
    
    g = z'*X; % Gradient
    if param.dofixedhessian
        if epoch == 1
            H = - 0.25 * X' * X - eta*diag([0,ones(1,NP-1)]); % Fixed Hessian
            %H = - X' * ( X .* repmat(a,1,size(X,2)) )...
            %    - eta*diag([0,ones(1,NP-1)]); % Fixed Hessian
            %invH = inv(H);
            r = chol(-H);
        else
        end
        g = g - eta * [0,W(2:end)];
        %u = - g / H;
        %u = - g * invH;
        u = solve_triu(r,solve_tril(r',g'));% Uses Lightspeed toolbox by Tom Minka
        u = u';
        switch 1
            case 1
         ug = u * g'; % scalor
         ux = X * u'; % T dimensional column vector
         uHu = a' * ux.^2 + eta * sum( u(2:end).^2 ); % scalor
         W = W + (ug/uHu) * u;
         if param.ispositive
             idx = find( W(2:end)<0 );
             W(idx+1)=0;
         end
            case 2
         W = W + u;
        end
    else
        % Newton Raphson update
        H = X' * ( X .* repmat(a,1,size(X,2)) ); % Hessian
        W = (W*H+g)/( H + eta*diag([0,ones(1,NP-1)]) ); % Omit W(1)
        %W = (W*H+g)/( H + eta*diag(ones(1,NP)) );
    end
end
W(1) =W(1)+W0;
betaR = W;
o.p = p;
