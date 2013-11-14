function [lfdr, o] = nci_calclfdr( R, R0 )
% [lfdr, o] = nci_calclfdr( R, R0, varargin )
% Calculates local false discovery rate
% Input
%  R, R0 are SRFs of original and shuffled

NN = length(R); 
NN2 = length(R(1).SRF);
if NN2>NN
    appendmode = true;
    disp( 'Append mode is detected' );
else
    appendmode = false;
end
K = length( R(1).SRF{1} );
X = zeros( NN*(NN-1), K );
X0 = zeros( NN*(NN-1), K );
count = 0;
for i=1:NN
    for c=1:NN
        if i~=c
            count = count+1;
            if appendmode
               X(count,:) = R(i).SRF{c};
               X0(count,:) = R(i).SRF{c+NN};          
            else
            X(count,:) = R(i).SRF{c};
            X0(count,:) = R0(i).SRF{c};
            end
            Xi(count) = i;
            Xc(count) = c;
        end
    end
end

%% Strength statistics
strength = sum( X.^2, 2 );
strength0 = sum( X0.^2, 2);

%% Null subspace
[U,S,V] = svd( X0 );
ds = diag(S);
K2 = sum( cumsum( ds.^2 )/sum( ds.^2 ) < 0.99 )
W0 = V*diag( 1./ds ); % Transformation matrix
W0 = W0(:,1:K2);
W0back = diag(ds(1:K2))*V(:,1:K2)'; % Pseudo inverse transformation matrix
%% Application of the subspace
XW0 = X*W0;
XW00 = X0*W0;
strengthw0  = sum( XW0.^2, 2 );
strengthw00 = sum( XW00.^2, 2 );

%%
figure
for i=1:2
    switch i
        case 1
            s = log(strength); s0 = log(strength0);
            nm = 'Strength';
        case 2
            s = log(strengthw0); s0 = log(strengthw00);
            nm = 'NullStrength';
    end
    subplot(2,1,i)
    [h,bin] = hist([s,s0],30);
    h = hist(s,bin);
    h0 = hist(s0,bin);
    plot( bin,h,'r-',bin,h0,'b-')
    xlabel(nm)
end
%%
if false
figure
x0 = X(1,:);
x00 = x0 * W0 * W0back;
plot( x0, 'b-' )
hold on
plot( x00, 'r--' )
end
o.U = U;
o.V = V;
lfdr = [];

o.strength = reconstruct( strength, Xi, Xc );
o.strength0 = reconstruct( strength0, Xi, Xc );
o.strengthw0 = reconstruct( strengthw0, Xi, Xc );
o.strengthw00 = reconstruct( strengthw00, Xi, Xc );

function xr = reconstruct(x, i, c)
for j=1:length(x)
    xr( i(j), c(j) ) = x(j);
end
