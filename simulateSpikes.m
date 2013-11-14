function N = simulateSpikes( R, T, E )
% N = simulateSpikes( R, T, E )
% Simulates and generates spike sequence
% R: Structure variable of SRFs
% T: time length
% E: External input (optional)
%  E.E(k,t) time course of  k  external input 
%  E.A(k,c) coefficient of the c-th neuron to the k-th input

NN = length( R ); % number of neurons
D = length( R(1).SRF{1} ); % length of spike response delay
N = sparse( NN, T );
if nargin == 3
    isexternal = true;
else
    isexternal = false;
end
for t=2:T
    idx = 1:(min( t-1, D ));
    Nh = N(:,t-idx);
    for i=1:NN
        if t==2
            dum = zeros(NN,1);
            for c=1:NN
                dum(c)=R(i).SRF{c}(1);
            end
            Rh{i} = dum;
        elseif t<=D+1
            for c=1:NN
                dum(c)=R(i).SRF{c}(t-1);
            end
            Rh{i} = [Rh{i}, dum];
        end
%        size(Nh)
%        size(Rh)
%        size(Rh{1})
        L(i) = R(i).A0 + Nh(:)'*Rh{i}(:);
        if isexternal
            L(i) = L(i) + E.E(:,t)'*E.A(:,i);
        end
    end
    L = exp( L );
    dt = 0.001;
    p = 1 - exp( -L*dt );
    N(:,t) = (rand(1,NN)<p);
    if mod(t,1000)==0
        disp( sprintf( '%d / %d', t, T ) )
    end
end

