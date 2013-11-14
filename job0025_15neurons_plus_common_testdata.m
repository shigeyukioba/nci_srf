% job0025_15neurons_plus_common_testdata.m
% modifying job0014_15neurons_testdata.m


Na = 5; % num of neurons in a group
RLa = [1,0,3;
       2,1,0;
       0,2,1];
%A0a = [4,3,3]; %case a
%A0a = [2,1,1];  %case b
%A0a = [1.5,1.2,1.2]; %case c
A0a = [3.5,2.5,2.5]; %case d
RL = [];A0 = [];
for i=1:size(RLa,1)
    RLl =[];
    for j=1:size(RLa,2)
        switch RLa(i,j)
            case 0
                RLb = zeros(Na,Na);
            case 1
                RLb = eye(Na);
            case 2
                RLb = zeros(Na,Na);
                for k=1:Na
                    RLb(k,randperm(Na,2))=2;
                end
            case 3
                RLb = zeros(Na,Na);
                for k=1:Na
                    RLb(k,randperm(Na,2))=3;
                end                
        end
        RLl =[RLl,RLb];
    end
    RL = [RL;RLl];
    A0 = [A0, A0a(i)*ones(1,Na)];
end
NN = length(RL) % num of neurons;

s = 1:100;
B1 = exp( -0.2*s );
B2 = exp( -0.4*s );
B3 = exp( -0.6*s );
Rt{0+1} = B1*0;
Rt{1+1} = -B2;
Rt{2+1} = (B1-B3)*2;
Rt{3+1} = -(B2-B3)*5;

clear R
for i=1:NN
    Rc.SRF = cell(1,NN);
    Rc.A0=A0(i);
    for c=1:NN
        Rc.SRF{c} = Rt{ RL(i,c)+1 };
    end
    R(i) = Rc;
end


T = 600000;
E.E = sin( (1:T)*2*pi/1000 )+1;  % Common Input
E.A = rand( 1, NN ); % Coefficient to the common input
datafilename = sprintf('rec/simu0025_N%dd.mat',NN)
if true
  N0 = simulateSpikes( R, T, E );
  save( datafilename, 'N0', 'R', 'RL' );
end
numspikes = sum(N0')

%%
figure
imagesc(RL)
xlabel('Pre neuron index')
ylabel('Post neuron index')

figure
plotSRF( R, 'xlim',[0,30], 'ylim',[-1,1] )
%%
