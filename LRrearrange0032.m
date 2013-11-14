function [LR,LR0,y] = LRrearrange0032( recLR, recRL )
% Sorts likelihood ratio
%  recLR(1:numtrials, 1:(C*(C+C0))) : likelihood ratio
%  recRL(1:numtrials, 1:(C*C)       : response label of 0,1,2,3 ... etc
numtrials = size(recLR, 1 );
C=15; C0=10;
LR=[]; y=[]; LR0 = [];
RL = recRL(1,:);
yorg = RL*0; yorg( RL>0 )=1;
for t = 1:numtrials
    a = recLR(t,:);
    a = reshape( a, C, C+C0 );
    a1 = a( 1:C, 1:C );
    a0 = a( 1:C, C+(1:C0) );
    LR = [LR; a1(:)];
    y = [y; yorg'];
    LR0 = [LR0; a0(:)];
end
