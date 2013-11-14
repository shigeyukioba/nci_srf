function [LR,LR0,y] = LRrearrange0030( recLR, recRL )
% Sorts likelihood ratio
%  recLR : likelihood ratio
%  recRL : response label of 0,1,2,3 ... etc
flag_good = ~isnan(recLR)&(recLR>0)&(recLR<1e10);
idx1 = ( flag_good & recRL>1 ); nidx1 = sum(idx1(:));
idx0 = ( flag_good & recRL==0 ); nidx0 = sum(idx0(:));
LR0 = recLR(idx0); LR = recLR( idx1|idx0 );
y = recLR*0; y(idx1)=1; y = y( idx1|idx0 );
[dum,idx] = sort(-LR); LR = LR(idx); y = y(idx);
