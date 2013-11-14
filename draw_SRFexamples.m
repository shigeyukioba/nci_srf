function draw_SRFexamples( R_est, q )
C = length(q);
CI = repmat( (1:C)', 1,C );
CJ = repmat( (1:C), C, 1 );
[dum, idx] = sort( q(:) );

col = {'b-','g-','r-'};
tcol = {[0,0,1], [0,1,0], [1,0,0]};
spi = 1; count = 1; txt = '';
for i0 = 1:36
   i = CI(idx(i0));
   j = CJ(idx(i0));
   M = 50;
   srf = R_est(i).SRF{j}(1:M);
   srf2 = R_est(j).SRF{i}(1:M);
   subplot( 4,3,spi)
   if count == 1
      plot([0,0],[-2,9],'-',...
          'Color',[1,1,1]*0.5,'LineWidth',2) 
      hold on
      plot([-1,1]*M,[0,0],'-',...
          'Color',[1,1,1]*0.5,'LineWidth',2) 
   end
   srf0 = (srf(1)+srf2(1))/2;
   dsrf = [srf2(end:-1:1),srf0,srf];
   plot(-M:M,dsrf , col{count}, 'LineWidth',1.5 )
   xlim([-M,M])
   hold on
   %txt = [txt, sprintf('(%d,%d)',i,j)];
   text(-45, 11-count*2, sprintf('(%d,%d)',i,j), 'Color', tcol{count})
   count = count + 1;
   if count == 4
      txt = {txt, sprintf('q<%5.3g',q(idx(i0)))};
      title(txt); txt = '';
      spi = spi+1;
      count = 1;
   end
   ylim([-2,10])
   set(gca,'YTick',[0,5],'XTick',[-50,0,50],'XTickLabel',[-500,0,500])
end
set(gcf,'Position', [200,200,800,800])
