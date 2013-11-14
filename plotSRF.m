function plotSRF( R, varargin )
% function plotSRF( R, varargin )
% Plots response function(s)
% Usage:
%  plotSRF( R )  plots all 
%  plotSRF( R, 'out', [1,2,3], 'in', [1:10] )
%     plots SRFs of 3*10 pairs of neurons
%  plotSRF( R, 'ylim', [-2,2])
%     sets ylim([-2,2]) for all panel


%% Default setting
if isstruct(R)
    idx_i = 1:length(R);
    idx_c = 1:size(R(1).SRF,2);
end
yl = [-1,1]*0.5;xl='auto';
stl = 'b-';
plusminus = false;
stackmode = false;
ispatch = false;
drawBG = false;
pvaluemode = false;
width=1;
%% Specific settings
for i=1:nargin-1
    if i<nargin & isstr(varargin{i})
      switch varargin{i}
          case 'out'
              idx_i = varargin{i+1};
          case 'in'
              idx_c = varargin{i+1};
          case 'ylim'
              yl = varargin{i+1};
          case 'xlim'
              xl = varargin{i+1};
          case 'patch'
              ispatch = true;
              patchcolor = varargin{i+1};
          case 'style'
              stl = varargin{i+1};
          case 'plusminus'
              plusminus = true;
          case 'width'
              width = varargin{i+1};
          case 'BG'
              drawBG = true;
          case 'stack'
              stackmode = true;
          case 'pvalue'
              pvalue = varargin{i+1};
              pvaluemode = true;
          case 'ClusterNo'
              cn = varargin{i+1};
              ccol{1} = [0,0,0];
              ccol{2} = [0,0,1];
              ccol{3} = [0,1,0];
              ccol{4} = [1,0,0];
              ccol{5} = [0.7,0.7,0];
              for i=5:(max(cn(:))+1)
                  dum=sqrt([3,5,7])*i;
                  ccol{i} = dum - floor(dum);
              end
      end
    end
end

%%
Ni = length(idx_i);
Nc = length(idx_c);

if isstr(xl)
  switch xl
  case 'auto'
    if plusminus
        xl2 = [-1,1]*length(R(1).SRF{1});
    else
        xl2 = [0,1]*length(R(1).SRF{1});
    end
  end
else
  xl2 = xl;
end


for i0=1:Ni
  for c0=1:Nc
      if ~stackmode
        subplot(Ni,Nc, (i0-1)*Nc + c0 )
        hold on
      end
      if drawBG | (stackmode & i0==1 & c0==1)
           plot( [0,0],yl,'k-', xl2,[0,0],'k-' )
      end
      if exist('cn')
          cc = ccol{cn(i0,c0)};
      else
          cc = [];
      end
      r = R(idx_i(i0)).SRF{idx_c(c0)};
      plotSRFsingle( r, stl, cc, width );
      if pvaluemode
        [dy,dx] = max(abs(r));
        h = text( dx, min(dy,yl(2))*sign(r(dx)),...
            pvaluetext( pvalue(idx_i(i0),idx_c(c0)) ));
        set(h,'FontSize',12)
      end

      if plusminus
        r = R(idx_c(c0)).SRF{idx_i(i0)};
        plotSRFsingleMinus( r, stl, cc, width );
      end
      if ~stackmode
      if i0==1
            h=title(sprintf('No:%d',idx_c(c0)));
            set(h,'FontSize',12)
      end
      if i0~=Ni
            set(gca,'XTick',[])
      end
      if c0==1
            h=ylabel(sprintf('No:%d',idx_i(i0)));
            set(h,'FontSize',12)
      else
            set(gca,'YTick',[])
      end
      end
      ylim(yl)
      xlim(xl)
  end
end

end

%%
function s = pvaluetext( p )
s = '';
if p<0.05
    s = '*';
end
if p<0.01
    s = '**';
end
if p<0.001
    s = '***';
end
end
%%

function plotSRFsingle( r, stl, cc, w )
 n = length(r);
 x = [1:n, n, 1];
 y = [r, 0, 0];
 if length(cc)==3
    plot( r, stl, 'Color', cc, 'LineWidth', w )
 else
    plot( r, stl, 'LineWidth', w )
 end
 hold on
end

%%
function plotSRFsinglepatch( r, cc )
 n = length(r);
 x = [1:n, n, 1];
 y = [r, 0, 0];
 patch( x, y, 'Color', cc )
end

%%
function plotSRFsingleMinus( r, stl, cc, w )
if length(cc)==3
  plot( -length(r):1:-1,r(end:-1:1), stl,'Color',cc,'LineWidth', w )
else
  plot( -length(r):1:-1,r(end:-1:1), stl,'LineWidth', w )
end
hold on
end
