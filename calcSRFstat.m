function o = calcSRFstat( R, varargin )

N = length( R ); % number of objective neurons
C = length( R(1).A ); % number of input neurons
M = length( R(1).SRF{1} );
K = length( R(1).A{1} );

%% setting options
doplot = false;
suridx = [];NS = 0;
for i=1:nargin-1
    if i<nargin && ischar(varargin{i})
      switch varargin{i}
          case 'plot'
              doplot = true;
          case 'surrogate'
              suridx = varargin{i+1};
              NS = length(suridx);
              N = N-NS;
      end
    end
end

%% initializing
z = zeros(N,C);
F = zeros(N,C,M);
A = zeros(N,C,K);
o.surface = z;
o.peak = z;
o.FS = z;
o.AS = z;
o.max = z;
o.min = z;
o.delay = z;

%% calc stats
for i=1:N
    for c=1:C
        f = R(i).SRF{c};
        F(i,c,:)=f;
        A(i,c,:)=R(i).A{c};
        o.surface(i,c) = sum( abs(f) );
        [o.peak(i,c), o.delay(i,c)] = max( abs(f) );
        o.max(i,c) = max( f );
        o.min(i,c) = min( f );
    end
end

%%
FS = reshape( F(:,suridx,:), N*NS, M );
iSigmaFS = inv( FS' * FS + eye(size(FS,2)) );
AS = reshape( A(:,suridx,:), N*NS, K );
iSigmaAS = inv( AS' * AS + eye(size(AS,2)) );
for i=1:N
    for c=1:C
        f = R(i).SRF{c};
        o.FS(i,c) = f*iSigmaFS*f';
        a = R(i).A{c};
        o.AS(i,c) = a*iSigmaAS*a';
    end
end

%% visualize
if doplot
    names = fieldnames( o );
    nn = length(names)
    figure
    for f=1:nn
       subplot( 2,4,f )
       x = o.(names{f});
       imagesc( x );
       title( names(f) )
       if length(suridx)>0
           hold on
           plot( N+[0,0], [0,N], 'r--')
       end
    end
    figure
    for f1=1:nn
        for f2=1:f1-1
            subplot(nn-1,nn-1,(f1-2)*(nn-1)+f2)
            x = o.(names{f2});
            y = o.(names{f1});
            x1 = x(:, 1:N ); y1 = y(:,1:N);
            x0 = x(:,suridx); y0 = y(:,suridx);
            plot( x1(:), y1(:), 'r.', ...
                x0(:), y0(:), 'b.')
            xlabel(names{f2})
            ylabel(names{f1})
        end
    end
end
