function N = shuffle( N, mode )
[C,T] = size( N );
if nargin==1
    mode = 1;
end

switch mode
    case 2
    for i=1:C
        N(i,:) = N(i,randperm(T));
    end
    case 1
    for i=1:T
        N(:,i) = N(randperm(C),i);
    end        
end

