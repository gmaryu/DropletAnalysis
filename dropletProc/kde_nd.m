function [density,xmesh]=kde_nd(data,bandwidth,n,MIN,MAX)
% Kernel density estimator for one-dimensional data, designed to be pretty
% fast and easy to understand.  Bandwidth must be set.
%
% Modified from Botev's kde. by Jason Levy
%
% INPUTS:
%   data    - N x dim matrix of data; Note: a row vector is treated as a 
%               single datapoint of whatever dimension.
%   bandwidth - the selected bandwidth, a 1xdim vector
%   n  - the number of mesh points used in the uniform discretization of the
%               interval [MIN, MAX]; also a 1xdim vector.  The entries of n
%               are rounded up to powers of two for maximum speed.
%   MIN, MAX  - defines the interval [MIN,MAX] on which the density estimate is constructed;
%               the default values of MIN and MAX are:
%               MIN=min(data)-Range/10 and MAX=max(data)+Range/10, where
%               Range=max(data)-min(data); data values outside this range
%               are ignored.
% OUTPUTS:
%     density - dim-dimensional array of length 'n' with the values of the density
%               estimate at the grid points;
%     xmesh   - the grid over which the density estimate is computed;
%
% Examples:
%
% [density, xmesh] = kde_nd([0; 1], 0.05); 
% figure, plot(xmesh,density)
%
% % Remember, a row vector would be interpreted as a single 2D entry
%
% [density, xmesh] = kde_nd([0 2; 1 3], 0.05); 
% figure;surf(xmesh{:},density'), xlabel('x')
%
% [density, xmesh]=kde_nd([1,2,3;2,5,4],[0.05 0.05 0.05],[2^7 2^9 2^8]);
% figure;isosurface(xmesh{:},permute(density,[2 1 3]),100), xlabel('x')

[~, dim] = size(data); 

if nargin<5 %define the default  interval [MIN,MAX]
    minimum=min(data,[],1); maximum=max(data,[],1); %=usual min/max unless data is a row vector
    Range=maximum-minimum;
    MIN=minimum-Range/10; MAX=maximum+Range/10;
elseif any(MIN>MAX)
    error('specified left endpoint is larger than specified right endpoint!')
end

if nargin<3 || numel(n)==0 % if n is not supplied switch to a default
    n=2^(ceil(16/dim))*ones(1,dim);
end
n=2.^ceil(log2(n)); % round up n to the next power of 2;
n(MIN==MAX)=1; %if only one data point in this direction, don't need smoothing.


%ensure that all entries are between MIN and MAX by deleting exceptions
%Only matters if MIN, MAX are manually set
data(any(data < repmat(MIN,[size(data,1) 1]),2),:)=[];
data(any(data > repmat(MAX,[size(data,1) 1]),2),:)=[];


% set up the grid over which the density estimate is computed; 
xmesh=cell(1,dim); imesh=xmesh;
for i=1:dim
    xmesh{i} = linspace(MIN(i),MAX(i),n(i))'; 
    if dim==1, xmesh=xmesh{1}; end  % a vector instead of a cell
    imesh{i} = pi * (0:n(i)-1) * bandwidth(i)/(MAX(i)-MIN(i));
    if n(i)==1, imesh{i}=0; end
end
[I{1:dim}]=ndgrid(imesh{:});

% bin the cropped data uniformly using the limits defined above
initial_data=ndhist(data,MIN,MAX,n);

a=dct_nd(initial_data); % discrete cosine transform of initial data

% build DCT of the Gaussian kernel; up to a scalar factor, h_hat is dct_nd
% of product in each direction of (skipping all the p's)
% exp(-1/2 ((1:n)-0.5)/n /sigma) with sigma=bandwidth/(MAX-MIN)
% The scalar factor ensures that density is properly normalized at the end,
% since a_t(1)=a(1).
%
% Verification:
% x=-1/2+(1:n)';
% s=n*(bandwidth/(MAX-MIN));
% h=exp(-1/2*(x/s).^2); h= h/sum(h);
% sum(abs(h_hat - dct_nd(h))) % is a tiny number, nonzero because of
%                             % discretization and rounding (not sure which
%                             % dominates.
%
P=cat(dim+1,I{:});
h_hat=prod(exp(-P.^2 / 2), dim+1);

a_t = a.* h_hat; 

% now apply the inverse discrete cosine transform
density=idct_nd(a_t) *prod(n)/prod(MAX(n>1)-MIN(n>1));
% the n>1 excludes those dimensions with MIN=MAX and no scaling is done.
end

function binned_data=ndhist(data,MIN,MAX,n)
% this function computes the histogram of an dim-dimensional data set;
% 'data' is N rows by dim columns, all lies between MIN and MAX
% xmesh is a 1xdim cell, listing the bin endpoints in each direction
% n is a row vector listing the number of bins used in each dimension

[N,dim]=size(data);

if dim==1, n=[n 1]; end %MATLAB is dumb with 1d data


% histc uses the endpoints of binranges, as opposed to centres, so we'll
% offset to have the centres of our binranges be the points in "xmesh" above.
dx = (MAX-MIN)./(n-1);
bins=zeros(size(data));
for i=1:dim
    if n(i)==1, bins(:,i)=1; continue, end
    [~,bins(:,i)] = histc(data(:,i), ...
        linspace(MIN(i) - dx(i)/2, MAX(i) + dx(i)/2, n(i)+1));
end

% Combine the vectors of 1D bin counts into a grid of dimD bin
% counts.
%-----------------edited by Shixuan Liu-----------------------
bins(bins==0) = 1;
%-----------------------------------------------------------------------
binned_data = accumarray(bins,1/N,n);
end


function Y=dct_nd(X)
%Compute the discrete cosine transform of arbitrary-dimensional data.  Note
%that the normalization is a little funny in the first index; the inverse
%normalization naturally takes place in the inverse transform idct_nd.

dim=length(size(X)); %will say 2 if X is a vector, otherwise correct.

Y = X;
for p = 1:dim
%hard to adjust things just in the pth dimension, so we'll change it to the
%first and then change back.

    Z=permute(Y,[p setdiff(1:dim, p)]);
    np = size(Z,1);  
 
    if np==1, continue, end %for example if MATLAB's confused about a row vector
    
    Z(:,:)=Z([1:2:end, end:-2:2],:); %re-order; note that end=np. 
    
    sz=size(Z); sz(1)=[];  %size of Z except for 1st dim.
    weight = repmat(exp(-1i* (0:np-1)' *pi/(2*np)), [1 sz]);
    
    Z = real(weight.*fft(Z));  %along 1st dimension
    Y=ipermute(Z,[p setdiff(1:dim, p)]);
end

end





function Y=idct_nd(X)

dim=length(size(X)); %will output 2 if X is a row vector, otherwise correct.

Y = X;
for p = 1:dim
%hard to adjust things just in the pth dimension, so we'll change it to the
%first and then change back.

    np = size(Y,p);  
    if np==1, continue, end %for example if MATLAB's confused about a row vector
    
    Z=permute(Y,[p setdiff(1:dim, p)]);
    
    sz=size(Z); sz(1)=[];  %size of Z except for 1st dim.
    weight = repmat([1;2*exp(1i* (1:np-1)' *pi/(2*np))], [1 sz]);
 
    Z = real(ifft(weight.*Z)); %along 1st dimension

    Z([1:2:np, np:-2:2],:) = Z(:,:); %permute back
    
    Y=ipermute(Z,[p setdiff(1:dim, p)]);
end

end
