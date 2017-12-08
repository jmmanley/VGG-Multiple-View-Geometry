% [P,X] = vgg_signsPX_from_x(P,X,x)  Finds signs of P and X in a projective reconstruction.
%
% Given a projective reconstruction, i.e. P, X, and x such that 
%   s_n^k x_n^k = P^k X_n,
% where 
%   - P^k is k-th camera matrix
%   - X_n is n-th scene point
%   - x_n^k is image projection of X_n in camera P^k
%   - s_n^k is scale,
% it will change signs of P^k and X_n such that all s_n^k are positive. Positivity of s_n^k
% determines the signs of P and X uniquely up to a single overall sign.
%
% Parameters:
%   P ... cell(K) of 3-by-4 matrices, camera matrices.
%     P also can be (3*K)-by-4 joint camera matrix.
%     P also can be 3-by-4-by-K array.
%   X ... double(4,N), scene points in homog. coordinates
%   x ... cell(K) of double(3,N), image points in homog. coordinates.
%     If an image point is missing, set x{k}(:,n) = [NaN;NaN;NaN].
%     x also can be joint (3*K)-by-N joint image point matrix, again with NaNs if a point is missing.
%     x also can be 3-by-N-by-K array.
%
% If it is not possible to change signs of P^k and X_n such that s_n^k are positive, it is P=X=[].
% This means that the projective reconstruction [P,X,x] does not correspond to any real scene.
%
% The function works for any dimension, ie, D-by-(D+1) camera matrices.
%
% See also vgg_selfcalib_qaffine.

function [P0,X] = vgg_signsPX_from_x(P0,X,u0)

% Re-arrange input data to joint camera matrix and joint image points.
if iscell(P0)
  P = vertcat(P0{:});
  u = vertcat(u0{:});
else
  if ndims(P0)==2 % joint camera matrix
    P = P0;
    u = u0;
  else % P(:,:,k) is k-th camera matrix
    P = [];
    u = [];
    for k = 1:size(P0,3)
      P = [P; P0(:,:,k)];
      u = [u; u0(:,:,k)];
    end
  end
end

[D N] = size(X);
K = size(P,1)/(D-1);

% Do sign swapping in joint image / joint camera matrix format
if any(isnan(u(:)))
  [P,X] = signsPX_from_x_occl(P,X,u); % slower code but can handle undefined points
else
  [P,X] = signsPX_from_x(P,X,u); % faster code if all image points are defined
end

if isempty(P)
  P0 = [];
  return
end

% Re-arrange back to original format
if iscell(P0)
  for k = 1:length(P0)
    P0{k} = P([1:D-1]+(D-1)*(k-1),:);
  end
else
  if ndims(P0)==2
    P0 = P;
  else
    for k = 1:size(P0,3)
      P0(:,:,k) = P([1:D-1]+(D-1)*(k-1),:);
    end
  end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Does the sign swapping if all image points are defined (ie, no nans are in u).
function [P,X] = signsPX_from_x(P,X,x)

[D N] = size(X);
K = size(P,1)/(D-1);

s = sign( reshape( sum(reshape(x,D-1,K*N).*reshape(P*X,D-1,K*N)), K,N) );

sP = s(:,1);
s = sP(:,ones(1,N)) .* s;

sX = s(1,:);
s = sX(ones(1,K),:) .* s;

if any(s(:)<0)
  P = [];
  X = [];
  return
end

aux = sP(:,ones(1,D-1))'; aux = aux(:);  P = P .* aux(:,ones(1,D));
X = X .* sX(ones(1,D),:);

return


% Does sign swapping if there are undefined points (nans) in x.
function [P,X] = signsPX_from_x_occl(P,X,u)

[D N] = size(X);
K = size(P,1)/(D-1);

PX = reshape( dot(reshape(u,D-1,K*N),reshape(P*X,D-1,K*N)), K,N);
for initp = 1:K
  p = NaN*ones(K,1);
  x = NaN*ones(1,N);
  p(initp) = 1;
  n = 1; 
  oldn = 0;
  while n-oldn > 0
    oldn = n;
    x = updatej(p,x,PX);
    p = updatej(x',p',PX')';
    n = nnz(~isnan([p' x]));
  end
  if all(all(isnan(PX) | ~isnan((p*x).*PX)))
    break
  end
end
P = (reshape(ones(D-1,1)*p',(D-1)*K,1)*ones(1,D)) .* P;
X = (ones(D,1)*x) .* X;

% check if the sign changing process was succesful
uPX = reshape( dot(reshape(u,D-1,K*N),reshape(P*X,D-1,K*N)), K,N);
if ~all(all( (uPX>0) | isnan(uPX) ))
  P = [];
  X = [];
end

return


% Auxiliary function used for swaping the signs of P^k, X_n.
function j = updatej(i,j,IJ)
mj = any(~isnan(i*ones(1,length(j))) &  isnan(ones(length(i),1)*j) & ~isnan(IJ));
nj = IJ(:,mj) .* (i*ones(1,nnz(mj)));
nj_ = nj(:); nj_(isnan(nj_)) = 0; nj(:) = nj_;
j(mj) = sign(sum(nj));
return
