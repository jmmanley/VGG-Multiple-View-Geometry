% vgg_mrdivs  Solves equation system Y*diag(s) = A*X with unkowns A, s.
%
% A = vgg_mrdivs(X,Y) solves (overdetermined) equation system Y*diag(s) = A*X
% by linear method (DLT algorithm).
% Parameters:
%   X ... double (N,K)
%   Y ... double (M,K)
%   A ... double (M,N)
%   s ... double (1,K)
%
% Preconditioning of the points not included in the function. Use vgg_conditioner_*.
%
% Typical usage:
%   1. Estimating an image homography from K pairs of corresponding points.
%      If 3-by-K matrices x and y are the points in homogeneous coordinates, the 3-by-3 homography
%      matrix is obtained as H = vgg_mrdivs(x,y).
%
%   2. Estimating 3-by-4 camera projection matrix P from corresponding pairs of image and scene points.
%      For image points x (3xK matrix) and scene points X (4xK matrix) do P = vgg_mrdivs(X,x).

% (c) werner@robots.ox.ac.uk

% Algorithm: 
% 
% For n-th point pair X(:,n) and Y(:,n) we have
%  s*X(:,n) = A*Y(:,n)
% We eliminate s what results in MY*(MY-1)/2 linear homogenous equations 
% for elements of A. We solve this system by svd or eig.

function A = vgg_mrdivs(X,Y)

[MX,N] = size(X);
[MY,NY] = size(Y);
if N ~= NY, error('Matrices A, B must have equal number of columns.'); end

 % construct the measurement matrix
W = zeros(MX*MY,MY*(MY-1)/2*N);
k = 1;
for i = 1:MY
  for j = 1:i-1
    W([[1:MX]+MX*(j-1) [1:MX]+MX*(i-1)],[1:N]+N*(k-1)) =  [(ones(MX,1)*Y(i,:)).*X; -(ones(MX,1)*Y(j,:)).*X];
    k = k+1;
  end
end

% solve the system || A'(:)' * W || ---> min
[dummy,s,A] = svd(W',0);
A = reshape(A(:,end),MX,MY)';

return