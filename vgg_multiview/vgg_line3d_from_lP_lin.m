% vgg_line3d_linear  Linear estimation of 3d line from image lines and camera matrices.
%
% SYNOPSIS
% L = vgg_line3d_from_lP_lin(s,P [,imsize]), where
%
% s ... cell(K) of double(3,3), inv. covariance matrices of the K image line segments:-
%   - If the segments are estimated from edges, it is s(:,k) = x*x',
%     where x (3-by-N) are homog. coordinates of the edgels with last components 1.
%   - If only end points are available, s(:,k) = d*x*y' where x, y (column 2-vectors)
%     are the segment's end points and d its length.
%
% P ... cell(K) of double(3,4), camera matrices
%
% imsize ... double(2,K), image sizes (for preconditioning).
%   Omit if s and P are already preconditioned.
%
% L ... double(4,2), 3d straight line; columns of L are two homogeneous 
%   points spanning the line.

function L = vgg_line3d_from_lP_lin(s,P,imsize)

if nargin < 3, imsize = []; end

K = length(P); % number of images

% l := straight lines in homog. coords
for k = 1:K
  l(k,:) = vgg_fit_hplane_to_x(s{k});
end

if ~isempty(imsize) & K>2 % Preconditioning; for 2 views is not needed
  for k = 1:K
    [H,invH] = vgg_conditioner_from_image(imsize(:,k));
    P{k} = H*P{k};
    l(k,:) = l(k,:)*invH;
  end
  l = norml(l);
end

M = [];
for k = 1:K
  M = [M; l(k,:)*P{k}];
end

[u,s,v] = svd(normx(M')',0);
L = v(:,end-1:end);

return


function x = normx(x)
x = x./(ones(size(x,1),1)*sqrt(sum(x.*x)));

function l = norml(l)
% l = norml(l)  Multiplies hyperplane l by scalar so that for each n, norm(l(1:end-1,n))==1. 
l = l./(sqrt(sum(l(:,1:end-1).^2,2))*ones(1,size(l,2)));
