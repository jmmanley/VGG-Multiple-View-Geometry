% XY = vgg_line3d_XY_from_pm(L)  Converts Pluecker matrix 3d line to a pair of homogeneous 3d points.
%
% L ... double(4,4), skew-symmetric Pluecker matrix of 3d line
% XY ... double(4,2), pair of 3d points spanning the line, in homogen. coordinates
%
% XY are obtained by svd, their homogeneous vectors are mutually orthogonal.

% T.Werner

function XY = vgg_line3d_XY_from_pm(L)

[u,s,v] = svd(vgg_contreps(L),0);
XY = u(:,3:4);

return