% XY = vgg_line3d_XY_from_pv(L)  Converts Pluecker vector 3d line to a pair of homogeneous points.
%
% L ... double(1,6), Pluecker vector
% XY ... double(4,2), pair of 3d points spanning the line
%
% XY are obtained by svd, their homogeneous vectors are mutually orthogonal.

% T.Werner

function XY = vgg_line3d_XY_from_pv(L)

XY = vgg_line3d_XY_from_pm(vgg_line3d_pm_from_pv(L));

return