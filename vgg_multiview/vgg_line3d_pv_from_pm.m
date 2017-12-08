% L = vgg_line3d_pv_from_pm(G)  Conversion of Pluecker matrix to Pluecker vector 3d line representation.
%
% G ... double(4,4), skew-symmetric Pluecker matrix
% L ... double(1,6), Pluecker vector

% T.Werner

function L = vgg_line3d_pv_from_pm(G)

L = [G(1:3,4); vgg_contreps(G(1:3,1:3))']';

return