% Q = vgg_line3d_Ppv(P)  Transforms a camera matrix to use it with Pluecker vector 3d line representation.
%
% P ... double(3,4), ordinary camera matrix
% Q ... double(6,3), transformed camera matrix (quadratic function of elements of P)
%   such that the projection of a 3d line by the camera is given by
%
%     l = Q*L
%
% where l is 1-by-3 image line vector, and L is 1-by-6 Pluecker vector of the 3d line.

% T.Werner

function Q = vgg_line3d_Ppv(P)
 
K = size(P,1)/3;
i1 = 1:3:3*K;
i2 = 2:3:3*K;
i3 = 3:3:3*K;

Q(:,i1) = vgg_line3d_pv_from_XY(P(i2,:)',P(i3,:)')';
Q(:,i2) = vgg_line3d_pv_from_XY(P(i3,:)',P(i1,:)')';
Q(:,i3) = vgg_line3d_pv_from_XY(P(i1,:)',P(i2,:)')';

return