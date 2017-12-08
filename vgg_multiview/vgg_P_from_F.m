%P = vgg_P_from_F(F)  Compute cameras from fundamental matrix.
%   F has size (3,3), P has size (3,4).
%
%   If x2'*F*x1 = 0 for any pair of image points x1 and x2,
%   then the camera matrices of the image pair are 
%   P1 = eye(3,4) and P2 = vgg_P_from_F(F), up to a scene homography.

% Tomas Werner, Oct 2001

function P = vgg_P_from_F(F)

[U,S,V] = svd(F);
e = U(:,3);
P = [-vgg_contreps(e)*F e];

return