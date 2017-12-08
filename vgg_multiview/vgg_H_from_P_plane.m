% vgg_H_from_P_plane  Given 3D plane and camera, returns 4-by-3 matrix mapping image points to points on the plane.
%
% H = vgg_H_from_P_plane(A,P), where
%   A ... size (4,1), scene plane
%   P ... size (3,4), camera matrix
%   H ... size (4,3), matrix such that X = x*H, where x (size (3,1)) is an image point
%     and X is a scene point lying on the scene plane A.
%     In other words, P*H = s*eye(3), where s is a scalar factor.
% 
% The overall sign of H is chosen such that s*A*vgg_wedge(P) > 0.

% T. Werner, Oct 2001

function H = vgg_H_from_P_plane(A,P)

A = A';
p1 = P(1,:);
p2 = P(2,:);
p3 = P(3,:);

H = [vgg_contreps(p2'*p3 - p3'*p2)*A ...
     vgg_contreps(p3'*p1 - p1'*p3)*A ...
     vgg_contreps(p1'*p2 - p2'*p1)*A];

return
