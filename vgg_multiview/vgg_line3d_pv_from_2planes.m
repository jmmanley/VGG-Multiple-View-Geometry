% L = vgg_line3d_pv_from_2planes(A,B) Pluecker vector of 3d line met by two planes.
%
% Syntax: L = vgg_line3d_pv_from_2planes(A,B) or 
%         L = vgg_line3d_pv_from_2planes(AB)
%
% A, B ... size (N,4), 3d planes
% AB ... size (2*N,4), N pairs of 3d planes. It is AB = [A1; B1; A2; B2; ... ; AN; BN].
% L ... size (N,6), line(s) in Pluecker vector.
%
% It is vectorized for N>1.

% T.Werner

function L = vgg_line3d_pv_from_2planes(A,B)

if nargin == 1
  L = vgg_line3d_pv_from_XY(A');
else
  L = vgg_line3d_pv_from_XY(A',B');
end
L = L(:,[4:6 1:3]);
 
return