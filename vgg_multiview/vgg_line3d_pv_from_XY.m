% L = vgg_line3d_pv_from_XY(X,Y)  Pluecker vector 3d line from two 3d points.
%
% Syntax: L = vgg_line3d_pv_from_XY(X,Y) or 
%         L = vgg_line3d_pv_from_XY(XY)
% 
% X, Y ... size (4,N), 3D points
% XY ... size (4,2*N), XY stacked as XY = [X1 Y1 ... XN YN].
% L ... size (N,6), Pluecker vector(s) of 3D line(s).
%
% It is vectorized for N>1.
%
% See vgg_line3d_*

% T.Werner

function L = vgg_line3d_pv_from_XY(X,Y)

if nargin == 1
  Y = X(:,2:2:end);
  X = X(:,1:2:end);
end

if size(X,2) == 1
  L = [(X(1:3)*Y(4)-Y(1:3)*X(4))' vgg_wedge([X(1:3) Y(1:3)])];
else % vectorized version
  L = [(X(1:3,:).*([1;1;1]*Y(4,:))-Y(1:3,:).*([1;1;1]*X(4,:)))' vgg_wedge(X(1:3,:),Y(1:3,:))];
end

return