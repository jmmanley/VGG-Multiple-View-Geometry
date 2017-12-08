%F = vgg_F_from_P(P)  Compute fundamental matrix from two camera matrices.
%   P is cell (2), P = {P1 P2}. F has size (3,3). It is x2'*F*x1 = 0
%
%   Overall scale of F is unique and such that, for any X, P1, P2, it is
%   F*x1 = vgg_contreps(e2)*x2, where
%   x1 = P1*X, x2 = P2*X, e2 = P2*C1, C1 = vgg_wedge(P1).

function F = vgg_F_from_P(P, P2)

if nargin == 1
  P1 = P{1};
  P2 = P{2};
else
  P1 = P;
end

X1 = P1([2 3],:);
X2 = P1([3 1],:);
X3 = P1([1 2],:);
Y1 = P2([2 3],:);
Y2 = P2([3 1],:);
Y3 = P2([1 2],:);

F = [det([X1; Y1]) det([X2; Y1]) det([X3; Y1])
     det([X1; Y2]) det([X2; Y2]) det([X3; Y2])
     det([X1; Y3]) det([X2; Y3]) det([X3; Y3])];

return