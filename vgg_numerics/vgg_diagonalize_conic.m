% H = vgg_diagonalize_conic(C)  Finds Euclidean transformation sending conic canonical position.
%
% For any symmetric matrix C, returns Euclidean transformation H such that 
%    H'*C*H
% is a diagonal matrix.
%
% Typical usage is to transform conics to canonical form, to classify or plot them.

function H = vgg_diagonalize_conic(C)

N = size(C,1);

C = C/C(end,end);
s = C(1:end-1,end);
Q = C(1:end-1,1:end-1);

[U,S,V] = svd(Q);
sw = diag([-1 ones(1,N-2)]);
if det(U) < 0
  U = U*sw;
  S = sw*S;
end
if det(V) < 0
  V = V*sw;
  S = S*sw;
end

H = [U -inv(Q)*s; zeros(1,N-1) 1];

return