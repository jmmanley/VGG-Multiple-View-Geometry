% VGG_COMMUT_MATRIX  Commutation matrix, to transpose a matrix(:)
%
% Classical matrix re-arrangement operator, see book Magnus-Neudecker.
%
% Useful for rearranging matrix equations. It is
%  vgg_vec(X') = vgg_commut_matrix(size(X))*vgg_vec(X).
%
% See also vgg_matrix_test

% Added by Tom Werner, originally from Tomas Minka.

function k = vgg_commut_matrix(n, m)

if nargin < 2
  m = n(2);
  n = n(1);
end

k = reshape(kron(vgg_vec(eye(n)), eye(m)), n*m, n*m);

return
