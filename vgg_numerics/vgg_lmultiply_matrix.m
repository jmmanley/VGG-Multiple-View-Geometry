function M = vgg_lmultiply_matrix(A, size_B)
% VGG_LMULTIPLY_MATRIX  Return M such that M*vec(B) = vec(A*B)
%      M = VGG_LMULTIPLY_MATRIX(A, size(B))
%
%      See also vgg_matrix_test, vgg_rmultiply_matrix

% Author: awf@robots.ox.ac.uk

if isfinite(size_B(1))
  if size_B(1) ~= size(A,2)
    error
  end
end

if ~issparse(A)
  E = eye(size_B(2));
else
  E = speye(size_B(2));
end

M = kron(E, A);
