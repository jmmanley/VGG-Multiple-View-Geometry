function M = vgg_lmultiply_matrix(A, size_B)
% VGG_LMULTIPLY_MATRIX  Return M such that M*vec(B) = vec(B*A)
%      M = VGG_RMULTIPLY_MATRIX(A, size(B))
%
%      See also vgg_matrix_test, vgg_lmultiply_matrix

% Author: awf@robots.ox.ac.uk

if isfinite(size_B(1))
  if size_B(2) ~= size(A,1)
    error
  end
end

if ~issparse(A)
  E = eye(size_B(1));
else
  E = speye(size_B(1));
end

M = kron(A', E);
