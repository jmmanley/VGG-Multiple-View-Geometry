% vgg_vec  De-/vectorization of a matrix.
%
% For a matrix X, vgg_vec(X) = X(:).
% For a N^2-vector x, vgg_vec(x) = reshape(x,N,N).
%
% Classical matrix re-arrangement operator, see book Magnus-Neudecker.
% Trivial function, included mainly for consistency with notation in literature.
%
% Useful for rearranging matrix equations. Matrix from the middle of a product can be put to the right as
%
%    vgg_vec(A*B*C) = kron(C',A)*vgg_vec(B)
%
% See also vgg_vec_swap, vgg_commut_matrix, vgg_duplic_matrix.

function v = vgg_vec(A)

if any(size(A)==1)
  v = reshape(A,[1 1]*sqrt(prod(size(A))));
else
  v = A(:);
end
  
return