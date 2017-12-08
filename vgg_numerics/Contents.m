% VGG numerics library
%
% Linear equation system solving
%  vgg_mrdivs             - solving for matrix A from equation system y =~ A*x (equality up to scale)
%
% Rearranging matrix equations
%  vgg_vec                - vectorize matrix, same as (:) operator, and its inverse for square matrix
%  vgg_vech               - as vec, takes only entries on and under diagonal
%  vgg_commut_matrix       - it is vec(X') = commut_matrix(size(X))*vec(X)
%  vgg_duplic_matrix       - for square X, it is duplic_matrix(n)*vech(X) = vec(X)
%  vgg_vech_swap           - for square A, it is x'*A*y = vec(A)'*vec_swap(x,y), vectorized
%  vgg_vec_swap            - for square symmetric A, it is x'*A*y = vech(A)'*vech_swap(x,y), vectorized
%
% Misc
%  vgg_gauss_mask          - univariate gaussian or any its derivative
%  vgg_diagonalize_conic   - transforms conic to canonical position
%  vgg_wedge              - wedge (Grassmann outer) product of N-1 N-vectors
