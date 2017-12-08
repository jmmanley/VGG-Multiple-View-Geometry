function h = vgg_vech(m)
% VGG_VECH   Vectorization ("flattening") of symmetric matrix.
%
% For square matrix X, vgg_vech(X) is the column vector of elements on or
% below the main diagonal of m.
%
% Also works inversely: for a N*(N+1)/2-vector x, it returns symmetric
% N-by-N matrix X = vgg_vech(x) such that vgg_vech(X) = x.
%
% Useful for solving linear matrix equations, see Magnus and Neudecker.
%
% See vgg_matrix_test, vgg_duplic_matrix, and also vgg_vech_swap,
% vgg_commut_matrix.

[M N] = size(m);

if M==1 | N==1
  N = (sqrt(8*M*N+1)-1)/2;
  r = (1:N)'*ones(1,N);
  c = r';
  h = zeros(N);
  h(find(c <= r)) = m;
  h = h+h'-diag(diag(h));  
else
  r = (1:M)'*ones(1,N);
  c = ones(M,1)*(1:N);
  h = m(find(c <= r));
end

return