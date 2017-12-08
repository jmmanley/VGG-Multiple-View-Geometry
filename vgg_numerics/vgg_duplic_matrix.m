% d = vgg_duplic_matrix(n)  Duplication matrix.
%
% Classical matrix re-arrangement operator, see book Magnus-Neudecker.
%
% Useful for rearranging equations with symmetric matrices. 
% For square symmetric X, it is
%
%   vgg_duplic_matrix(n)*vgg_vech(X) = vgg_vec(X)
%
% See also vgg_matrix_test, vgg_vech_swap.

% Added by Tom Werner, originaly from T.Minka.


function d = vgg_duplic_matrix(n)

a = tril(ones(n));
i = find(a);
a(i) = 1:length(i);
a = a + tril(a,-1)';
j = a(:);

m = n*(n+1)/2;
d = zeros(n*n,m);
for r = 1:size(d,1)
  d(r, j(r)) = 1;
end
