% S = vgg_vech_swap(x,y)  For symmetric matrix A it is x'*A*y = vgg_vech(A)'*vgg_vech_swap(x,y). 
%
% x ... matrix N-by-K
% y ... matrix N-by-K
% S ... matrix K*(K+1)/2-by-N
%
% Examples :-
%
%  - Estimating symmetric matrix A that should satisfy x(:,k)'*A*y(:,k) = 0 from
%    given points x(:,k), y(:,k). We solve homogeneous system vech(A)*vgg_vec_swap(x,y) = 0.
%
%  - Evaluating values v(k) of x(:,k)'*A{k}*y(:,k) for all k, A{k} are square symmetric matrices.
%    If column vectors vech(A{k}) are stacked to a single matrix AA, then
%        v = sum(AA'.*vgg_vech_swap(x,y)).
%
% See also vgg_matrix_test, vgg_lineseg_from_x, vgg_vec_swap, vgg_commut_matrix, vgg_duplic_matrix.

% Tomas Werner, Oct 2001

function M = vgg_vech_swap(x,y)

if nargin==1
  y = x;
end

[i j] = find(ones(size(x,1)));
d = i>=j;
i = i(d);
j = j(d);

M = x(i,:).*y(j,:) + x(j,:).*y(i,:).*((i~=j)*ones(1,size(x,2)));

return