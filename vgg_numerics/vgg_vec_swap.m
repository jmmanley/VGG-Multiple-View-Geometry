% S = vgg_vec_swap(x,y)  For square matrix A and vectors x, y it is x'*A*y = vgg_vec(A)'*vgg_vec_swap(x,y). 
%
% x ... matrix N-by-K
% y ... matrix N-by-K
% S ... matrix K^2-by-N
%
% Examples :-
%
%  - Estimating fundamental matrix F that should satisfy x(:,k)'*F*y(:,k) = 0 from
%    given points x(:,k), y(:,k). We solve homogeneous system vec(F)*vgg_vec_swap(x,y) = 0.
%
%  - Evaluating values v(k) of x(:,k)'*A{k}*y(:,k) for all k, A{k} are square matrices.
%    If column vectors vec(A{k}) are stacked to a single matrix AA, then
%        v = sum(AA'.*vgg_vec_swap(x,y)).
%
% See also vgg_vech_swap

% Tomas Werner, Oct 2001

function M = vgg_vec_swap(x,y)

if nargin==1
  y = x;
end

[i j] = find(ones(size(x,1)));
M = x(i,:).*y(j,:);

return

