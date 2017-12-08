function vgg_matrix_test
% Test the matrix "flattening" functions from Magnus and Neudecker
%  vgg_vec
%  vgg_vech
%  vgg_commut_matrix
%  vgg_duplic_matrix
%  vgg_lmultiply_matrix

% Author: awf@robots.ox.ac.uk

A = randn(4,5);
B = randn(5,2);
C = randn(2,4);

fprintf('vgg_matrix_test: BEGIN\n');
assert('vgg_vec(A*B*C)', 'kron(C'',A)*vgg_vec(B)');
assert('vgg_vec(A*B)', 'vgg_lmultiply_matrix(A,size(B))*vgg_vec(B)');
assert('vgg_vec(B*C)', 'vgg_rmultiply_matrix(C,size(B))*vgg_vec(B)');


% Now try to solve A X + X B = C:
X_true = randn(4,4);
A = randn(4,4);
B = randn(4,4);
C = A*X_true + X_true*B;

% Re-express in terms of x = vec(X):
% fA x + gB x  = vec(C)
% so
% x = (fA + gB) \ vec(C)
size_X = [4 4]; % need to know size(X) to convert the matrix multiplys
fA = vgg_lmultiply_matrix(A, size_X);
gB = vgg_rmultiply_matrix(B, size_X);
x = (fA + gB)\vgg_vec(C);

% Now check we got it right...
assert('x', 'vgg_vec(X_true)')



fprintf('vgg_matrix_test: END\n');

function assert(sx,sy)
x = evalin('caller', sx);
y = evalin('caller', sy);
err = norm(x-y);
if err > 1e-12
  dbstack
  fprintf('vgg_matrix_test: FAILED %s = %s\n', sx, sy);
  keyboard
else
  fprintf('vgg_matrix_test: PASSED %s = %s\n', sx, sy);
end
