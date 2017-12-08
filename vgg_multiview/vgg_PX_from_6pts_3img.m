function [P,X] = vgg_PX_from_6pts_3img(x1,x2,x3)
% vgg_PX_from_6pts_3img  Computes camera matrices and world points 
% from 6 points across 3 images.
%
%   [P,X] = vgg_PX_from_6pts_3img(x), where
%      x ... double(3,6,3) or cell{3} of double(3,6), 6 homogeneous points in 3 images
%      P ... double(3,4,3), P(:,:,k) is k-th camera matrix
%      X ... double(4,6), homogeneous world points
%   There are 0 to 3 solutions for (P,X). Solutions are pruned by requirement that
%   scalars s in all projective equations s*x==P*X are positive.
%   In case of multiple solutions, P and X have one dimension
%   more such that P(:,:,:,n) and X(:,:,n) is the n-th solution.
%
%   Also the form [P,X] = vgg_PX_from_6pts_3img(x1,x2,x3) is accepted.

% Algorithm in Hartley-Zisserman, Alg 19.1 page 493 in 1st edition,
%                                 Alg 20.1 page 511 in 2nd edition
% Coded by werner@robots.ox.ac.uk, Nov 2002.

if nargin==3
  x = cat(3,x1,x2,x3);
else
  if iscell(x1)
    x = cat(3,x1{:});
  else
    x = x1;
  end
end
if any(size(x)~=[3 6 3])
  error('Wrong size of input points.');
end

for k = 1:3
  % Find homographies H_k mapping first 4 pts in each image to standard projective basis.
  % Now, x(:,1:4,k) = H(:,:,k)*[eye(3) [1;1;1]].
  H(:,:,k) = H_from_4x(x(:,3:6,k));
  
  % Form transformed points xs(:,1:2,k)
  xs(:,:,k) = inv(H(:,:,k))*x(:,:,k);
end

% Compute dual fundamental matrix
Fd = Fdual_from_x(xs(:,1:2,:));

% Retrieve (non-dual) cameras P and world points X from (each solution for) dual fund. matrix.
P = [];
X = [];
for i = 1:size(Fd,3)
  
  % Compute canonical Xi and Pi
  Xi = X_from_Fdual(Fd(:,:,i));
  Pi = P_from_Xx_canonical(Xi,xs);
  
  % return from canonical to original image bases
  for k = 1:3
    Pi(:,:,k) = H(:,:,k)*Pi(:,:,k);
  end
  
  % Compute signs of P and X, if possible. If impossible, Pi==Xi==[].
  [Pi,Xi] = vgg_signsPX_from_x(Pi,Xi,x);
  if isempty(Pi), continue, end
  %for k = 1:3, Pi(:,:,k)*Xi ./ x(:,:,k), end  % test code

  P = cat(4,P,Pi);
  X = cat(3,X,Xi);

end

return


%%%%%%%%%%%%%%%% auxiliary functions


% Solve for dual fundamental matrix.
function F = Fdual_from_x(x)

% Linear step. After that,
% for k=1:3 and i=1:2, it is x(:,1,k)'*F_i*x(:,2,k)==0.
A = [];
for k = 1:3
  x1 = x(1,1,k);  y1 = x(2,1,k);  z1 = x(3,1,k);
  x2 = x(1,2,k);  y2 = x(2,2,k);  z2 = x(3,2,k);
  A = [A; [x1*y2-z1*y2,...
           x1*z2-z1*y2,...
           y1*x2-z1*y2,...
           y1*z2-z1*y2,...
           z1*x2-z1*y2] ];
end
[u,s,v] = svd(A,0);
v1 = v(:,end-1);
v2 = v(:,end);
FF{1} = [0 v1(1:2)'; v1(3) 0 v1(4); v1(5) -sum(v1) 0];
FF{2} = [0 v2(1:2)'; v2(3) 0 v2(4); v2(5) -sum(v2) 0];

% Non-linear step. Find linear combination of F_i having zero determinant.
% Dual fund. matrix is now  F = a*FF{1} + (1-a)*FF{2}, for each element of a.
a = vgg_singF_from_FF(FF);
for i = 1:length(a)
  F(:,:,i) = a(i)*FF{1} + (1-a(i))*FF{2};
end
return


% Given dual fundamental matrix Fd, computes cameras P and world points X such that
% for each k, P(:,:,k)*X ~ x(:,:,k).
function X = X_from_Fdual(Fd)

% Retrieve second dual camera matrix Pd2.
% It is done by considering Fdual in form
%[0 b*(d-c) -c*(d-b)
% -a*(d-c) 0 c*(d-a)
% a*(d-b) -b*(d-a) 0].

Fd_aux = reshape( Fd([4 7 1 2 1 8 1 3 6]), [3 3] );
[u,s,v] = svd(Fd_aux,0);
ABC = v(:,3); % It is a:b:c = A:B:C
[u,s,v] = svd(Fd',0);
KLM = v(:,3); % It is K:L:M = ((d-a):(d-b):(d-c)

G = [0 -ABC(3) ABC(2) 0
     ABC(3) 0 -ABC(1) 0
     -ABC(2) ABC(1) 0 0
     KLM(2) -KLM(1) 0 KLM(1)-KLM(2)
     0 KLM(3) -KLM(2) KLM(2)-KLM(3)
     -KLM(3) 0 KLM(1) KLM(3)-KLM(1) ];
[u,s,v] = svd(G,0);
abcd = v(:,end);

% The test code:
%Pd1 = [eye(3) [1;1;1]];
%Pd2 = [diag(abcd(1:3)) [1;1;1]*abcd(4)];
%vgg_F_from_P(Pd1,Pd2) ./ Fd

% Retrieve 6 world points X(:,1:6)
X = [ abcd [1;1;1;1] eye(4) ];

return


% Computes camera P from world points X and image points x, everything in canonical form.
function P = P_from_Xx_canonical(X,x)
X = X(:,1);
for k = 1:3
  A = [contreps(x(:,1,k))*[X(1) 0 0 X(4)
                           0 X(2) 0 X(4)
                           0 0 X(3) X(4)]
       contreps(x(:,2,k))*[eye(3) [1;1;1]]];
  [u,s,v] = svd(A,0);
  a = v(:,end);
  P(:,:,k) = [a(1) 0 0 a(4)
              0 a(2) 0 a(4)
              0 0 a(3) a(4)];
end
return


% A bit faster vgg_contreps.
function X = contreps(x)
X = [0 x(3) -x(2)
     -x(3) 0 x(1)
     x(2) -x(1) 0];
return


% H_from_4x  Having four point matches, the homography relating them is given by
% H = H_from_4x(x2)*inv(H_from_4x(x1)).
function H = H_from_4x(x)
H = x(:,1:3) * diag(inv(x(:,1:3))*x(:,4));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Test code:
P = randn(3,4,3); X = randn(4,6);
for k = 1:3, x(:,:,k) = P(:,:,k)*X; end
[Pn,Xn] = vgg_PX_from_6pts_3img(x);

% check image points x predicted up to scale
for k=1:3, Pn(:,:,k)*Xn(:,:) ./ x(:,:,k), end

