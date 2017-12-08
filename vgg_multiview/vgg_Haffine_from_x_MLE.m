function H = vgg_Haffine_from_x_MLE(xs1,xs2)
% H = vgg_Haffine_from_x_MLE(xs1,xs2)
%
% Compute MLE for affine H, i.e. find H and xhat1 such that
% d^2(xs1,xhat1) + d^2(xs2,xhat2) minimized where xhat2 is affine transf of xhat1.
%
% Parameters:
%   xs1, xs2 ... double(3,N), N pairs of corresponding points (homogeneous)
%   H ... double(3,3), affine transformation
%
% See HZ page 115 1st edition, page 130 2nd edition
% az 17/11/2001

if any(size(xs1) ~= size(xs2))
 error ('Input point sets are different sizes!')
end

% condition points

nonhomg = vgg_get_nonhomg(xs1);
means = mean(nonhomg');
maxstds = max(std(nonhomg'));
C1 = diag([1/maxstds 1/maxstds 1]);  % only similarity 
C1(:,3) = [-means/maxstds 1]';

nonhomg = vgg_get_nonhomg(xs2);
means = mean(nonhomg');
C2 = C1;            % nb must use same scaling for both point sets
C2(:,3) = [-means/maxstds 1]';

xs1 = vgg_condition_2d(xs1,C1);
xs2 = vgg_condition_2d(xs2,C2);

% NB conditioned points have mean zero, so translation
% part of affine transf is zero 2-vector

xs1nh = vgg_get_nonhomg(xs1);
xs2nh = vgg_get_nonhomg(xs2);

A = [xs1nh;xs2nh]';

% Extract nullspace
[u,s,v] = svd(A); 
s = diag(s);
 
nullspace_dimension = sum(s < eps * s(2) * 1e3);
if nullspace_dimension > 2
  fprintf('Nullspace is a bit roomy...');
end

% compute affine matrix from two largest singular vecs

B = v(1:2,1:2);
C = v(3:4,1:2);

H = [ C * pinv(B) , zeros(2,1); 0 0 1];

%decondition
H = inv(C2) * H * C1;

H = H/H(3,3);
