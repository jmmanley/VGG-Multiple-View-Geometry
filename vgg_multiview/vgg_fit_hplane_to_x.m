% vgg_fit_hplane_to_x  Fitting hyperplane to set of points.
%
% SYNOPSIS
% [A,e] = vgg_fit_hplane_to_x(s), where
%
%   s ... double(N+1,N+1), inv. covariance matrix s = x'*x, the columns of which
%      are (N+1)-vectors of homog. coordinates of points in N-space (s(end,:)==1)
%   A ... double(1,N+1), fitted hyperplane (homog. coordinates).
%     Homogeneous part of A is normalized, norm(A(1,1:N))==1.
%   e ... double(N,1), eignevalues of the fit
%
% Hyperplane A minimizes sum of squared orthogonal distances of points to it.
%
% EXAMPLE
% Let e be (2,?) vector of detected edgels in nonhomog. coords. Then fitting a straight
% line l to this set of edgels is:
%   e(end,:) = 1;
%   l = vgg_fit_hplane_to_x(e*e');

function [A,e] = vgg_fit_hplane_to_x(s)

N = size(s,1) - 1; % space dimension
c = s(1:N,end)/s(end,end); % centroid

[U,e,V] = svd(s(1:N,1:N)-c*s(end,1:N),0);

A = U(:,end)';
A = [A -A*c];
e = diag(e);

return
