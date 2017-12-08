% vgg_selfcalib_metric_vansq  Metric selfcalibration from 3 orthogonal principal directions and square pixels.
%
% DESCRIPTION
% Given projective camera matrices P and 3 scene points V, it computes 3D-to-3D
% homography H which upgrades the old reconstruction to metric one, i.e.,
% differing from the true one only by isotropic scaling.
% H is computed from the following constraints :-
%   (1) points H*V are at infty and mutually orthogonal (i.e., H*V==eye(4,3)),
%   (2) cameras P*inv(H) have square pixels,
% Constraint (1) is hard, (2) is soft. I.e., if [P,V] are not consistent
% (2) will be satisfied only partially (in linear least squares sense).
%
% SYNSOPSIS
% [H,sv] = vgg_selfcalib_metric_vansq(P,V), where
%   P ... cell{K} of double(3,4), projective cameras
%   V ... double(4,3), 3 projective scene points (homog. coords.)
%   H ... double(4,4), upgrading 3D-to-3D homography
%   sv ... 2-vector, last 2 singular values of linear system solving for square pixels.
%     In healthy situation, sv(2) must be tiny and sv(1) reasonably large.
%
% NOTE: To get correct handedness (= non-mirroring) of the reconstruction,
% make sure that V satisfies
%   vgg_wedge(V)*[X C] > 0
% for all scene points X and camera centers C(:,k) = vgg_wedge(P{k}). This can be
% achieved by swapping signs of P and X using vgg_signsPX_from_x. The thing
% requires also positive handedness of image and scene coord. systems.
%
% SEE ALSO vgg_selfcalib_qaffine.

% T.Werner, Feb 2002

function H = vgg_selfcalib_metric_vansq(P,V)

K = length(P);

%%%%%%%%%%
% Step 1: 
%   Find 3D homography H1 sending V to eye(4,3).
%   This results in reconstruction differing from metric one only in scaling in axis directions.
%%%%%%%%%%%


H1 = inv(normx([V vgg_wedge(V)']));
for k = 1:K
  P{k} = P{k}/H1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2:
%   Find scaling in axis directions from square pixel assumption. 
%   This scaling is represented by diagonal 3D homography H2.
%
% Algorithm: 
% Consider input camera matrix P = [M v] and output matrix Q = K*R*[eye(3) -t]. It is Q =~ P*H where H = diag([d 1]). 
% Let O = diag([1 1 1 0]) be absolute quadric. Then Q*O*Q' = DIAC = inv(IAC). For square pixels, IAC(1,1)=IAC(2,2) and IAC(1,2)=0.
% Substitution gives Q*O*Q' =~ P*H*O*H'*P' where P*H*O*H'*P' = M*diag(d.^2)*M' =~ inv(IAC). I.e., 
%   IAC =~ inv(M)*diag(d.^(-2))*inv(M)'.
% The RHS of the last expression can be rearranged as 
%   vech(IAC) =~ pinv(duplication(3))*kron(inv(M)',inv(M)')*diagonalize(3)*d, where size(d)=[3 1].
% This is used to build a system of linear equations. We showed things only for onen camera - more cameras can be added to the system easily.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compose the linear system
A = [];
for k = 1:K
  aux = inv(P{k}(:,1:3))';
  aux = pinv(vgg_duplic_matrix(3))*kron(aux,aux)*diagonalize(3);
  A = [A; [aux(1,:)-aux(4,:); aux(2,:)]];
end

% solve it
[dummy,sv,d] = svd(A,0);
d = d(:,end);
sv = diag(sv);
sv = sv(2:3)/sv(1); % normalize sing values

% form H
d = 1./sqrt(abs(d));
H2 = diag([d;1]);
H = inv(H2)*H1;

return

%%%%%%%%%%%%%%%%%%%%


% G = diagonalize(n)  Diagonalization matrix. It is vec(diag(x)) = diagonalize(length(x))*x.
function G = diagonalize(n)
G = zeros(n^2,n);
i = [];
for j = 0:n-1
  i = [i 1+n*j+j];
end
G(i,:) = eye(n);
return


% x = normx(x)  Normalize MxN matrix so that norm of each its column is 1.
function x = normx(x)
if ~isempty(x)
  x = x./(ones(size(x,1),1)*sqrt(sum(x.*x)));
end
return