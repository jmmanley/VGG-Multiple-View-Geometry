% vgg_selfcalib_qaffine  Upgrading projective to quasi-affine reconstruction.
%
% Given projective reconstruction [P,X] with correct signs of P and X
% (the output of vgg_signsPX_from_x), it finds homography H transforming
% [P,X] to quasi-affine reconstruction [Pq,Xq] = [P*inv(H),H*X].
% Let Ainf=[0 0 0 1] be plane at infinity, then [Pq,Xq] has the property that :-
%
%   Ainf * Xq > 0                (all scene points in front of plane at infty)
%   Ainf * vgg_wedge(Pq{k}) > 0  (all camera centers in front of plane at infty)
%
% H = vgg_selfcalib_qaffine(P,X), where
%   P ... cell(K) of double(3,4), camera matrices. K is number of cameras.
%     P also can be 3x4xK array.
%   X ... double(4,N), scene points in homog. coordinates.
%   H ... cell{I} of double(4,4), homographies upgrading [P,X] to quasi-affine reconstruction.
%     There can be 0, 1, or 2 solution classes (corresponding to I=0,1,2) :-
%       - I==0 ... no solution, ie [P,X] cannot be transformed to any affine scene.
%       - I==1 ... 1 solution, ie camera centers and scene points
%           are not separable by a plane in the true scene.
%       - I==2 ... 2 solutions, ie camera centers and scene points are separable
%           by a plane in the true scene. Then there are two solutions for plane at infinity,
%           differing by sign(det(H{i})). The two reconstruction corresponding to H{1} and H{2}
%           have oppposite handedness and we cannot say which handedness is that of the true scene.
% (Note: by 'solution' we mean rather 'class of solutions' - indeed there are infinitely many
% solutions if I>0, and linear programming chooses a single solution out of them.)
%
% EXAMPLE: Let [P,X] be a projective reconstruciton from homogeneous image points x.
% Upgrade to quasi-affine reconstruction is done as follows:
%   [P,X] = vgg_signsPX_from_x(P,X,x);
%   H = vgg_selfcalib_qaffine(P,X);
%   H = H{1}; % single solution assumed
%   P = P*inv(H);
%   X = H*X;
% If either of rows 1 and 2 returns no solution, there's something wrong with
% the reconstruction, eg an outlier.

% T.Werner, Feb 2002, werner@robots.ox.ac.uk

function H = vgg_selfcalib_qaffine(P,X)

if ndims(P)==3
  for k = 1:size(P,3)
    Q{k} = P(:,:,k);
  end
  P = Q;
end

[D N] = size(X);
K = length(P);

for k = 1:K
  C(:,k) = vgg_wedge(P{k}); % oriented camera centers
end

% Solve chiral equalities:
%
% A := found plane at infinity
% detH := required det(H)
% (A and detH can be none, one, or two according to the number of solution classes)
detH = [];
A = [];
for detHa = [-1 1]
  Aa = sephplane([X detHa*C]);
  if ~isempty(Aa)
    A = [A; Aa];
    detH = [detH; detHa];
  end
end

if isempty(A)
  H = {};
  return
end


% compose final homography H
for i = 1:size(A,1)

  % find H{i} such that H{i}(4,:)==A
  [dummy,dummy,H{i}] = svd(A(i,:),0);
  H{i} = H{i}(end:-1:1,:);

  % make det(H{i}) the same sign as detH(i)
  if det(H{i})*detH(i) < 0
    H{i} = H{i}([2 1 3 4],:);
  end

  % 'beautifier' of X: 
  % Do singular value equalization on the set X,
  % i.e., make mean(nhom(H{i}*X),2)==[0;0;0] and svd(nhom(H{i}*X))==[1;1;1].
  Xi = vgg_get_nonhomg(H{i}*X);
  c = mean(Xi,2); % centroid
  Xi = Xi - c*ones(1,N);
  [U,S] = eig(Xi*Xi'); % sv equalization
  S = diag(1./sqrt(diag(S)));
  K = S*U';
  if det(K) < 0 % we want the sv equalization to be parity-preserving
    K = -K;
  end
  H{i} = [ K -K*c; 0 0 0 1 ]*H{i};
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% A = sephplane(X)  Finds separating hyperplane A such that all(A*X)>0.
% If no solution exists, A = [].
% Works for any dimension of X.
function A = sephplane(X)

[D,N] = size(X);

X = X ./ (ones(D,1)*sqrt(sum(X.*X)));
A = [-X' ones(N,1)];
b = zeros(size(A,1),1);
f = [zeros(1,D) -1]';
LB = [-ones(1,D) 0];
UB = [ones(1,D) Inf];
fprintf('vgg selfcalib_qaffine: linprog for %d %dd pts ... ', size(X,2), D);
options = optimset('linprog');
options.Display = 'off';
[res,FVAL,EXITFLAG] = linprog(f,A,b,[],[],LB,UB, [], options);
if isempty(res)
  fprintf('no feasible plane\n');
  A = [];
  return
end
A = res(1:D)';

if ~all(A*X > 0)
  fprintf('feasible plane returned, but is not in fact feasible\n');
  A = [];
end

fprintf('Got plane [%.2f %.2f %.2f %.2f]\n', A);

return


%i = k2i(k)
% Computes indices of joint point matrix rows corresponding to views k.
%% function i = k2i(k,step)
%% k = k(:)';
%% i = [1:3]'*ones(size(k)) + 3*(ones(3,1)*k-1);
%% i = i(:);
