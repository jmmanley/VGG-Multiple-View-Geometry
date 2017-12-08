% vgg_line3d_from_lP_nonlin  Non-linear estimation of (possibly constrained) 3D line segment from image line segments.
%
% SYNOPSIS
% L = vgg_line3d_from_lP_nonlin(s,P [,imsize] [,L0] [,X] [,nonlin_opt])
%
% s ... cell(K) of double(3,3), inv. covariance matrices of the K image line segments:-
%   - If the segments are estimated from edges, it is s(:,k) = x*x',
%     where x (3-by-N) are homog. coordinates of the edgels with last components 1.
%   - If only end points are available, s(:,k) = d*x*y' where x, y (column 2-vectors)
%     are the segment's end points and d its length.
%
% P ... K-cell with 3-by-4 camera matrices
%
% imsize ... size (2,K), image size(s) for preconditiong.
%   Omit if s and P are already preconditioned.
%
% L0 ... double(4,2), initial scene line (optional). Homogeneous points L0(:,i) span 
%   the line. If omitted, linear estimation is done first.
%
% X ... constraint on L :-
%   - if X is omitted:     no constraint on L
%   - if X is double(4,1): L goes thru point X
%   - if X is double(4,2): L goes thru 3D line spanned by X(:,i)
%
% nonlin_opt ... options for Levenberg-Marquardt. It is comma-separated list
%   of pairs ['option',value]. Possible options are :-
% opt ... options structure with possible fields :-
%   - verbose ... 1 or 0 (default: 0)
%   - niter_term ... maximum number of iterations
%   - rmsstep_term ... terminating step of rms of residuals
%   - lambda_term ... terminating value of lambda (default: 1e10)
%   - lambda_init ... initial value of lambda
% E.g., vgg_line3d_from_lP_nonlin(...,'lambda_init',1e-9,'niter_term',5).
%
% L ... double(4,2), estimated 3D line. Points L(:,i) span the line.
%
% Note: use [] if you want to omit a parameter and use a later one, e.g.
%   vgg_line3d_from_lP_nonlin(s,P,imsize,[],[],'verbose',1,'lam_init',1e-9)
%
% ALGORITHM
% - Minimization is done by Levenberg-Marquardt.
% - 3D line L is parameterized by image lines in the first two images.
%   The positions of these image lines are possibly constrained by X.

% T.Werner, Feb 2002

function L = vgg_line3d_from_lP_nonlin(s,P,imsize,L,X,varargin)

if nargin < 3, imsize = []; end
if nargin < 4, L = []; end
if nargin < 5, X = []; end

if isempty(L)
  L = vgg_line3d_from_lP_lin(s,P,imsize);
end
K = length(P); % number of images
if K<2
  error('Cannot reconstruct 3D line from 1 image');
end
if isempty(X) & K==2 % no need for non-linear minimization
  return
end

% Prepare square root of covariance matrices; now s{k}(:,n) has meaning of 3 homogeneous image points
for k = 1:K
  [us,ss,vs] = svd(s{k},0);
  s{k} = us*sqrt(ss);
end

% Preconditioning
if ~isempty(imsize)
  for k = 1:K
    H = vgg_conditioner_from_image(imsize(:,k));
    P{k} = H*P{k};
    s{k} = H*s{k};
    scale(k) = H(1,1); % save the scales for evaluating objective function
  end
else
  scale = ones(1,K);
end

switch size(X,2)
 case 0 
  % Scene line L is unconstrained, having thus 4 DOF.
  % L is parameterized by two image lines in images 1 and 2, each having 2 DOF, as follows:
  %   l1 = l0(1,:) + p(1:2)'*ldelta{1}
  %   l2 = l0(2,:) + p(3:4)'*ldelta{2}
  % where row 4-vector p represents 4 DOF of L.
  for k = 1:2
    l0(k,:) = normx(vgg_wedge(P{k}*L)')';
    ldelta{k} = null(l0(k,:))';
  end

  % optimization
  p = levmarq(@F, {vertcat(P{:}),s,scale,l0,ldelta},...
              @normsolve,...
              [0;0;0;0],...
              varargin{:});
  l = l12_from_p(p,l0,ldelta);

 case 1 
  % Scene line L is constrained to intersect the scene point X, having thus 2 DOF.
  % L is parameterized by two image lines in images 1 and 2, each having 1 DOF, as follows:
  %   l1 = l0(1,:) + p(1)*ldelta{1}
  %   l2 = l0(2,:) + p(2)*ldelta{2}
  % where 2-vector p represents 2 DOF of L.
  for k = 1:2
    l0(k,:) = normx(vgg_wedge(P{k}*L)')';
    x = P{k}*X;
    
    % Since L might not intersect X, move l0(k,:) 'as little as possible' to intersect x.
    l0(k,:) = l0(k,:) - (l0(k,:)*x)/(x'*x).*x';
    
    Q = null(x')';
    ldelta{k} = null(l0(k,:)*pinv(Q))'*Q;
  end

  % optimization
  p = levmarq(@F, {vertcat(P{:}),s,scale,l0,ldelta},...
              @normsolve,...
              [0;0],...
              varargin{:});
  l = l12_from_p(p,l0,ldelta);

 case 2
  % Scene line L is constrained to intersect the scene line given by 2 points X.
  % This constraint is given by 
  %   l1*G*l2' = 0
  % where G is 3x3 rank 2 matrix (analogical in fact to fundamental matrix) and
  %   G = P{1}*M*P{2}'
  % where M is Pluecker matrix of line given by X, M = X(:,1)*X(:,2)'-X(:,2)*X(:,1).
  %
  % This constraint can be written as
  %   p(1:2)*D*p(3:4)' + p(1:2)*d{2}' + d{1}*p(3:4)' + c = 0
  % where D, d, c are given below and 4-vector p are 4 parameters of L.
  %
  % L is parameterized by two image lines in images 1 and 2, each having 2 DOF, as follows:
  %   l1 = l0(1,:) + p(1:2)'*ldelta{1}
  %   l2 = l0(2,:) + p(3:4)'*ldelta{2}
  % where p(1:3) are chosen freely and p(4) is computed from the above formula as
  %   p(4) = -(p(1:2)'*(D(:,1)*p(3)+d{1})+p(3)*d{2}(1)+c)/(p(1:2)'*D(:,2)+d{2}(2)).
  Lpm = X(:,1)*X(:,2)' - X(:,2)*X(:,1)';
  G = P{1}*Lpm*P{2}';
  for k = 1:2
    l0(k,:) = normx(vgg_wedge(P{k}*L)')';
  end  

  % As L might not intersect line X, move l0(2,:) 'as little as possible' to enforce l0(1,:)*G*l0(2,:)'==0.
  x = (l0(1,:)*G)';
  l0(2,:) = l0(2,:) - (l0(2,:)*x)/(x'*x).*x';
  
  for k = 1:2
    ldelta{k} = null(l0(k,:))';
  end  
  
  D = ldelta{1}*G*ldelta{2}';
  d{1} = ldelta{1}*G*l0(2,:)';
  d{2} = ldelta{2}*G'*l0(1,:)';
  c = l0(1,:)*G*l0(2,:)';
    
  % optimization
  p = levmarq(@F, {vertcat(P{:}),s,scale,l0,ldelta,D,d,c},...
              @normsolve,...
              [0;0;0],...
              varargin{:});
  l = l12_from_p(p,l0,ldelta,D,d,c);

end
if all(~isnan(l(:)))
  L = null([l(1,:)*P{1}; l(2,:)*P{2}]);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% objective function


% Objective function of Levenberg-Marquardt
function [y,w,J] = F(p,P,s,scale,varargin)
K = length(s);

% l := lines in images 1, 2 from p
l = l12_from_p(p,varargin{:});

% l := reprojection of lines in images 1, 2 to all images
if all(abs(p) < inf)
  [dummy,dummy,L]= svd([l(1,:)*P(1:3,:); l(2,:)*P(4:6,:)],0);
else
  L = inf*ones(4);
end
l = norml(vgg_wedge(reshape(P*L(:,3),[3 K]),reshape(P*L(:,4),[3 K])));

% compute residual function
y = [];
for k = 1:K
  y = [y l(k,:)*s{k}];
end
y = y';

w = [1;1;1] * (1./scale);
w = w(:);

% else, compute also jacobian
if nargout < 2
  return
end
dif = 1e-6;
J = zeros(length(y),length(p));
for i = 1:length(p)
  pdif = p;
  pdif(i) = pdif(i) + dif;
  J(:,i) = (F(pdif,P,s,scale,varargin{:}) - y)/dif;
end
return


% The following function computes lines in the first two images
% from parameters p. Explanation see above.
function l = l12_from_p(p,l0,ldelta,D,d,c)
switch length(p)
 case 4 % unconstrained
  l = [l0(1,:) + p(1:2)'*ldelta{1}
       l0(2,:) + p(3:4)'*ldelta{2}];
 case 2 % going thru X
  l = [l0(1,:) + p(1)*ldelta{1}
       l0(2,:) + p(2)*ldelta{2}];
 case 3 % going thru L
  p(4) = -(p(1:2)'*(D(:,1)*p(3)+d{1})+p(3)*d{2}(1)+c)/(p(1:2)'*D(:,2)+d{2}(2));
  l = [l0(1,:) + p(1:2)'*ldelta{1}
       l0(2,:) + p(3:4)'*ldelta{2}];
end
return


function dp = normsolve(J,Y,w,lambda)
OLDWARN = warning('off');
dp = -( J'*diag(w)*J + lambda*eye(size(J,2)) ) \ ( J'*(Y.*w) );
warning(OLDWARN);
return

function x = normx(x)
if ~isempty(x)
  x = x./(ones(size(x,1),1)*sqrt(sum(x.*x)));
end

function l = norml(l)
% l = norml(l)  Multiplies hyperplane l by scalar so that for each n, norm(l(1:end-1,n))==1. 
l = l./(sqrt(sum(l(:,1:end-1).^2,2))*ones(1,size(l,2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a = levmarq(@RES,PARAMS,@NORMSOLVE,a [,opt])  Non-linear least-squares by Levenberg-Marquardt.
%
% Minimizes f(a)'*W*f(a) over a.
%
% @RES ... residual function f called like [e,w,J] = RES(a,PARAMS{:}), where
%    - a ... double(M,1), parameter vector
%    - e ... double(N,1), residual vector
%    - J ... double(N,M), derivative of e wrt a
%    - w ... double(N,1), weights of e; covariance matrix of e is diag(1/e.^2).
%            Use 1 instead of ones(N,1).
% For efficiency, RES should not compute the jacobian if called with two output parameters only.
%
% @NORMSOLVE ... function solving normal equations, called like da = NORMSOLVE(J,e,W,lambda).
% a ... initial parameter vector
% opt ... options structure, see code

function [a,w] = levmarq(RES,PARAMS,NORMSOLVE,a0,varargin)

% options
[opt,rem_opt] = vgg_argparse( { 'niter_term',     +inf,...
                                'drmsrel_term',    0,...
                                'loglambda_term',  6,...
                                'loglambda_init', -4,...
                                'verbose',         0 },...
                              varargin );
if ~isempty(rem_opt), if ~isempty(fieldnames(rem_opt))
  error(['Unknown option(s) ' fieldnames(rem_opt)]);
end, end


% Initial statistics
if opt.verbose
  [e0,w0] = feval(RES,a0,PARAMS{:});
  ssd0 = sum( (e0.*w0).^2 );
  fprintf( '                         [rms=%14.12g] [maxabs=%14.12g]\n',...
          sqrt(ssd0/length(e0)),...
          max(abs(e0.*w0)) );
end


loglambda = opt.loglambda_init;
niter = 0;
while 1

  % Compute actual residual and jacobian
  [e0,w0,J] = feval(RES,a0,PARAMS{:});  

  % Update a as a := a0 + da, by finding
  % optimal lambda and solving normal equations for da.
  nfail = 1;
  while (loglambda < opt.loglambda_term)
    niter = niter + 1;
  
    a = a0 + feval(NORMSOLVE,J,e0,w0,10^loglambda);
    [e,w] = feval(RES,a,PARAMS{:});

    if sum((e.*w).^2) < sum((e0.*w0).^2) % success
      a0 = a;
      loglambda = loglambda - 1;
      break
    end

    if opt.verbose
      fprintf('%4i.%.2i: [loglambda=%3i] [REJECTED]\n',niter,nfail,loglambda);
    end

    loglambda = loglambda + 1;
    nfail = nfail + 1;
  end

  % Print statistic after successful iteration
  ssd0 = sum( (e0.*w0).^2 );
  ssd = sum( (e.*w).^2 );
  if opt.verbose
    fprintf( '%4i   : [loglambda=%3i] [rms=%14.12g] [maxabs=%14.12g] [drmsrel=%4g%%]\n',...
             niter,...
             round(loglambda),...
             sqrt(ssd/length(e)),...
             max(abs(e.*w)),...
             100*(1-sqrt(ssd/ssd0)) );
  end

  % Termination criteria
  test(1) = loglambda <  opt.loglambda_term;
  test(2) = ssd0-ssd  >= opt.drmsrel_term^2*ssd0;
  test(3) = niter     <  opt.niter_term;
  if any(test==0)
    break
  end
end

if opt.verbose
  onoff = {'YES','no'};
  fprintf( ' Levenberg-Marquardt finished succesfully.\n Reason for termination:\n' );
  fprintf( '   lambda  = %s\n', onoff{test(1)+1} );
  fprintf( '   drmsrel = %s\n', onoff{test(2)+1} );
  fprintf( '   niter   = %s\n', onoff{test(3)+1} );
end

return



function print_statistics(niter,loglambda,e0,w0,e,w,opt)
if opt.verbose
  ssd0 = sum( (e0.*w0).^2 );
  ssd = sum( (e.*w).^2 );
  fprintf( '%4i   : [loglambda=%3i] [rms=%14.12g] [maxabs=%14.12g] [drmsrel=%11.5g]\n',...
           niter,...
           round(loglambda),...
           sqrt(ssd/length(e)),...
           max(abs(e.*w)),...
           sqrt(1-ssd/ssd0) );
end
return

