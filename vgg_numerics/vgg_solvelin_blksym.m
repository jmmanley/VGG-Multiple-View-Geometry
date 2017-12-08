%VGG_SOLVELIN_BLKSYM  Solves M*x==y where M is (typically huge sparse) symmetric 4-block matrix.
%   It solves the system much more efficiently than a general (sparse) linear system solver.
%   Typical usage to solve normal equations in Levenberg-Marquardt in bundle adjustment.
% 
%   X = VGG_SOLVELIN_BLKSYM(A,B,C,p,q [,'nocheck']) solves M*x==Y where 
%   M=[A B; B' C] and y=[p;q].
%
%   X = VGG_SOLVELIN_BLKSYM(M,y,sideA [,'nocheck']) does the same but takes M and splits it,
%   M=[A B; B' C] and size(A)==[sideA sideA].
%
%   The function checks whether C is really near to diagonal; if not, a warning is printed.
%   Parmeter 'nocheck' switches off this test.

% (c) {awf,werner}@robots.ox.ac.uk, March 2002


function x = vgg_solvelin_blksym(varargin)

switch nargin
 case 6
  [A,B,C,p,q,check] = deal(varargin{:});
 case 5
  [A,B,C,p,q] = deal(varargin{:});
  check = '';
 case 4
  [M,y,sideA,check] = deal(varargin{:});
 case 3
  [M,y,sideA] = deal(varargin{:});
  check = '';
 otherwise
  error('Bad number of parameters');
end

if nargin<5
  Rtop = 1:sideA;
  Rbot = sideA+1:size(M,1);
  C = M(Rbot, Rbot); 
  B = M(Rtop, Rbot);
  A = M(Rtop, Rtop);
  p = y(Rtop);
  q = y(Rbot);
end

if ~strcmp(check,'nocheck') & nnz(C(1,:))/size(C,1) > .1
  warning('It is likely that M was splitted incorrectly. Check diagonality of C and/or value of sideA.');
end

%tic

% Surprisingly, branch 1 is 2x slower than branch 2. We don't know why.
switch 2
 case 1
  invC_Btq = C \ [B' q];
 case 2
  invC_Btq = inv(C) * [B' q];
end

invC_Bt = invC_Btq(:,1:end-1);
invC_q = invC_Btq(:,end);
      
u = (A - B * invC_Bt) \ (p - B * invC_q);
v = invC_q - invC_Bt*u;

x = [u;v];

%toc

return