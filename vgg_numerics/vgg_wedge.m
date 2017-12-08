% vgg_wedge  Wedge product of N-1 N-vectors (generalization of cross product).
%
% Y = vgg_wedge(X)  Wedge product of columns/rows of X.
%   Y ... double (1,N).
%   X ... double (N,N-1).
% It is Y = X(:,1) \wedge X(:,2) \wedge ... \wedge X(:,N-1). For N=3,
% wedge product is the same as cross (vector) product. E.g., for N=4,
% wedge product of three 4-vectors is in fact computation of a plane in
% 3-D projective space from 3 points in the plane.
% The sign is chosen so that for any square matrix X it is: det(X) == wedge(X(:,1:end-1))*X(:,end)
% Works also dually for
%   Y ... double (N,1).
%   X ... double (N-1,N).
%
% Y = vgg_wedge(X_1,X_2,...,X_{N-1})  Wedge product for each (N-1)-tuple of corresponding columns of X_n.
%   X_n ... double (N,K)
%   Y ... double (K,N)
% Equivalent to
%   for k = 1:K, Y(k,:) = wedge([X_1(:,k) ... X_{N-1}(:,k)]); end
% E.g.: wedge(X1,X2) is the same as cross(X1,X2)' but faster.
% Dual form is not available.

function Y = vgg_wedge(varargin)

if nargin == 1

  X = varargin{1};
  
  [N,Nm1] = size(X);
  if Nm1>N
    Y = vgg_wedge(X')';
    return
  end

  switch N
   case 3 % make it faster for special case N==3
    Y = [X(2,1).*X(3,2)-X(3,1).*X(2,2),...
         X(3,1).*X(1,2)-X(1,1).*X(3,2),...
         X(1,1).*X(2,2)-X(2,1).*X(1,2)];
   otherwise
    for n = 1:N
      Y(n) = (-1)^(n+N)*det(X([1:n-1 n+1:N],:));
    end
  end

else
  
  N = nargin + 1;
  switch N
   case 3 % make it faster for special case N==3
    X1 = varargin{1};
    X2 = varargin{2};
    Y = [(X1(2,:).*X2(3,:)-X1(3,:).*X2(2,:))',...
         (X1(3,:).*X2(1,:)-X1(1,:).*X2(3,:))',...
         (X1(1,:).*X2(2,:)-X1(2,:).*X2(1,:))'];
   otherwise
    for k = 1:size(varargin{1},2)
      for n = 1:N-1
	X(:,n) = varargin{n}(:,k);
      end
      Y(k,:) = vgg_wedge(X);
    end
  end

end

  
return