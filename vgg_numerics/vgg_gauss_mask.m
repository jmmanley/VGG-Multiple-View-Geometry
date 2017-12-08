function Y = vgg_gauss_mask(a,der,X,flag)

%VGG_GAUSS_MASK Univariate Gaussian function or its derivatives.
%   Y = vgg_gausmask(sigma,der,X) yields der-th derivative of Gaussian.
%   For der==0, plain Gaussian is yielded.
%   If X is omitted, X=[-3*sigma:3*sigma] is used, giving reasonably small tail.
%   If der is omitted too, der=0 is used.
%
% Try  plot(vgg_gauss_mask(10,6,-100:100))

% Tom Werner, Oct 2001


%%
 % For computing derivatives, the following recurrent formulas are used:
 %   gausmask(sigma,0,X) = exp(a*X.^2)/(sigma*sqrt(2*pi)),
 %   gausmask(sigma,n,X) = 2*a*( (der-1)*gausmask(sigma,der-2,X) + X.*gausmask(sigma,der-1,X) ),
 % where
 %   a = -0.5/(a*a).
 %%

if nargin == 4
  if der > 0
    Y = 2*a*( (der-1)*vgg_gauss_mask(a,der-2,X,0) + X.*vgg_gauss_mask(a,der-1,X,0) );
  elseif der == 0
    Y = ones(size(X));
  else
    Y = zeros(size(X));
  end
else
  if nargin < 3
    X = ceil(3*a);
    X = -X:X;
  end
  if nargin < 2
    der = 0;
  end
  aa = -0.5/(a*a);
  Y = exp(aa*X.^2)/(a*sqrt(2*pi)) .* vgg_gauss_mask(aa,der,X,0);
end

return