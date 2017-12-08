function [C,invC] = vgg_conditioner_from_image(c,r)
%
% function [C,invC] = vgg_conditioner_from_image(image_width, image_height)
%
%   Makes a similarity metric for conditioning image points.
%
%   Also can be called as vgg_conditioner_from_image([image_width image_height])
%
%   invC is inv(C), obtained more efficiently inside the function.

if nargin<2
  r = c(2);
  c = c(1);
end

f = (c+r)/2;
C = [1/f 0 -c/(2*f) ;
     0 1/f -r/(2*f) ;
     0 0 1];

if nargout > 1
  invC = [f 0 c/2 ;
          0 f r/2 ;
          0 0 1];
end