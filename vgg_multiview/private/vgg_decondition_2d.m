function pd = vgg_decondition_2d(p,C);
% function pd = vgg_decondition_2d(p,C);
%
% Decondition a set of 2D homogeneous or nonhomogeneous points using (original) conditioner C

[r,c] = size(p);
if r == 2
  pd = vgg_get_nonhomg(inv(C) * vgg_get_homg(p));
elseif r == 3
  pd = inv(C) * p;
else
  error ('rows != 2 or 3');
end
