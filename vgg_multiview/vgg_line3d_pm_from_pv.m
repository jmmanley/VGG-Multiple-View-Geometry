% G = vgg_line3d_pm_from_pv(L) Conversion of Pluecker vector to Puecker matrix 3d line representation.
%
% L ... double(1,6), Pluecker vector of the 3d line
% G ... double(4,4), Puecker matrix

% T.Werner

function G = vgg_line3d_pm_from_pv(L)

G = [vgg_contreps(L(4:6)) L(1:3)'; -L(1:3) 0];

return