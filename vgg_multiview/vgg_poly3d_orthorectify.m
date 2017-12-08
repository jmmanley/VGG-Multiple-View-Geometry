%vgg_poly3d_orthorectify  Frontoparallel rectification of image of 3D polygon.
%   [H,imsize] = vgg_poly3d_orthorectify(Q,u) finds homography H that
%   removes perspective distortion of image Q*hom(u) of a 3D polygon. Parameters:
%     u ... double(2,N), inhomog. coordinates of polygon vertices measured in an orthonormal
%       coordinate frame in the 3D polygon plane.
%     Q ... double(3,3), homography that maps u to image plane.
%       Image of the polygon is thus x = Q*hom(u) (homog. coords.).
%     H ... double(3,3), rectifying homography
%     imsize ... double(2,1), size of the target image
%   The rectification is done as follows:
%     Irect = vgg_warp_H(I{k},inv(H),'linear',[1 imsize(2) 1 imsize(1)]).
%
%   [H,imsize] = vgg_poly3d_orthorectify(Q,u,smax) allows to specify maximum ratio of
%   area element in output and input image (if omitted, smax=5). This is to prevent
%   the output image from becoming huge if the perspective distortion is large.
%
%   If Q is K-cell of double(3,3), rectification is done simultaneously for K images
%   of the 3D polygon, so that pixels with the same coordinates correspond in the
%   rectified images. Then H is K-cell of double(3,3) and imsize is double(2,K).
%
%   [H,imsize] = vgg_poly3d_orthorectify(Q,u,smax,in_imsize) specifies sizes of input
%   images. This is needed if some vertices project outside some input images.

function [H,imsize] = vgg_poly3d_frontorectify(H,u,smax,imsize)

if ~iscell(H)
  H = {H};
end

if nargin<3
  smax = 5;
end

% find scales of H in all u
for k = 1:length(H)
  if nargin>3 % consider only vertices inside input image
    i = boxclipu( [[1;1] imsize(:,k)], vgg_get_nonhomg(H{k}*vgg_get_homg(u)) );
  else
    i = logical(ones(1,size(u,2)));
  end
  s(k,~i) = nan;
  if nnz(i)==0, continue, end
  s(k,i) = sqrt(abs( dethomog(H{k},u(:,i)) ));
end

% determine the scale;
% forbid scaling polygons more than smax-times
s = s(:);
s(isnan(s)) = [];
if isempty(s) % polygon is not in the image
  imsize = [0;0];
  return
end
if nargin < 5
  smax = inf;
end
s = min(s) * min(smax,max(s)/min(s));

% scale H and u consistently by s
u = u*s;
Hs = diag([1/s 1/s 1]);
for k = 1:length(H)
  H{k} = H{k}*Hs;
end

% translate the polygon to the origin
b = bboxu(u);
imsize = ceil( b(:,2) - b(:,1) );
Ht = eye(3);
Ht(1:2,3) = b(:,1) + 1;

for k = 1:length(H)
  H{k} = H{k}*Ht;
end

if length(H)==1
  H = H{1};
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% D = dethomog(H,q)  Determinant of Jacobian of homography transformation r = nhom(H*hom(q)), D = det(dr/dq). Works vector-wise for size(q)=[2 ?].
%
function D = dethomog(H,q)
% Derivative of homography r = h(q) = nhom(H*[q;1]) is:
%   dh = [eye(2) -r]/y(3)*H(:,1:2), where y = H*[q;1].
q(3,:) = 1;
r = vgg_get_nonhomg(H*q);
D = ((H(1,1)-r(1,:)*H(3,1)).*(H(2,2)-r(2,:)*H(3,2)) - (H(1,2)-r(1,:)*H(3,2)).*(H(2,1)-r(2,:)*H(3,1))) ./ (H(3,:)*q).^2;
return


function i = boxclipu(b,u)
% i = boxclipu(b,u)  Returns logical indices of points inside box b.
i = all( u >= b(:,row1(u)) & u <= b(:,2*row1(u)), 1 );
return


function ub = bboxu(u)
% ub = bboxu(u)  Bounding box of point cloud u.
ub = [min(u')' max(u')'];
return