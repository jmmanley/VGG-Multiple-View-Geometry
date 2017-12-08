function PM = vgg_quat_matrix(p)
% VGG_QUAT_MATRIX  Quaternion multiplication matrix
%               PM = vgg_quat_matrix(q) returns PM such that
%               PM * q = vgg_quat_mul(p,q)

% awf@robots.ox.ac.uk, 30/01/05

PM = [
  [  p(1), -p(2), -p(3), -p(4)]
  [  p(2),  p(1), -p(4),  p(3)]
  [  p(3),  p(4),  p(1), -p(2)]
  [  p(4), -p(3),  p(2),  p(1)]
  ];
