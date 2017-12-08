function qout = vgg_quat_mul(r,q)
% VGG_QUAT_MUL  Quaternion multiplication
%               pq = vgg_quat_mul(p,q)

% awf@robots.ox.ac.uk, 30/01/05

qout = ...
     [(r(1)*q(1) - r(2)*q(2) - r(3)*q(3) - r(4)*q(4))  
      (r(1)*q(2) + r(2)*q(1) + r(3)*q(4) - r(4)*q(3))  
      (q(3)*r(1) - q(4)*r(2) + q(1)*r(3) + q(2)*r(4))  
      (q(4)*r(1) + q(3)*r(2) - q(2)*r(3) + q(1)*r(4)) ];
