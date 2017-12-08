function M = vgg_rotmat_from_exp(r)
% VGG_ROTMAT_FROM_EXP  Convert from exponential to matrix rotation parameterization.
%            R = vgg_rotmat_from_exp([r1, r2, r3]) generates the rotation
%            matrix with axis along r, angle = norm(r).
%            This is equivalent to expm(cross_matrix(r)) but more stable.

% Andrew Fitzgibbon <awf@robots.ox.ac.uk>

H = [0, -r(3), r(2);  r(3), 0, -r(1); -r(2), r(1), 0];

if 1
  angle = norm(r);
  if (angle < eps)
    M=eye(3,3);
  else
    ef = sin(angle)/angle;
    gee = (1.0 - cos(angle))/ (angle*angle);
    M = (H*H)*gee + H*ef + eye(3,3);
  end
else 
  M = expm(H);
end
