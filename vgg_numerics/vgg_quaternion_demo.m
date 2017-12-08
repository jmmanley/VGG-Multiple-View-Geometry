

% Make some synthetic quaternions
syms   q0 q1 q2 q3   p0 p1 p2 p3   real

q = [q0 q1 q2 q3]';
p = [p0 p1 p2 p3]';

vgg_assert_equal('vgg_quat_mul(p,q)', 'vgg_quat_matrix(p) * q', 1e-8,1);

% Generate a random unit quaternion.
% This represents a 3D rotation.
q = randn(4,1);
q = q / norm(q);

% Now convert it to a rotation matrix,
% and convert back to ensure these routines
% are self-inverse
R = vgg_quat_rotation_matrix(q);
qback = vgg_quat_from_rotation_matrix(R);

vgg_assert_equal(q, qback, -1e-8,1);

