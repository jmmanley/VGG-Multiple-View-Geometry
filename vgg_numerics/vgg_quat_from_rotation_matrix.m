function q = vgg_quat_from_rotation_matrix( R )
% vgg_quat_from_rotation_matrix Generates quaternion from rotation matrix 
%            q = vgg_quat_from_rotation_matrix(R)

q = [	(1 + R(1,1) + R(2,2) + R(3,3))
	(1 + R(1,1) - R(2,2) - R(3,3))
	(1 - R(1,1) + R(2,2) - R(3,3))
	(1 - R(1,1) - R(2,2) + R(3,3)) ];

if ~issym(q)
  % Pivot to avoid division by small numbers
  [b I] = max(abs(q));
else
  % For symbolic quats, just make sure we're nonzero
  for k=1:4
    if q(k) ~= 0
      I = k;
      break
    end
  end
end

q(I) = sqrt(q(I)) / 2 ;

if I == 1 
	q(2) = (R(3,2) - R(2,3)) / (4*q(I));
	q(3) = (R(1,3) - R(3,1)) / (4*q(I));
	q(4) = (R(2,1) - R(1,2)) / (4*q(I));
elseif I==2
	q(1) = (R(3,2) - R(2,3)) / (4*q(I));
	q(3) = (R(2,1) + R(1,2)) / (4*q(I));
	q(4) = (R(1,3) + R(3,1)) / (4*q(I));
elseif I==3
	q(1) = (R(1,3) - R(3,1)) / (4*q(I));
	q(2) = (R(2,1) + R(1,2)) / (4*q(I));
	q(4) = (R(3,2) + R(2,3)) / (4*q(I));
elseif I==4
	q(1) = (R(2,1) - R(1,2)) / (4*q(I));
	q(2) = (R(1,3) + R(3,1)) / (4*q(I));
	q(3) = (R(3,2) + R(2,3)) / (4*q(I));
end
