%VGG_SINGF_FROM_FF  Linearly combines two 3x3 matrices to a singular one.
%
%   a = vgg_singF_from_FF(F)  computes scalar(s) a such that given two 3x3 matrices F{1} and F{2},
%   it is det( a*F{1} + (1-a)*F{2} ) == 0.

function a = vgg_singF_from_FF(F)

% precompute determinants made from columns of F{1}, F{2}
for i1 = 1:2
  for i2 = 1:2
    for i3 = 1:2
      D(i1,i2,i3) = det([F{i1}(:,1) F{i2}(:,2) F{i3}(:,3)]);
    end
  end
end

% Solve The cubic equation for a
a = roots([-D(2,1,1)+D(1,2,2)+D(1,1,1)+D(2,2,1)+D(2,1,2)-D(1,2,1)-D(1,1,2)-D(2,2,2)
            D(1,1,2)-2*D(1,2,2)-2*D(2,1,2)+D(2,1,1)-2*D(2,2,1)+D(1,2,1)+3*D(2,2,2)
            D(2,2,1)+D(1,2,2)+D(2,1,2)-3*D(2,2,2)
            D(2,2,2)]);
a = a(abs(imag(a))<10*eps);

return