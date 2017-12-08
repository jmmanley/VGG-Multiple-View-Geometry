function e = vgg_rms_rrror(M)
% e = vgg_rms_rrror(M)
%
% Get RMS diff from zero of matrix or vector M

e = sqrt(sum(sum(M.*M)) / prod(size(M)));
