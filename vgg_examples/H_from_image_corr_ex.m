% Example on computing H from a set correspondences for planar points (on the
% floor) using vgg_H_from_x_lin and displaying result

% Note points are read from stored example in vgg_example_scene

% Read in all detected interest points and images
[view, Xi, X, Li, L] = vgg_example_scene(2);

n=1; k=2; 

% select points for which there are correspondences
i=all(Xi([n k],:)>0); i=find(i);

% compute H from points on the floor
H=vgg_H_from_x_lin(view(n).x(:,Xi(n,i([406 399 405 133 115 109]))),...
             view(k).x(:,Xi(k,i([406 399 405 133 115 109]))));

% Display
vgg_gui_H(view(n).I, view(k).I, H)
disp('Computed homography of the floor. Move the mouse to verify')
