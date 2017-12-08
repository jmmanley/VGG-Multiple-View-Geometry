% Example on computing F from two camera matrices using vgg_F_from_P
% and displaying result

% Note P's are read from stored example in vgg_example_scene

% Read in P's and images
[view] = vgg_example_scene(2);

% Compute F from P's
F = vgg_F_from_P(view(1).P, view(2).P);

% Display
vgg_gui_F(view(1).I, view(2).I, F');
disp('Computed epipolar geometry. Move the mouse to verify')
