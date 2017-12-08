% Example of using vgg_gui_F

% read in chapel images

im1 = imread('chapel00.png');
im2 = imread('chapel01.png');

% read in fundamental matrix
F = load('chapel.00.01.F');

% view epipolar geometry in GUI, move mouse with button down
% in either window to see transferred point
% NB function uses F transpose
vgg_gui_F(im1, im2, F')
