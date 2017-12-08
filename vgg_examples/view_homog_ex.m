% Example of using vgg_gui_H

% read in basement images

im1 = imread('bt.000.png');
im2 = imread('bt.002.png');

% read in homography for ground plane
H = load('bt.00.02.H');

% view homography in GUI, move mouse with button down
% in either window to see transferred point
vgg_gui_H(im1, im2, H)
