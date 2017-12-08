% Demonstration of feature matching via simple correlation, and then using
% RANSAC to estimate the homography matrix and at the same time identify
% (mostly) inlying matches

% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% February 2004

% Adapted to use vgg functions by Peter Kovesi and Andrew Zisserman
% 2012 Updated by Alexander Khanin

function H = testhomog_vgg
    
    % prepare workspace
    close all
    clear all
%    clc
    
    % block of settings
    thresh = 100;   % Harris corner threshold
    dmax = 100;     % 
    w = 11;         % window size for correlation matching
    
    % try to read a pair of images
    try
        im1 = imread('keble.000.png');
        im2 = imread('keble.003.png');
    catch ME
        warning(ME.message);
        return;
    end

    % imput images must be grayscale
    if ndims(im1) == 3, im1 = rgb2gray(im1); end
    if ndims(im2) == 3, im2 = rgb2gray(im2); end
    
    % Find Harris corners in image1 and image2
    [cim1, r1, c1] = harris(im1, 1, thresh, 3);
    [cim2, r2, c2] = harris(im2, 1, thresh, 3);

    % create figure and fit it size to screen size
    scrsz = get(0,'ScreenSize');
    figure('Name', 'Example', 'Position',[1 1 scrsz(3) scrsz(4)])
    
    % display detected corners
    subplot(3,4,1),imshow(im1),title('Corners on 1st image');
    hold on, plot(c1,r1,'r+');
    subplot(3,4,5),imshow(im2),title('Corners on 2nd image');
    hold on, plot(c2,r2,'r+');

    % find correlation
    [m1,m2] = matchbycorrelation(im1, [r1';c1'], im2, [r2';c2'], w, dmax);

    % Display putative matches
    subplot(3,4,9),imshow(im1),title('Putative matches');
    for n = 1:length(m1);
        line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)])
    end

    % Assemble homogeneous feature coordinates for fitting of the
    % homography matrix, note that [x,y] corresponds to [col, row]
    x1 = [m1(2,:); m1(1,:); ones(1,length(m1))];
    x2 = [m2(2,:); m2(1,:); ones(1,length(m1))];    
    
    t = .001;  % Distance threshold for deciding outliers
    [H, inliers] = ransacfithomography_vgg(x1, x2, t);

    fprintf('Number of inliers was %d (%d%%) \n', ...
	    length(inliers),round(100*length(inliers)/length(m1)))
    fprintf('Number of putative matches was %d \n', length(m1))        
    
    % Display both images overlayed with inlying matched feature points
    subplot(3,4, [2 3 4 6 7 8 10 11 12]);
    imagesc(double(im1)+double(im2)), title('Inlying matches'), hold on 
    plot(m1(2,inliers),m1(1,inliers),'r+');
    plot(m2(2,inliers),m2(1,inliers),'g+');    

    % Step through each matched pair of points and display the
    % line linking the points on the overlayed images.
    for n = inliers
        line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 0 1])
    end
   
    return
end
    
