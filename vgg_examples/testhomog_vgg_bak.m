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

function H = testhomog_vgg

    close all    
    
    thresh = 500;   % Harris corner threshold
    nonmaxrad = 3;  % Non-maximal suppression radius
    dmax = 100;
    w = 11;    % Window size for correlation matching
    
    im1 = rgb2gray(imread('keble.000.png'));
    im2 = rgb2gray(imread('keble.003.png'));

    % Find Harris corners in image1 and image2
    [cim1, r1, c1] = harris(im1, 1, thresh, 3);
    show(im1,1), hold on, plot(c1,r1,'r+');

    [cim2, r2, c2] = harris(im2, 1, thresh, 3);
    show(im2,2), hold on, plot(c2,r2,'r+');

    drawnow

tic
    [m1,m2] = matchbycorrelation(im1, [r1';c1'], im2, [r2';c2'], w, dmax);
toc
    % Display putative matches
    show(im1,3), set(3,'name','Putative matches'), hold on    
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
    show(double(im1)+double(im2),4), set(4,'name','Inlying matches'), hold on 
    plot(m1(2,inliers),m1(1,inliers),'r+');
    plot(m2(2,inliers),m2(1,inliers),'g+');    

    % Step through each matched pair of points and display the
    % line linking the points on the overlayed images.

    for n = inliers
	line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 0 1])
    end
   
    return
    
    
