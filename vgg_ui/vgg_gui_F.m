function fig=vgg_gui_F(i1,i2,F)
%
%	fig=vgg_gui_F(i1,i2,F)
%
%
% Visualizes the fundamental matrix of two views
%
%IN:
%	i1 - Matlab image
%	i2 - Matlab image
%	F - Fundamental matrix (p1'*F*p2=0). Assumes that image coordiantes
%		are 1..width where pixel centers are at integer locations. 
%
%OUT:
%	fig - handle to the figure

% $Id: vgg_gui_F.m,v 1.8 2002/05/20 21:35:42 wexler Exp $
% Yoni, Tue Mar 27 19:31:16 2001


if nargin==3
   action='start';
else
   action=i1;
   ud = get(gcf, 'UserData');
end

if strcmp(action,'start'),
   if nargin ~= 3
      error('Must give 3 arguments... read the docs.\n');
   end

   h0 = figure('Color',[0.8 0.8 0.8], ...
	       'NumberTitle','off', ...
	       'Name','Play With Fundamental Matrix', ...
	       'ButtonDownFcn', 'disp(''Click on images'')',...
	       'WindowButtonUpFcn', 'vgg_gui_F(''none'');',...
	       'WindowButtonMotionFcn', 'vgg_gui_F(''move'')', ...
	       'Pointer', 'crosshair', ...
	       'DoubleBuffer', 'on',...
	       'Units','normalized');
   [pointerShape, pointerHotSpot] = CreatePointer;
   set(h0, 'Pointer', 'custom', ...
	   'PointerShapeCData', pointerShape, ...
	   'PointerShapeHotSpot', pointerHotSpot);
   
   m0=uimenu('Label', '&Color');
   uimenu(m0, 'Label', 'blac&K', 'ForegroundColor', [0 0 0], ...
	  'Accelerator', 'k', 'Callback', 'vgg_gui_F(''ck'');');
   uimenu(m0, 'Label', '&Red', 'ForegroundColor', [1 0 0], ...
	  'Accelerator', 'r', 'Callback', 'vgg_gui_F(''cr'');');
   uimenu(m0, 'Label', '&Green', 'ForegroundColor', [0 1 0], ...
	  'Accelerator', 'g', 'Callback', 'vgg_gui_F(''cg'');');
   uimenu(m0, 'Label', '&Blue', 'ForegroundColor', [0 0 1], ...
	  'Accelerator', 'b', 'Callback', 'vgg_gui_F(''cb'');');
   m1=uimenu('Label', '&Size');
   uimenu(m1, 'Label', '&Increase', 'Callback', 'vgg_gui_F(''s+'');', 'Accelerator', '+');
   uimenu(m1, 'Label', '&Decrease', 'Callback', 'vgg_gui_F(''s-'');', 'Accelerator', '-');
   uimenu(m1, 'Label', '&1', 'Callback', 'vgg_gui_F(''s1'');', 'Accelerator', '1');
   uimenu(m1, 'Label', '&2', 'Callback', 'vgg_gui_F(''s2'');', 'Accelerator', '2');
   uimenu(m1, 'Label', '&3', 'Callback', 'vgg_gui_F(''s3'');', 'Accelerator', '3');
   uimenu(m1, 'Label', '&4', 'Callback', 'vgg_gui_F(''s4'');', 'Accelerator', '4');
   uimenu(m1, 'Label', '&5', 'Callback', 'vgg_gui_F(''s5'');', 'Accelerator', '5');
   uimenu(m1, 'Label', '&6', 'Callback', 'vgg_gui_F(''s6'');', 'Accelerator', '6');
   uimenu(m1, 'Label', '&7', 'Callback', 'vgg_gui_F(''s7'');', 'Accelerator', '7');
   uimenu(m1, 'Label', '&8', 'Callback', 'vgg_gui_F(''s8'');', 'Accelerator', '8');
   uimenu(m1, 'Label', '&9', 'Callback', 'vgg_gui_F(''s9'');', 'Accelerator', '9');
   
   ah1 = axes('Parent', h0, ...
	      'Position',[0 0 .5 1]);
   h1=imshow(i1); hold on; title('Image 1');
   set(h1, 'ButtonDownFcn','vgg_gui_F(''b1'');');

   ah2 = axes('Parent',h0, ...
	      'Position',[.5 0 .5 1], ...
	      'Tag','Axes2');
   h2=imshow(i2); hold on; title('Image 2');
   set(h2, 'ButtonDownFcn','vgg_gui_F(''b2'');');

   point=plot(-1000, -1000,'EraseMode','xor');
   l=plot([-1000, -1001], [-1000 -1000], 'r-','EraseMode','xor');

   s1=size(i1); s2=size(i2);
   t(:,:,1)=F';  t(:,:,2)=F;  F=t;

   ud=struct('h0', h0, 'h',[h1 h2], 'ah', [ah1, ah2], ...
	     'sizes', [s1(1:2); s2(1:2)], ...
	     'current', -1, 'color', 'k', 'size', 1, ...
	     'p', point, 'F', F, 'l', l );

   set(h0,'UserData',ud);

   if nargout > 0, fig = h0; end
elseif strcmp(action, 'move')
   if ud.current<0 return; end;
   pt=get(ud.ah(ud.current),'CurrentPoint');
   pts=get_line_points(ud.F(:,:,ud.current)*pt(1,:)', ...
		       ud.sizes(ud.current,:));

   set(ud.l, 'XData', pts(1,:), 'YData', pts(2,:))
   set(ud.p, 'XData', pt(1,1), 'YData', pt(1,2))

elseif action(1)=='b'
   if action(2)=='1' ud.current=1;
   elseif action(2)=='2' ud.current=2;
   else return;
   end
   pt=get(ud.ah(ud.current),'CurrentPoint');

   pts=get_line_points(ud.F(:,:,ud.current)*pt(1,:)', ud.sizes(1,:));
   delete(ud.p);
   delete(ud.l);
   axes(ud.ah(ud.current));
   ud.p=plot(pt(1,1), pt(1,2), [ud.color '+'], ...
	     'MarkerSize', 8+2*ud.size, 'LineWidth', ud.size,...
	     'EraseMode','xor');
   axes(ud.ah(3-ud.current));
   ud.l=plot(pts(1,:), pts(2,:), [ud.color '-'], ...
	     'LineWidth', ud.size, 'EraseMode','xor');
elseif action(1)=='c'
   ud.color=action(2);
   %get(ud.l)
   set(ud.l, 'Color', ud.color);
   set(ud.p, 'Color', ud.color);
elseif action(1)=='s'
   if action(2)=='+' ud.size=ud.size+1;
   elseif action(2)=='-' ud.size=max(1, ud.size-1);
   else
      ud.size = str2num(action(2));
   end
   set(ud.p, 'LineWidth', ud.size, 'MarkerSize', 8+2*ud.size);
   set(ud.l, 'LineWidth', ud.size, 'MarkerSize', 8+2*ud.size);
elseif strcmp(action, 'none')
   ud.current = -1;
else
   error(['Unknown command: ' action]);
end

set(ud.h0, 'UserData',ud);



function pts=get_line_points(l,sz)
a=l(1); b=l(2);c=l(3);
h=sz(1); w=sz(2);

% This might cause 'divide by zero' warning:
ys=c/-b ;
yf=-(a*w+c)/b;
xs=c/-a;
xf=-(b*h+c)/a;

m1 = [[xs;1] [xf;h] [1;ys] [w;yf]];
w2 = [(xs<=w & xs>=1) (xf<=w & xf>=1) (ys<=h & ys>=1) (yf<=h & yf>=1)];
v = w2>0;
pts = [m1(:,v)];



% Taken from pixval.m:
function [pointerShape, pointerHotSpot] = CreatePointer

pointerHotSpot = [8 8];
pointerShape = [ ...
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
    1   1   1   1   1    1   2  NaN  2   1   1   1   1   1   1   1
    2   2   2   2   2    2  NaN NaN NaN  2   2   2   2   2   2   2
   NaN NaN NaN NaN NaN  NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    2   2   2   2   2    2  NaN NaN NaN  2   2   2   2   2   2   2
    1   1   1   1   1    1   2  NaN  2   1   1   1   1   1   1   1
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN   1   2  NaN  2   1 NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN NaN NaN NaN  NaN NaN NaN NaN NaN NaN NaN NaN NaN];
