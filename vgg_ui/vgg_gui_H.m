function fig=vgg_gui_H(i1,i2,H)
%
%	fig=vgg_gui_H(i1,i2,H)
%
%
% Visualizes an homography matrix of two views
%
%IN:
%	i1 - Matlab image
%	i2 - Matlab image
%	H - 3x3 Homography matrix. Assumes that image coordiantes are 1..width
%		where pixel centers are at integer locations.
%
%OUT:
%	fig - handle to the figure

% $Id: vgg_gui_H.m,v 1.3 2001/10/30 14:42:41 wexler Exp $
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
	       'Name','Play With Homography Matrix', ...
	       'ButtonDownFcn', 'disp(''Click on images'')',...
	       'WindowButtonUpFcn', 'vgg_gui_H(''none'');',...
	       'WindowButtonMotionFcn', 'vgg_gui_H(''move'')', ...
	       'Pointer', 'crosshair', ...
	       'DoubleBuffer', 'on',...
	       'Units','normalized');
   m0=uimenu('Label', '&Color');
   uimenu(m0, 'Label', '&Red', 'ForegroundColor', [1 0 0], ...
	  'Accelerator', 'r', 'Callback', 'vgg_gui_H(''cr'');');
   uimenu(m0, 'Label', '&Green', 'ForegroundColor', [0 1 0], ...
	  'Accelerator', 'g', 'Callback', 'vgg_gui_H(''cg'');');
   uimenu(m0, 'Label', '&Blue', 'ForegroundColor', [0 0 1], ...
	  'Accelerator', 'b', 'Callback', 'vgg_gui_H(''cb'');');
   m1=uimenu('Label', '&Size');
   uimenu(m1, 'Label', '&Increase', 'Callback', 'vgg_gui_H(''s+'');', 'Accelerator', '+');
   uimenu(m1, 'Label', '&Decrease', 'Callback', 'vgg_gui_H(''s-'');', 'Accelerator', '-');
   uimenu(m1, 'Label', '&1', 'Callback', 'vgg_gui_H(''s1'');', 'Accelerator', '1');
   uimenu(m1, 'Label', '&2', 'Callback', 'vgg_gui_H(''s2'');', 'Accelerator', '2');
   uimenu(m1, 'Label', '&3', 'Callback', 'vgg_gui_H(''s3'');', 'Accelerator', '3');
   uimenu(m1, 'Label', '&4', 'Callback', 'vgg_gui_H(''s4'');', 'Accelerator', '4');
   uimenu(m1, 'Label', '&5', 'Callback', 'vgg_gui_H(''s5'');', 'Accelerator', '5');
   uimenu(m1, 'Label', '&6', 'Callback', 'vgg_gui_H(''s6'');', 'Accelerator', '6');
   uimenu(m1, 'Label', '&7', 'Callback', 'vgg_gui_H(''s7'');', 'Accelerator', '7');
   uimenu(m1, 'Label', '&8', 'Callback', 'vgg_gui_H(''s8'');', 'Accelerator', '8');
   uimenu(m1, 'Label', '&9', 'Callback', 'vgg_gui_H(''s9'');', 'Accelerator', '9');
   
   ah1 = axes('Parent', h0, ...
	      'Position',[0 0 .5 1]);
   h1=imshow(i1); hold on; title('Image 1');
   set(h1, 'ButtonDownFcn','vgg_gui_H(''b1'');');

   ah2 = axes('Parent',h0, ...
	      'Position',[.5 0 .5 1], ...
	      'Tag','Axes2');
   h2=imshow(i2); hold on; title('Image 2');
   set(h2, 'ButtonDownFcn','vgg_gui_H(''b2'');');

   point=plot(-1000, -1000,'EraseMode','xor');
   point2=plot(-1000, -1000,'EraseMode','xor');

   s1=size(i1); s2=size(i2);
   t(:,:,1)=H;  t(:,:,2)=inv(H);  H=t;

   ud=struct('h0', h0, 'h',[h1 h2], 'ah', [ah1, ah2], ...
	     'sizes', [s1(1:2); s2(1:2)], ...
	     'current', -1, 'color', 'r', 'size', 1, ...
	     'p', point, 'H', H, 'p2', point2 );

   set(h0,'UserData',ud);

   if nargout > 0, fig = h0; end

elseif strcmp(action, 'move')
   if ud.current<0 return; end;
   pt=get(ud.ah(ud.current),'CurrentPoint');
   pt2=ud.H(:,:,ud.current)*pt(1,:)';
   pt2=pt2/pt2(3);

   set(ud.p2, 'XData', pt2(1,1), 'YData', pt2(2,1))
   set(ud.p, 'XData', pt(1,1), 'YData', pt(1,2))

elseif action(1)=='b'
   if action(2)=='1' ud.current=1;
   elseif action(2)=='2' ud.current=2;
   else return;
   end
   pt=get(ud.ah(ud.current),'CurrentPoint');

   p2=ud.H(:,:,ud.current)*pt(1,:)';
   p2=p2/p2(3);

   delete(ud.p);
   delete(ud.p2);
   axes(ud.ah(ud.current));
   ud.p=plot(pt(1,1), pt(1,2), [ud.color '+'], ...
	     'MarkerSize', 8+2*ud.size, 'LineWidth', ud.size,...
	     'EraseMode','xor');
   axes(ud.ah(3-ud.current));
   ud.p2=plot(p2(1,1), p2(2,1), [ud.color '+'], ...
	     'MarkerSize', 8+2*ud.size, 'LineWidth', ud.size,...
	     'EraseMode','xor');

elseif action(1)=='c'
   ud.color=action(2);
   %get(ud.l)
   set(ud.p2, 'Color', ud.color);
   set(ud.p, 'Color', ud.color);

elseif action(1)=='s'
   if action(2)=='+' ud.size=ud.size+1;
   elseif action(2)=='-' ud.size=max(1, ud.size-1);
   else
      ud.size = str2num(action(2));
   end
   set(ud.p, 'LineWidth', ud.size, 'MarkerSize', 8+2*ud.size);
   set(ud.p2, 'LineWidth', ud.size, 'MarkerSize', 8+2*ud.size);

elseif strcmp(action, 'none')
   ud.current = -1;

else
   error(['Unknown command: ' action]);
end

set(ud.h0, 'UserData',ud);

