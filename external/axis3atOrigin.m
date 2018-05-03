function r = axis3atOrigin(limits)
% Plots 3D axes through the origin.
% Function has be adapted from Plot3AxisAtOrigin.m by Michael Robbins [1].
%
% authors:
%   Michael Robbins (michael.robbins@us.cibc.com, robbins@bloomberg.net)
%   Matthias Wolff, BTU Cottbus-Senftenberg (some modifications to my taste :)
%
% references:
% [1] https://de.mathworks.com/matlabcentral/fileexchange/3245-plot3axisatorigin
%     visited Mar. 30, 2018.

axis off;
hold on;

% GET AXIS START/END AND OFFSETS
x0 = limits(1);
x1 = limits(2);
y0 = limits(3);
y1 = limits(4);
z0 = limits(5);
z1 = limits(6);
dx = (x1-x0)/50;
dy = (y1-y0)/50;
dz = (z1-z0)/50;

% DRAW AXIS LINEs
xlim([x0 x1]);
ylim([y0 y1]);
zlim([z0 z1]);
plot3([x0 x1],[0 0],[0 0],'k');
plot3([0 0],[y0 y1],[0 0],'k');
plot3([0 0],[0 0],[z0 z1],'k');

% DRAW AXES ARROWS
plot3([x1-dx x1 x1-dx], [0 0 0], [-dz/2 0 dz/2], 'k');
plot3([0 0 0], [y1-dy y1 y1-dy], [-dz/2 0 dz/2], 'k');
plot3([-dx/2.5 0 dx/2.5], [0 0 0], [z1-1.5*dz z1 z1-1.5*dz], 'k');

% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');
Z=get(gca,'Ztick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');
ZL=get(gca,'ZtickLabel');

% REMOVE TICKS
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'Ztick',[]);

% DRAW TICKS
%%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
for i=1:length(X)
  plot3([X(i) X(i)],[0 0],[-dz dz],'k');
end
for i=1:length(Y)
  plot3([0 0],[Y(i) Y(i)],[-dz dz],'k');
end
for i=1:length(Z)
  plot3([-dx/1.5 dx/1.5],[0 0],[Z(i) Z(i)],'k');
end

% DRAW LABELS
text(X,zeros(size(X)),zeros(size(X))-1.5*dz,XL, ...
  'HorizontalAlignment','center','VerticalAlignment','top');
text(zeros(size(Y)),Y,zeros(size(Y))-1.5*dz,YL, ...
  'HorizontalAlignment','center','VerticalAlignment','top');
text(zeros(size(Z))-1.5*dx,zeros(size(Z)),Z,ZL, ...
  'HorizontalAlignment','right','VerticalAlignment','middle');

% DRAW AXES LABELS
r.xlabel = text(x1,0,-1.5*dz,'x',...
  'HorizontalAlignment','left','VerticalAlignment','top');
r.ylabel = text(0,y1,-1.5*dz,'y',...
  'HorizontalAlignment','right','VerticalAlignment','top');
r.zlabel = text(-1.5*dx,0,z1,'z',...
  'HorizontalAlignment','right','VerticalAlignment','bottom');

% EOF
