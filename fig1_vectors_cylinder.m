clear all;
figure;

% parameters
slice = 750;
arrowScale = 0.07;
    
% read data
f = flow('data/cylinder2d.nc');
[U,V] = f.components();
[X,Y] = f.grid();

% select subset of arrows
sub = 6;
X = X(1:sub:end, 1:sub:end);
Y = Y(1:sub:end, 1:sub:end);
U = U(1:sub:end, 1:sub:end, slice);
V = V(1:sub:end, 1:sub:end, slice);

% plot reflectance
R = f.Reflectance(:,:,slice);
x = linspace(f.AxisX(1), f.AxisX(2), size(R,1));
y = linspace(f.AxisY(1), f.AxisY(2), size(R,2));
imagesc(x, y, R');
colormap(gray());
caxis([0 1]);
hold all;

% arrow plot
q = quiver(X, Y, U*arrowScale, V*arrowScale, 'yellow', 'LineWidth', 1.2);
q.AutoScale = 'off';
axis equal;
eps = 0.1;
axis([min(X(:))-eps max(X(:))+eps min(Y(:))-eps max(Y(:))+eps]);

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'white');