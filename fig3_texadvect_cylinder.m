clear all;
figure;

% parameters
res = [2560, 320];
stepSize = 0.01;
duration = 1.0;
slice = 1000;

% read data
f = flow('data/cylinder2d.nc');

% create texture and advect
tex = zeros(res);
for ix=1:res(1)
    for iy=1:res(2) 
        tex(ix,iy) = mod(floor(ix/40) + floor(iy/40), 2) * 0.25 + 0.5;
    end
end
adv = f.texadvect(stepSize, duration, tex, res, slice);
    
% drawing
x = linspace(f.AxisX(1), f.AxisX(2), res(1));
y = linspace(f.AxisY(1), f.AxisY(2), res(2));
imagesc(x, y, adv');
set(gca,'YDir','normal');
colormap('gray');
caxis([0 1])
axis equal;
axis([f.AxisX(1) f.AxisX(2) f.AxisY(1) f.AxisY(2)]);
hold all;

% draw contour
plot(f.Contour(:,1), f.Contour(:,2), 'black');