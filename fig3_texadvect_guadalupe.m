clear all;
figure;

% parameters
res = [700, 640];   % resolution of the output texture
stepSize = 100;     % integration step size in seconds
duration = 10000;   % integration duration in seconds
slice = 70;         % time slice

% read data
f = flow('data/guadalupe.nc');

% create texture and advect
tex = zeros(res);
for ix=1:res(1)
    for iy=1:res(2) 
        tex(ix,iy) = mod(floor(ix/10) + floor(iy/10), 2) * 0.25 + 0.5;
    end
end
adv = f.texadvect(stepSize, duration, tex, res, slice);
ilat = imresize(f.Latitude, size(adv));
ilon = imresize(f.Longitude, size(adv));
worldmap(f.AxisY, f.AxisX);
geoshow(ilat, ilon, double(adv));

% draw contour
geoshow(f.Contour(:,2), f.Contour(:,1), 'color', [0 0 0]);