function [mlat_final, mlon, apexA, geolont, idx1, idx2, mlat, lon, alt, geolat, geolon] = importgrid()

ncfile = '/glade/work/chunyen/matlabcode/IPE_Grid.nc';
h5file = '/glade/work/chunyen/matlabcode/IPE_Grid.h5';

colat = ncread(ncfile, 'm_colat');
lon = ncread(ncfile, 'longitude');
alt = ncread(ncfile, 'altitude');
alt = alt/1000;

mlat = 90-colat*180/pi;
mlon = 0:4.5:355.5;
mlat_final = mlat(1,:);

apexA = nan(170,1);
idx1 = nan(170,1);
for i = 1 : 170
    [apexA(i), idx1(i)] = max(alt(:,i));
end

geolon = nan(170,80);
for j = 1 : 170
    geolon(j, :) = lon(idx1(j), j, :);
end
geolon = geolon * 180/pi;
geolon = geolon(end, :);
[geolont, idx2] = sort(geolon);

lon = lon*180/pi;
cogeolat = h5read(h5file, '/apex_grid/colatitude');
geolat = 90-cogeolat*180/pi;
