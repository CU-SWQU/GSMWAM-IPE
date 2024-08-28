function [apexlat, mlon, apex_height, glon_meq2, idx1, idx2, mlat, glon3d, alt, glat3d, glon_meq1] = importgrid()

% The size of IPE parameter is 1115x170x80. The three dimensions represent
% magnetic point along a field line, magnetic latitude (qd-lat), magnetic
% longitude.
%
% apexlat: qd-lat
% mlon: qd-lon
% apex_height: magnetic fied-line apex point
% glon_meq2: geographic longitude at magnetic equator (0~360E)
% idx1: apex point index
% idx2: apex point index
% mlat: magnetic latitude
% glon3d: 3d geographic longitude
% alt: 2d altitude (all field lines should be the same)
% glat3d: 3d geographic latitude
% glon_meq1: geographic longitude at magnetic equator (same as the glon_meq2, 
% but haven't sorted)


h5file = 'D:\WAMIPE\IPE_Grid.h5';

colat = h5read(h5file, '/apex_grid/m_colat');
lon   = h5read(h5file, '/apex_grid/longitude');
alt   = h5read(h5file, '/apex_grid/altitude');
alt = alt/1000;

mlat = 90-colat*180/pi;
mlon = 0:4.5:355.5;
apexlat = mlat(1,:);

apex_height = nan(170,1);
idx1 = nan(170,1);
for i = 1 : 170
    [apex_height(i), idx1(i)] = max(alt(:,i));
end

glon_meq1 = nan(170,80);
for j = 1 : 170
    glon_meq1(j, :) = lon(idx1(j), j, :);
end
glon_meq1 = glon_meq1 * 180/pi;
glon_meq1 = glon_meq1(end, :);
[glon_meq2, idx2] = sort(glon_meq1);

glon3d = lon*180/pi;
cogeolat = h5read(h5file, '/apex_grid/colatitude');
glat3d = 90-cogeolat*180/pi;