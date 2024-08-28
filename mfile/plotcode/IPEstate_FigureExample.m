clear
clc

[apexlat, mlon, apex_height, glon_meq2, idx1, idx2, mlat, glon3d, alt, glat3d, glon_meq1] = importgrid;

% The size of IPE parameter is 1115x170x80. The three dimensions represent
% magnetic point along a field line, magnetic latitude (qd-lat), and magnetic
% longitude.

filename = 'IPE_State.apex.201303180000.h5';
a = h5info(filename);
op = h5read(filename, '/apex/o_plus_density');

%% Example 1: QDlat-hA plot

longrid = 80;  % 284E
x = mlat(:);
y = alt(:);
z = op(:, :, longrid);
z = z(:);

scatter(x, y, 20, z, 'filled')
box on
cb = colorbar;
colormap jet
ylim([100 1000])
title('IPE O^+ density at 284E')
xlabel('Quasi-Dipole latitude')
ylabel('apex height (km)')
ylabel(cb, 'm^-^3')
set(gca, 'fontsize', 20)

%% Example 2: geolon-geolat plot

h = 302;  % at 302 km
alt3d = repmat(alt, [1 1 80]); % the size should be same as glon3d (glat3d)
u = find(alt3d == h);

x = glon3d(u);
y = glat3d(u);
z = op(u);

scatter(x, y, 20, z, 'filled')
box on
cb = colorbar;
colormap jet
ylim([-90 90])
xlim([0 360])
title('IPE O^+ density at 302 km')
xlabel('geographic longitude')
ylabel('geographic latitude')
ylabel(cb, 'm^-^3')
set(gca, 'fontsize', 25)
set(gca, 'xtick', 0:60:360, 'ytick', -90:30:90)

%% Example 3: conductance profile

filename90km = 'Elydn90km.apex.201303180000.h5';
filename150km = 'Elydn150km.apex.201303180000.h5';
b = h5info(filename90km);
% the definition of integrals 513~520 can be found in Richmond.(1995), JGR,
% and integrals 1 and 2 are Equation 28 and 31 in Huang et al. (2024), JGR.
% The two dimensions of the matrix are hemisphere, magnetic latitude 
% (qd-lat), and magnetic longitude.

sigmay90km = h5read(filename90km, '/apex/integral514');
sigmay90km = squeeze(sum(sigmay90km));

sigmay150km = h5read(filename150km, '/apex/integral514');
sigmay150km = squeeze(sum(sigmay150km));

x1 = sigmay90km(:, longrid); 
x2 = sigmay150km(:, longrid);
y = apex_height;

plot(x1, y, 'k', 'linewidth', 2)
hold on
grid on
plot(x2, y, 'r--', 'linewidth', 2)
ylim([90 1000])
xlabel('\Sigma_y (S)')
ylabel('apex height')
set(gca, 'fontsize', 20)
title('Pedersen Conductance at 284E at UT 0 ')
legend({'from 90 km'; 'from 150 km'})



