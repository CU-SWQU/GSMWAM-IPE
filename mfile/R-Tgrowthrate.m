clear
clc

addpath('/glade/work/chunyen/matlabcode')

input = '';
output = '';

if ~exist(output, 'dir')
    mkdir(output)
end

[apxmlat, mlon, apxalt, glon,...
    idx1, idx2, mlatorg, lonorg, altorg, geolat, geolon] = importgrid;

ncid = netcdf.open('/glade/work/chunyen/matlabcode/IPE_Grid.nc');
varid = netcdf.inqVarID(ncid, 'apex_be3');
Be3 = netcdf.getVar(ncid, varid);

apxmlat = apxmlat';
coslm = cosd(apxmlat);
coslmr = repmat(coslm, [1, 80]);

Re = 6371.2e3; % meters
href = 90e3; % meters
R = Re + href;

nu = 2 + sqrt( 4 - 3 * cosd(apxmlat).^2);
de = 3 * cosd(apxmlat);

y = R*log(nu./de) + href;  % meters
%y = apxalt;

wavelength = 500;
kwaveno = 6.283/(wavelength*1.e3);
% ;dyincr in m.  Choose at least 10*pi increments per wavelength.
dyincr = 1000;%3.5e3; %< .2/kwaveno;
hh = .5*dyincr;
ho6 = dyincr/6.;
% ;nincr from 90km to at least 1000km
nincr = fix(910.e3/dyincr) + 1;
% ;yreg,ymid in meters

yreg = (0 : nincr)*dyincr + href;
yunnormal = (yreg-90.e3)/6461.2e3;

fn90 = dir([input 'Cndcty90km*']);
fn150 = dir([input 'Cndcty150km*']);
% %
%  fn90 = dir([input 'conducti*']);
%  fn150 = dir([input 'Cndcty*']);
fnelydn = dir([input 'Elydn*']);

for fnidx = 1 : numel(fn90)

    disp(fnidx)
    % a = h5info([input fn90.name]);
    % integral1
    % integral2   = gravity driven current
    % integral513 = Pedersen d1 (sigma phi)
    % integral514 = Pedersen d2 (sigma lambda)
    % integral517 = hall
    % integral518 = Pedersen e3 (sigma c)
    % integral519 = zonal current (KDm phi)
    % integral520 = meridional current (KDm lambda)

    % from 90 km
    sp1t = h5read([input fn90(fnidx).name], '/apex/integral513');
    sp2t = h5read([input fn90(fnidx).name], '/apex/integral514');
    sht  = h5read([input fn90(fnidx).name], '/apex/integral517');
    sct  = h5read([input fn90(fnidx).name], '/apex/integral518');
    kdmposimt = h5read([input fn90(fnidx).name], '/apex/integral519');
    kdmlt = h5read([input fn90(fnidx).name], '/apex/integral520');
    kgft = h5read([input fn90(fnidx).name], '/apex/integral2');

    sp1t = squeeze(sum(sp1t));
    sp2t = squeeze(sum(sp2t));
    sht = squeeze(sum(sht));
    sct = squeeze(sum(sct));
    kdmposimt = squeeze(sum(kdmposimt));
    kdmlt = squeeze(sum(kdmlt));
    kgft = squeeze(sum(kgft));

    % from 150 km
    sp1t2 = h5read([input fn150(fnidx).name], '/apex/integral513');
    sp2t2 = h5read([input fn150(fnidx).name], '/apex/integral514');
    sht2  = h5read([input fn150(fnidx).name], '/apex/integral517');
    sct2  = h5read([input fn150(fnidx).name], '/apex/integral518');
    kdmposimt2 = h5read([input fn150(fnidx).name], '/apex/integral519');
    kdmlt2 = h5read([input fn150(fnidx).name], '/apex/integral520');
    kgft2 = h5read([input fn150(fnidx).name], '/apex/integral2');

    sp1t2 = squeeze(sum(sp1t2));
    sp2t2 = squeeze(sum(sp2t2));
    sht2 = squeeze(sum(sht2));
    sct2 = squeeze(sum(sct2));
    kdmposimt2 = squeeze(sum(kdmposimt2));
    kdmlt2 = squeeze(sum(kdmlt2));
    kgft2 = squeeze(sum(kgft2));

    % electric field

    % b = h5info([input fnelydn.name]);
    ed1 = h5read([input fnelydn(fnidx).name], '/apex/Ed1');
    ed2 = h5read([input fnelydn(fnidx).name], '/apex/Ed2');
    pot = h5read([input fnelydn(fnidx).name], '/apex/electric_potential');

    % compute total FLI currents for E+F region
    kpmposimt = sp1t.*ed1;
    khmposimt = -(sht-sct).*ed2;
    kpmlt     = -sp2t.*ed2;
    khmlt     = -(sht+sct).*ed1;
    kpmposimt2 = sp1t2.*ed1;
    khmposimt2 = -(sht2-sct2).*ed2;
    kpmlt2     = -sp2t2.*ed2;
    khmlt2     = -(sht2+sct2).*ed1;
    kmposimt  = kpmposimt + khmposimt + kdmposimt + kgft2;

    % compute total FLI currents only for F region
    kmposimt2  = kpmposimt2 + khmposimt2 + kdmposimt2 + kgft2;
    %      kmposimt2  = khmposimt2 +kdmposimt2 + kgft2;
    kmlt      = kpmlt2 + khmlt2 + kdmlt2;

    ex = ed1.*coslmr;
    ey = -ed2.*coslmr;
    kxgf = kgft.*coslmr;
    kxgf2 = kgft2.*coslmr;
    kxt = kmposimt.*coslmr;
    kxt2 = kmposimt2.*coslmr;
    kyt = kmlt.*coslmr;
    
    kpx = kpmposimt2.*coslmr;
    khx = khmposimt2.*coslmr;
    kdm = (kdmposimt2 + kgft2).*coslmr;

    gamma = cell(1, 80);
    ff = cell(1, 80);
    bb = cell(1, 80);
    cc = cell(1, 80);
    dd = cell(1, 80);
    ee = cell(1, 80);
    gg = cell(1, 80);
    t1 = cell(1, 80);
    t2 = cell(1, 80);
    t3 = cell(1, 80);
    for j = 1 : 80

        sp1t19b = sp1t(:, j);
        sp2t19b = sp2t(:, j);
        sht19b  = sht(:, j);
        sct19b  = sct(:, j);
        kxt19p  = kxt2(:, j);
        kyt19p  = kyt(:, j);
        be319   = Be3(:, j);
        kxgf19 = kgft2(:, j);

        kpx19   = kpx(:, j);
        khx19   = khx(:, j);
        kdx19   = kdm(:, j);
    
        kpxreg  = spline(y(8:end), kpx19(8:end), yreg);
        khxreg  = spline(y(8:end), khx19(8:end), yreg);
        kdxreg  = spline(y(8:end), kdx19(8:end), yreg);


	kgxreg = spline(y(8:end), kxgf19(8:end), yreg);

        be3 = Be3(:, j);
        be3reg = spline(y(8:end), be3(8:end), yreg);

%         test = kgxreg(1:211);
% 
%         new = nan(211, 1);
%         for i = 6 : 211-5
% 
%             a = test(i-5:i+5);
%             new(i) = mean(a, 'omitnan');
% 
%         end
%         new(1 : 5) = kgxreg(1:5);
%         new(206 : 211) = kgxreg(206 : 211);
% 
%         new2 = nan(211, 1);
%         for i = 6 : 211-5
% 
%             a = new(i-5:i+5);
%             new2(i) = mean(a, 'omitnan');
% 
%         end
%         new2(1 : 5) = new(1:5);
%         new2(206 : 211) = new(206 : 211);
% 
%         kgxreg = [new2', kgxreg(212:end)];

        n = length(yreg);
        kgxdy = zeros(1, n);
        for i = 2 : n-2

            y0 = kgxreg(i-1);
            yy = kgxreg(i);
            y1 = kgxreg(i+1);
            y2 = kgxreg(i+2);

            x0 = yreg(i-1);
            xx = yreg(i);
            x1 = yreg(i+1);
            x2 = yreg(i+2);

            a = y0*(2*xx-x1-x2)/((x0-x1)*(x0-x2));
            b = y1*(2*xx-x0-x2)/((x0-x1)*(x1-x2));
            c = y2*(2*xx-x0-x1)/((x0-x2)*(x1-x2));
            kgxdy(i) = a-b+c;
        end
        kgxdy(kgxdy==0) = nan;

        f = kgxdy./(kgxreg.*be3reg);

        %         f2 = kgxdy./(kxgf19.*be3);
        %         scale = kgxreg./kgxdy;
        %         scale = scale/1000;
        %         cc = kgxdy./kgxreg;

        sp2reg = spline(y(8:end), sp2t19b(8:end), yreg);
        sp1reg = spline(y(8:end), sp1t19b(8:end), yreg);
        shtreg = spline(y(8:end), sht19b(8:end), yreg);
        sctreg = spline(y(8:end), sct19b(8:end), yreg);
        kxtfreg = spline(y(8:end), kxt19p(8:end), yreg);
        kytfrep = spline(y(8:end), kyt19p(8:end), yreg);
        fkxtfreg = f.*kxtfreg;
        fkytfreg = f.*kytfrep;
        fkp = f.*kpxreg;
        fkh = f.*khxreg;
        fkd = f.*kdxreg;
        garr = fkxtfreg./sp1reg;
        
        garrp = fkp./sp1reg;
        garrh = fkh./sp1reg;
        garrkdx = fkd./sp1reg;
        gamma{j} = garr';
        ff{j} = f';
        bb{j} = kgxdy';
        cc{j} = kgxreg';

        dd{j} = kpxreg';

        ee{j} = be3reg';
   
        t1{j} = garrp';
        t2{j} = garrh';
        t3{j} = garrkdx';

        gg{j} = sp1reg';

    end
    gamma = cell2mat(gamma);
    factor = cell2mat(ff);
    factor = factor(:, idx2);

    sp1reg = cell2mat(gg);
    sp1reg = sp1reg(:, idx2);

    be3reg = cell2mat(ee);
    be3reg = be3reg(:, idx2);

    kgxdy = cell2mat(bb);
    kgxdy = kgxdy(:, idx2);

    kgxreg = cell2mat(cc);
    kgxreg = kgxreg(:, idx2);

    kpxreg = cell2mat(dd);
    kpxreg = kpxreg(:, idx2);

    garrkpx = cell2mat(t1);
    garrkdx = cell2mat(t3);
    garrkhx = cell2mat(t2);

    garrkpx = garrkpx(:, idx2);
    garrkdx = garrkdx(:, idx2);
    garrkhx = garrkhx(:, idx2);
    
    %     gamma(gamma<=0) = nan;
    %     u = find(yreg/1000<150);
    %     gamma(u, :) = nan;
    %     contourf(glon, yreg/1000, gamma, 50, 'linestyle', 'none')

    gamma = gamma(:, idx2);

    %     ff = ff(:, idx2);
    %     contourf(glon, yreg(62:end)/1000, ff(62:end, :), 30, 'linestyle', 'none')
    %     clim([-0.001 0.001])
    %     ylim([90 300])

    filename = ['garr.apex.' fnelydn(fnidx).name(12:23) '.mat'];

    save([output filename], 'gamma', 'yreg', 'kgxreg', 'kgxdy',...
        'kxt', 'kxt2', 'kyt', 'kxgf', 'kxgf2', 'kpxreg', 'sp1reg',...
        'garrkpx', 'garrkhx', 'garrkdx', 'factor', 'be3reg')
end

% 
% for i = 22
%     figure
%     subplot(221)
%     plot(kgxreg(:, i), yreg/1000, 'k', 'linewidth', 1.5)
%     ylim([0 300])
%     title('Kgx')
%     set(gca, 'fontsize', 15)
%     ylabel('hA (km)')
%     xlabel('Kgx (A/m)')
%     grid on
% 
%     subplot(222)
%     plot(kgxdy(:, i), yreg/1000, 'k', 'linewidth', 1.5)
%     ylim([0 300])
%     title('Kgx/dy')
%     set(gca, 'fontsize', 15)
%      ylabel('hA (km)')
%     xlabel('Kgx/dy (A/m^2)')
% grid on
%     subplot(223)
%     plot(gamma(:, i), yreg/1000, 'k', 'linewidth', 1.5)
%     ylim([0 300])
%     title('\gamma')
%     set(gca, 'fontsize', 15)
%        ylabel('hA (km)')
%     xlabel('growth rate (sec^-^1)')
%     grid on
%     xlim([0 8e-4])
% end
% set(gcf, 'color', 'w')

% gamma(gamma<0) = nan;
% contourf(glon, yreg/1000, gamma, 100, 'linestyle', 'none')
% ylim([150 1000])
% colormap jet


% plot kgxreg and zoom in to 100-500 km
% before and after the smoothing process

