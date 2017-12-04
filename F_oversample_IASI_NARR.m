function output_oversample = F_oversample_IASI_NARR(inp,output_subset)
% very similar to F_regrid_IASI_NARR.m, but instead of using 2-D Gaussian
% SRF, it uses the neighbors-within-a-circle oversampling

% Modified from F_regrid_IASI_NARR.m by Kang Sun on 2017/11/28

output_oversample = [];

if ~isfield(inp,'NARR_download_dir')
    error('You have to know where NARR data are!')
else
    NARR_download_dir = inp.NARR_download_dir;
end

wd_bin = inp.wd_bin;

if ~isfield(inp,'ws_bin')
    ws_bin = [0 inf];
else
    ws_bin = inp.ws_bin;
end

if ~isfield(inp,'T_bin')
    T_bin = [0 inf];
else
    T_bin = inp.T_bin;% in K
end

nwd = length(wd_bin);
nws = length(ws_bin)-1;
nT = length(T_bin)-1;

which_wind = inp.which_wind;

% IASI pixel area in km2/pi
A_IASI = [176.112319946289,175.997787475586,197.295578002930,197.296386718750,134.232391357422,134.307479858398,147.597915649414,147.511016845703,106.910308837891,106.910293579102,115.835662841797,115.835952758789,88.1102142333984,88.1099700927734,94.3334960937500,94.3336944580078,74.6136627197266,74.6135177612305,79.0981597900391,79.0983047485352,64.6809158325195,64.6807861328125,68.0252914428711,68.0253448486328,57.2129020690918,57.1962089538574,59.7423744201660,59.7241744995117,51.4984931945801,51.5112380981445,53.4395408630371,53.4396286010742,47.1394691467285,47.1294631958008,48.6127166748047,48.6127738952637,43.7793922424316,43.7793731689453,44.9149169921875,44.9149055480957,41.2354736328125,41.2354354858398,42.0903854370117,42.0904350280762,39.3517951965332,39.3517684936523,39.9780311584473,39.9780654907227,38.0385055541992,38.0384712219238,38.4638633728027,38.4605789184570,37.2202949523926,37.2202720642090,37.4690208435059,37.4670410156250,36.8619995117188,36.8632621765137,36.9433822631836,36.9425735473633,36.9458961486816,36.9475517272949,36.8649864196777,36.8641242980957,37.4771308898926,37.4791526794434,37.2297134399414,37.2281494140625,38.4804687500000,38.4804649353027,38.0526008605957,38.0526084899902,40.0019073486328,40.0066566467285,39.3770217895508,39.3770751953125,42.1226539611816,42.1291275024414,41.2644309997559,41.2644348144531,44.9572257995606,44.9657287597656,43.8254928588867,43.8177566528320,48.6675682067871,48.6785697937012,47.1893501281738,47.1893539428711,53.5103225708008,53.5103416442871,51.5756797790527,51.5885810852051,59.8340911865234,59.8340873718262,57.2961006164551,57.3128280639648,68.1210937500000,68.1691513061523,64.7891693115234,64.7891540527344,79.2890853881836,79.2891082763672,74.7851409912109,74.7851104736328,94.5494537353516,94.5494995117188,88.3023147583008,88.3407821655273,116.077011108398,116.198059082031,107.175537109375,107.228660583496,147.949081420898,148.125152587891,134.687866210938,134.687728881836,198.105438232422,198.105773925781,176.683227539063,176.683029174805];

Startdate = inp.Startdate;
Enddate = inp.Enddate;

res = inp.res;
max_x = inp.max_x;
max_y = inp.max_y;
clon = inp.clon;
clat = inp.clat;
R = inp.R;

if ~isfield(inp,'do_weight')
    do_weight = false;
else
    do_weight = inp.do_weight;
end

max_lon = clon+max_x*1.2/110/cos((abs(clat)+max_y/111)/180*pi);
min_lon = clon-max_x*1.2/110/cos((abs(clat)+max_y/111)/180*pi);
max_lat = clat+max_y*1.2/110;
min_lat = clat-max_y*1.2/110;

% define x y grids
xgrid = (-max_x+0.5*res):res:max_x;
ygrid = (-max_y+0.5*res):res:max_y;
nrows = length(ygrid);
ncols = length(xgrid);

% define x y mesh
[xmesh, ymesh] = meshgrid(single(xgrid),single(ygrid));

f1 = output_subset.utc >= single(datenum([Startdate 0 0 0])) ...
    & output_subset.utc <= single(datenum([Enddate 23 59 59]));
f2 = output_subset.lat >= min_lat-0.5 & output_subset.lat <= max_lat+0.5...
    & output_subset.lon >= min_lon-0.5 & output_subset.lon <= max_lon+0.5;
f3 = ~isnan(output_subset.colnh3);
f4 = ~isnan(output_subset.colnh3error);

validmask = f1&f2&f3&f4;

nL2 = sum(validmask);
disp([num2str(nL2),' pixels to be regridded...'])

Lon = output_subset.lon(validmask);
Lat = output_subset.lat(validmask);
Ifov = output_subset.ifov(validmask);
ColNH3 = output_subset.colnh3(validmask);
ColNH3e = output_subset.colnh3error(validmask);
UTC = output_subset.utc(validmask);

[UTC, I] = sort(UTC);
Lat = Lat(I);
Lon = Lon(I);
Ifov = Ifov(I);
ColNH3 = ColNH3(I);
ColNH3e = ColNH3e(I);

% allocate space for wind, associated with each satellite pixel
U_sat_pbl = nan(size(Lat),'single');
V_sat_pbl = U_sat_pbl;
T_sat_pbl = U_sat_pbl;

% allocate space for wind, associated with each rotation center
U_c_pbl = nan(size(Lat),'single');
V_c_pbl = U_sat_pbl;
T_c_pbl = U_sat_pbl;

day_1 = datevec(UTC(1));
day_end = datevec(UTC(end));

% saved NARR data in which hours
narr_hour_vec = 12:3:24;
% work out NARR's map projection
mstruct_narr = defaultm('lambertstd');
mstruct_narr.origin = [50 360-107 0];
mstruct_narr.mapparallels = 50;
mstruct_narr.nparallels = 1;
mstruct_narr.scalefactor = 6367470/1e3;
mstruct_narr.falseeasting = 0;
mstruct_narr.falsenorthing = 0;
mstruct_narr = defaultm(mstruct_narr);
% map center lat lon to narr projection x y
[x_c, y_c] = mfwdtran(mstruct_narr,clat,clon);

totalindex = 1:nL2;

% work out the wind interpolation day-by-day
disp('Interpolating NARR wind to satellite pixel locations...')

for narr_date = datenum(day_1(1:3)):1:datenum(day_end(1:3))
    int = UTC > datenum(narr_date) & UTC <= datenum(narr_date)+1;
    if sum(int) > 0
        [N,~,bin] = histcounts(UTC(int),narr_date+[narr_hour_vec-1.5,narr_hour_vec(end)+1.5]/24);
        lon_interp_day = Lon(int);
        lat_interp_day = Lat(int);
        idx_interp_day = totalindex(int);
        for ibin = 1:length(N)
            if N(ibin) > 0
                tmp = datevec(narr_date+narr_hour_vec(ibin)/24);
                narr_year = num2str(tmp(1));
                narr_month = num2str(tmp(2),'%02d');
                narr_day = num2str(tmp(3),'%02d');
                narr_hour = num2str(tmp(4),'%02d');
                disp(['Interpolating NARR for ',datestr(tmp)])
                inp_interp_narr = [];
                inp_interp_narr.narr_data_3d = load([NARR_download_dir,'/',narr_year,'/',narr_month,'/',...
                    'subset_3D_',narr_day,'_',narr_hour,'.mat']);
                inp_interp_narr.narr_data_sfc = load([NARR_download_dir,'/',narr_year,'/',narr_month,'/',...
                    'subset_sfc_',narr_day,'_',narr_hour,'.mat']);
                inp_interp_narr.lon_interp_hour = lon_interp_day(bin == ibin);
                inp_interp_narr.lat_interp_hour = lat_interp_day(bin == ibin);
                inp_interp_narr.x_c = x_c;
                inp_interp_narr.y_c = y_c;
                inp_interp_narr.mstruct_narr = mstruct_narr;
                inp_interp_narr.max_x = max_x;
                inp_interp_narr.max_y = max_y;
                inp_interp_narr.P_pblmax = 500;% maximum pbl thickness in hPa
                outp_interp_narr = F_interp_narr(inp_interp_narr);
                
                idx_interp_hour = idx_interp_day(bin == ibin);
                U_sat_pbl(idx_interp_hour) = outp_interp_narr.u_sat_pbl;
                V_sat_pbl(idx_interp_hour) = outp_interp_narr.v_sat_pbl;
                T_sat_pbl(idx_interp_hour) = outp_interp_narr.T_sat_pbl;
                U_c_pbl(idx_interp_hour) = outp_interp_narr.u_c_pbl;
                V_c_pbl(idx_interp_hour) = outp_interp_narr.v_c_pbl;
                T_c_pbl(idx_interp_hour) = outp_interp_narr.T_c_pbl;
            end
        end
    end
end

disp('Converting pixel lat lon to x y in km...')
% call function F_latlon2xy to convert lat lon to x y.
inp_xy = [];
inp_xy.clon = clon;
inp_xy.clat = clat;
inp_xy.lon = Lon(:);
inp_xy.lat = Lat(:);
outp_xy = F_latlon2xy(inp_xy);

disp('Categorize data based on WS/WD/T...')
switch which_wind
    case 'pbl_pixel'
        u_vec_rot = U_sat_pbl;
        v_vec_rot = V_sat_pbl;
    case 'pbl_center'
        u_vec_rot = U_c_pbl;
        v_vec_rot = V_c_pbl;
    case 'pbl_average'
        u_vec_rot = 0.5*(U_sat_pbl+U_c_pbl);
        v_vec_rot = 0.5*(V_sat_pbl+V_c_pbl);
end

ws_vec_rot = sqrt(u_vec_rot.^2+v_vec_rot.^2);
wd_vec_rot = u_vec_rot;
wd_vec_rot(v_vec_rot >= 0) = acos(u_vec_rot(v_vec_rot >= 0)./ws_vec_rot(v_vec_rot >= 0));
wd_vec_rot(v_vec_rot < 0) = 2*pi-acos(u_vec_rot(v_vec_rot < 0)./ws_vec_rot(v_vec_rot < 0));

[~,~,ws_idx] = histcounts(ws_vec_rot,ws_bin);
[~,~,wd_idx] = histcounts(wd_vec_rot,[0 wd_bin 2*pi]);
wd_idx(wd_idx == nwd+1) = 1;
[~,~,T_idx] = histcounts(T_sat_pbl,T_bin);

O = zeros(nwd,nws,nT,nrows,ncols,'single');
D = zeros(nwd,nws,nT,nrows,ncols,'single');

O_r = zeros(nwd,nws,nT,nrows,ncols,'single');
D_r = zeros(nwd,nws,nT,nrows,ncols,'single');

disp([num2str(nL2),' pixels to be regridded...'])

for iwd = 1:nwd
    for iws = 1:nws
        for iT = 1:nT
            use_idx = wd_idx == iwd & ws_idx == iws & T_idx == iT;
            nl2 = sum(use_idx);
            disp([num2str([iwd iws iT]),' has ',num2str(nl2),' pixels'])
            if nl2 > 0
                x_c_inp = outp_xy.x(use_idx);
                y_c_inp = outp_xy.y(use_idx);
                
                %ws_inp = ws_vec_rot(use_idx);
                wd_inp = wd_vec_rot(use_idx);
                
                vcd_inp = ColNH3(use_idx);
                vcde_inp = ColNH3e(use_idx);
                xtrack_inp = Ifov(use_idx);
                
                x_cr_inp = x_c_inp;
                y_cr_inp = y_c_inp;
                
                % rotate pixel according to wind direction
                for il2 = 1:nl2
                   x_c = x_c_inp(il2);
                   y_c = y_c_inp(il2);
                   
%                    ws = ws_inp(il2);
                   wd = wd_inp(il2);
                   
                   xy_rot = [cos(wd) sin(wd);-sin(wd) cos(wd)]*[x_c;y_c];
%                    x_rr = xy_rot(1,1:4);
%                    y_rr = xy_rot(2,1:4);
                   x_cr_inp(il2) = xy_rot(1,:);
                   y_cr_inp(il2) = xy_rot(2,:);
                end
                
                tmpO = zeros(nrows,ncols,'single');
                tmpD = zeros(nrows,ncols,'single');
                
                tmpOr = zeros(nrows,ncols,'single');
                tmpDr = zeros(nrows,ncols,'single');
                % sorry for many nested loops
                for irow = 1:nrows
                    for icol = 1:ncols
                        xcenter = xgrid(icol);ycenter = ygrid(irow);
                        
                        % for each grid center, find the closeby unrotated
                        % pixels
                        x = x_c_inp; y = y_c_inp;
                        tmp_ind = find(abs(x-xcenter) <= R & abs(y-ycenter) <= R);
                        Distance = sqrt((x(tmp_ind)-xcenter).^2+(y(tmp_ind)-ycenter).^2);
                        ind = tmp_ind(Distance <= R);
                        if do_weight
                            weight = 1./vcde_inp(ind)./(A_IASI(xtrack_inp(ind)));
                            weight = weight/nansum(weight);
                            tmpO(irow,icol) = nansum(vcd_inp(ind).*weight);
                        else
                            tmpO(irow,icol) = nanmean(vcd_inp(ind));
                        end
                        tmpD(irow,icol) = length(ind);
                        
                        % for each grid center, find the closeby rotated
                        % pixels
                        x = x_cr_inp; y = y_cr_inp;
                        tmp_ind = find(abs(x-xcenter) <= R & abs(y-ycenter) <= R);
                        Distance = sqrt((x(tmp_ind)-xcenter).^2+(y(tmp_ind)-ycenter).^2);
                        ind = tmp_ind(Distance <= R);
                        if do_weight
                            weight = 1./vcde_inp(ind)./(A_IASI(xtrack_inp(ind)));
                            weight = weight/nansum(weight);
                            tmpOr(irow,icol) = nansum(vcd_inp(ind).*weight);
                        else
                            tmpOr(irow,icol) = nanmean(vcd_inp(ind));
                        end
                        tmpDr(irow,icol) = length(ind);
                    end
                    
                end
                          
                O(iwd,iws,iT,:,:) = tmpO;
                D(iwd,iws,iT,:,:) = tmpD;
                
                O_r(iwd,iws,iT,:,:) = tmpOr;
                D_r(iwd,iws,iT,:,:) = tmpDr;
            end
        end
    end
end

output_oversample.O = O;
output_oversample.D = D;

output_oversample.O_r = O_r;
output_oversample.D_r = D_r;

output_oversample.xgrid = xgrid;
output_oversample.ygrid = ygrid;
output_oversample.xmesh = single(xmesh);
output_oversample.ymesh = single(ymesh);

disp('Transfer x y coordinates back to lat lon...')
[tmplat,tmplon] = minvtran(outp_xy.mstruct,xmesh(:),ymesh(:));
output_oversample.latmesh = single(reshape(tmplat,size(xmesh)));
output_oversample.lonmesh = single(reshape(tmplon,size(xmesh)));

output_oversample.mstruct = outp_xy.mstruct;
output_oversample.max_lon = max_lon;
output_oversample.min_lon = min_lon;

output_oversample.max_lat = max_lat;
output_oversample.min_lat = min_lat;
return

function outp_interp_narr = F_interp_narr(inp_interp_narr)
% matlab function to interpolate narr met field to given satellite pixel
% locations (lat_interp_hour, lon_interp_hour) and central locations (x_c,
% y_c, which should already been transformed to km using narr projection)

% outputs are surface and pbl-integrated meteorology at satellite locations
% and central location

% written by Kang Sun on 2017/11/09

narr_data_3d = inp_interp_narr.narr_data_3d;
narr_data_sfc = inp_interp_narr.narr_data_sfc;
lon_interp_hour = inp_interp_narr.lon_interp_hour;
lat_interp_hour = inp_interp_narr.lat_interp_hour;
x_c = inp_interp_narr.x_c;
y_c = inp_interp_narr.y_c;
mstruct_narr = inp_interp_narr.mstruct_narr;
max_x = inp_interp_narr.max_x;
max_y = inp_interp_narr.max_y;

intx = narr_data_3d.x > x_c - max_x-500 & narr_data_3d.x < x_c + max_x+500;
inty = narr_data_3d.y > y_c - max_y-500 & narr_data_3d.y < y_c + max_y+500;

x = narr_data_3d.x(intx);
y = narr_data_3d.y(inty);

u = narr_data_3d.u(:,inty,intx);
v = narr_data_3d.v(:,inty,intx);
T = narr_data_3d.T(:,inty,intx);

P_surf = narr_data_sfc.P_surf(inty,intx);
T_surf = narr_data_sfc.T_surf(inty,intx);
PBLH = narr_data_sfc.PBLH(inty,intx);

[x_sat, y_sat] = mfwdtran(mstruct_narr,lat_interp_hour,lon_interp_hour);

% pressure grid of narr data, in Pa
P = [550 600 650 700 725 750 775 800 825 850 875 900 925 950 975 1000]'*100;
% maximum pbl thickness in Pa
if ~isfield(inp_interp_narr,'P_pblmax')
    P_pblmax = 200*100;
else
    P_pblmax = inp_interp_narr.P_pblmax*100;% from hPa to Pa
end
%% work out the pbl pressure weighted wind
% % find the surface layer pressure index, not worth it
% sfc_idx = zeros(size(narr_data_sfc.P_surf));
% for i = length(P):-1:1
%     mask = narr_data_sfc.P_surf-P(i) > 0 & sfc_idx == 0;
%     sfc_idx(mask) = i;
%     if sum(mask(:)) == 0;
%         break;
%     end
% end
nx = sum(intx);
ny = sum(inty);
u_pbl = nan(ny,nx,'single');
u_sfc = nan(ny,nx,'single');
v_pbl = nan(ny,nx,'single');
v_sfc = nan(ny,nx,'single');
T_pbl = nan(ny,nx,'single');

% pbl top pressure
P_pbl = nan(ny,nx,'single');
for ix = 1:nx
    for iy = 1:ny
        scaleH = 287/9.8*(T_surf(iy,ix)-10);
        P_pbl(iy,ix) = max([P_surf(iy,ix)*exp(-PBLH(iy,ix)/scaleH) P_surf(iy,ix)-P_pblmax]);
        
        int = P >= P_pbl(iy,ix) & P <= P_surf(iy,ix) & P >= P_surf(iy,ix)-P_pblmax;
        % sometimes the pbl is so thin...
        if sum(int) == 0
            int = find(P <= P_surf(iy,ix),1,'last');
        end
        %         sum(int)
        localu = squeeze(u(int,iy,ix));
        localv = squeeze(v(int,iy,ix));
        u_pbl(iy,ix) = sum(localu.*P(int)/sum(P(int)));
        v_pbl(iy,ix) = sum(localv.*P(int)/sum(P(int)));
        
        u_sfc(iy,ix) = localu(end);
        v_sfc(iy,ix) = localv(end);
        
        T_pbl(iy,ix) = sum(squeeze(T(int,iy,ix)).*P(int)/sum(P(int)));
    end
end
%%
outp_interp_narr = [];
outp_interp_narr.u_sat_pbl = interp2(x,y,u_pbl,x_sat,y_sat);
outp_interp_narr.v_sat_pbl = interp2(x,y,v_pbl,x_sat,y_sat);
outp_interp_narr.T_sat_pbl = interp2(x,y,T_pbl,x_sat,y_sat);

nz = length(outp_interp_narr.u_sat_pbl);

outp_interp_narr.u_c_pbl = ones(nz,1,'single')*interp2(x,y,u_pbl,x_c,y_c);
outp_interp_narr.v_c_pbl = ones(nz,1,'single')*interp2(x,y,v_pbl,x_c,y_c);
outp_interp_narr.T_c_pbl = ones(nz,1,'single')*interp2(x,y,T_pbl,x_c,y_c);

outp_interp_narr.u_sat_sfc = interp2(x,y,u_sfc,x_sat,y_sat);
outp_interp_narr.v_sat_sfc = interp2(x,y,v_sfc,x_sat,y_sat);
outp_interp_narr.T_sat_sfc = interp2(x,y,T_surf,x_sat,y_sat);

outp_interp_narr.u_c_sfc = ones(nz,1,'single')*interp2(x,y,u_sfc,x_c,y_c);
outp_interp_narr.v_c_sfc = ones(nz,1,'single')*interp2(x,y,v_sfc,x_c,y_c);
outp_interp_narr.T_c_sfc = ones(nz,1,'single')*interp2(x,y,T_surf,x_c,y_c);
return
