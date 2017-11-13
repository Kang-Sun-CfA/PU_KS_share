function output_regrid = F_regrid_IASI_NARR(inp,output_subset)
% function to take in the output from F_subset_IASI.m and regrid these L2
% data to a L3 grid, centered at clon and clat and bounded by max_x and 
% max_y with resolution res in km.

% ... above is the description of F_regrid_IASI_km.m. This function is far
% more powerful. It also reads the NARR met data, finds the nearest time of
% satellite observation, and interpolate NARR met field to OMI pixel
% locations.

% As a result you need to know where the NARR data are. They were saved to
% SAO unix systems using a proprietary code fetch_NARR_monthlydownload.m

% then it categorize omi pixel according to wind direction, wind speed, and
% pixel pbl temperature

% it also rotates the omi pixel and align the wind direction. both the
% non-rotated and rotated sums are output.


% Modified from F_regrid_IASI_km.m by Kang Sun on 2017/11/13

output_regrid = [];

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

% the pixel geometry is much more elegant than the lat lon version
u_km = inp.u_km;
v_km = inp.v_km;
t_km = inp.t_km;
pixel_left = inp.pixel_left;
pixel_down = inp.pixel_down;

Startdate = inp.Startdate;
Enddate = inp.Enddate;

res = inp.res;
max_x = inp.max_x;
max_y = inp.max_y;
clon = inp.clon;
clat = inp.clat;

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

% construct a rectangle envelopes the orginal pixel
xmargin = 3; % how many times to extend zonally
ymargin = 2; % how many times to extend meridonally

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
Lat = Lat(I,:);
Lon = Lon(I,:);
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

for narr_date = datenum(day_1(1:3)):1:datenum(day_end(1:3));
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
                inp_interp_narr.P_pblmax = 100;% maximum pbl thickness in hPa
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

A = zeros(nwd,nws,nT,nrows,ncols,'single');
B = zeros(nwd,nws,nT,nrows,ncols,'single');
D = zeros(nwd,nws,nT,nrows,ncols,'single');

A_r = zeros(nwd,nws,nT,nrows,ncols,'single');
B_r = zeros(nwd,nws,nT,nrows,ncols,'single');
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
               
               Sum_Above = zeros(nrows,ncols,'single');
               Sum_Below = zeros(nrows,ncols,'single');
               Sum_SG = zeros(nrows,ncols,'single');

               Sum_Abover = zeros(nrows,ncols,'single');
               Sum_Belowr = zeros(nrows,ncols,'single');
               Sum_SGr = zeros(nrows,ncols,'single');
               
               parfor il2 = 1:nl2
                   x_c = x_c_inp(il2);
                   y_c = y_c_inp(il2);
                   
                   %ws = ws_inp(il2);
                   wd = wd_inp(il2);
                   
                   xy_rot = [cos(wd) sin(wd);-sin(wd) cos(wd)]*[x_c;y_c];
                   x_cr = xy_rot(1);
                   y_cr = xy_rot(2);
                   
                   ifov = xtrack_inp(il2);
                   vcd = vcd_inp(il2);
                   vcd_unc = vcde_inp(il2);
                   
                   u = u_km(ifov);
                   v = v_km(ifov);
                   t = t_km(ifov);
                   
                   pixel_edge = max(abs([pixel_left(ifov),pixel_down(ifov)]));
                   pixel_edge_inflate = max([xmargin ymargin]);
                   
                   local_left = x_c-pixel_edge_inflate*pixel_edge;
                   local_right = x_c+pixel_edge_inflate*pixel_edge;
                   
                   local_bottom = y_c-pixel_edge_inflate*pixel_edge;
                   local_top = y_c+pixel_edge_inflate*pixel_edge;
                   
                   x_index = xgrid >= local_left & xgrid <= local_right;
                   y_index = ygrid >= local_bottom & ygrid <= local_top;
                   
                   x_local_mesh = xmesh(y_index,x_index);
                   y_local_mesh = ymesh(y_index,x_index);
                   
                   SG = F_2D_SG(x_local_mesh,y_local_mesh,x_c,y_c,2*v,2*u,2,2,-t);
                   
                   sum_above_local = zeros(nrows,ncols,'single');
                   sum_below_local = zeros(nrows,ncols,'single');
                   D_local = zeros(nrows,ncols,'single');
                   
                   D_local(y_index,x_index) = SG;
                   sum_above_local(y_index,x_index) = SG/(u*v)/vcd_unc*vcd;
                   sum_below_local(y_index,x_index) = SG/(u*v)/vcd_unc;
                   Sum_Above = Sum_Above + sum_above_local;
                   Sum_Below = Sum_Below + sum_below_local;
                   Sum_SG = Sum_SG+D_local;
                   
                   % repeat it with rotated x y
                   local_left = x_cr-pixel_edge_inflate*pixel_edge;
                   local_right = x_cr+pixel_edge_inflate*pixel_edge;
                   
                   local_bottom = y_cr-pixel_edge_inflate*pixel_edge;
                   local_top = y_cr+pixel_edge_inflate*pixel_edge;
                   
                   x_index = xgrid >= local_left & xgrid <= local_right;
                   y_index = ygrid >= local_bottom & ygrid <= local_top;
                   
                   x_local_mesh = xmesh(y_index,x_index);
                   y_local_mesh = ymesh(y_index,x_index);
                   
                   SG = F_2D_SG(x_local_mesh,y_local_mesh,x_cr,y_cr,2*v,2*u,2,2,-t+wd);
                   
                   sum_above_local = zeros(nrows,ncols,'single');
                   sum_below_local = zeros(nrows,ncols,'single');
                   D_local = zeros(nrows,ncols,'single');
                   
                   D_local(y_index,x_index) = SG;
                   sum_above_local(y_index,x_index) = SG/(u*v)/vcd_unc*vcd;
                   sum_below_local(y_index,x_index) = SG/(u*v)/vcd_unc;
                   Sum_Abover = Sum_Abover + sum_above_local;
                   Sum_Belowr = Sum_Belowr + sum_below_local;
                   Sum_SGr = Sum_SGr+D_local;
               end
               
               A(iwd,iws,iT,:,:) = Sum_Above;
               B(iwd,iws,iT,:,:) = Sum_Below;
               D(iwd,iws,iT,:,:) = Sum_SG;
               
               A_r(iwd,iws,iT,:,:) = Sum_Abover;
               B_r(iwd,iws,iT,:,:) = Sum_Belowr;
               D_r(iwd,iws,iT,:,:) = Sum_SGr;
            end
        end
    end
end

output_regrid.A = A;
output_regrid.B = B;
output_regrid.D = D;

output_regrid.A_r = A_r;
output_regrid.B_r = B_r;
output_regrid.D_r = D_r;

output_regrid.xgrid = xgrid;
output_regrid.ygrid = ygrid;
output_regrid.xmesh = single(xmesh);
output_regrid.ymesh = single(ymesh);

disp('Transfer x y coordinates back to lat lon...')
[tmplat,tmplon] = minvtran(outp_xy.mstruct,xmesh(:),ymesh(:));
output_regrid.latmesh = single(reshape(tmplat,size(xmesh)));
output_regrid.lonmesh = single(reshape(tmplon,size(xmesh)));

output_regrid.mstruct = outp_xy.mstruct;
output_regrid.max_lon = max_lon;
output_regrid.min_lon = min_lon;

output_regrid.max_lat = max_lat;
output_regrid.min_lat = min_lat;
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
