function output_regrid = F_regrid_IASI_rotate_wind(inp,output_subset)
% updated from F_regrid_IASI_NARR.m by Kang Sun on 2019/02/06 to work with
% l2g data with wind data already there

output_regrid = [];

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

if ~isfield(inp,'k')
    k = 2;
else
    k = inp.k;
end

nwd = length(wd_bin);
nws = length(ws_bin)-1;
nT = length(T_bin)-1;

if isfield(inp,'which_wind')
which_wind = inp.which_wind;
else
    which_wind = 'U_pbl_NARR';
end
if ~isfield(output_subset,which_wind)
    error(['Wind source ',which_wind,' does not exist in your L2g data!']);
end

% % the pixel geometry is much more elegant than the lat lon version
% u_km = inp.u_km;
% v_km = inp.v_km;
% %t_km = inp.t_km;% t_km is not a constant field. need case-by-case calc
% pixel_left = inp.pixel_left;
% pixel_down = inp.pixel_down;

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
xmargin = 2; % how many times to extend zonally
ymargin = 2; % how many times to extend meridonally

f1 = output_subset.utc >= single(datenum([Startdate 0 0 0])) ...
    & output_subset.utc <= single(datenum([Enddate 23 59 59]));
f2 = output_subset.lat >= min_lat & output_subset.lat <= max_lat...
    & output_subset.lon >= min_lon & output_subset.lon <= max_lon;
f3 = ~isnan(output_subset.colnh3);
f4 = ~isnan(output_subset.colnh3error);

validmask = f1&f2&f3&f4;

nL2 = sum(validmask);
disp([num2str(nL2),' pixels to be regridded...'])

Lon = output_subset.lon(validmask);
Lat = output_subset.lat(validmask);
% Ifov = output_subset.ifov(validmask);
ColNH3 = output_subset.colnh3(validmask);
ColNH3e = output_subset.colnh3error(validmask);
UTC = output_subset.utc(validmask);

% "e" stands for "ellipse"
Ue = output_subset.u(validmask);
Ve = output_subset.v(validmask);
Te = output_subset.t(validmask);

u_vec_rot = output_subset.(which_wind)(validmask);
v_vec_rot = output_subset.(['V',which_wind(2:end)])(validmask);
T_vec_rot = output_subset.(['T',which_wind(2:end)])(validmask);

disp('Converting pixel lat lon to x y in km...')
% calculate four anchor points for each oval
Lon_r = nan(nL2,4,'single');
Lat_r = nan(nL2,4,'single');
for il2 = 1:nL2
    X = F_construct_ellipse([Lon(il2);Lat(il2)],Ve(il2),Ue(il2),Te(il2),5,0); 
    Lon_r(il2,:) = X(1,1:end-1);
    Lat_r(il2,:) = X(2,1:end-1);
end
% call function F_latlon2xy to convert lat lon to x y.
inp_xy = [];
inp_xy.clon = clon;
inp_xy.clat = clat;
inp_xy.lon = [Lon_r,Lon(:)];
inp_xy.lat = [Lat_r,Lat(:)];
outp_xy = F_latlon2xy(inp_xy);

disp('Categorize data based on WS/WD/T...')

ws_vec_rot = sqrt(u_vec_rot.^2+v_vec_rot.^2);
wd_vec_rot = u_vec_rot;
wd_vec_rot(v_vec_rot >= 0) = acos(u_vec_rot(v_vec_rot >= 0)./ws_vec_rot(v_vec_rot >= 0));
wd_vec_rot(v_vec_rot < 0) = 2*pi-acos(u_vec_rot(v_vec_rot < 0)./ws_vec_rot(v_vec_rot < 0));

[~,~,ws_idx] = histcounts(ws_vec_rot,ws_bin);
[~,~,wd_idx] = histcounts(wd_vec_rot,[0 wd_bin 2*pi]);
wd_idx(wd_idx == nwd+1) = 1;
[~,~,T_idx] = histcounts(T_vec_rot,T_bin);

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
               x_r_inp = outp_xy.x(use_idx,1:4);
               y_r_inp = outp_xy.y(use_idx,1:4);
               x_c_inp = outp_xy.x(use_idx,5);
               y_c_inp = outp_xy.y(use_idx,5);
               
               %ws_inp = ws_vec_rot(use_idx);
               wd_inp = wd_vec_rot(use_idx);
               
               vcd_inp = ColNH3(use_idx);
               vcde_inp = ColNH3e(use_idx);
%                xtrack_inp = Ifov(use_idx);
               
               Sum_Above = zeros(nrows,ncols,'single');
               Sum_Below = zeros(nrows,ncols,'single');
               Sum_SG = zeros(nrows,ncols,'single');

               Sum_Abover = zeros(nrows,ncols,'single');
               Sum_Belowr = zeros(nrows,ncols,'single');
               Sum_SGr = zeros(nrows,ncols,'single');
               
               parfor il2 = 1:nl2
                   x_r = x_r_inp(il2,:);
                   y_r = y_r_inp(il2,:);
                   x_c = x_c_inp(il2);
                   y_c = y_c_inp(il2);
                   
%                    ws = ws_inp(il2);
                   wd = wd_inp(il2);
                   
                   xy_rot = [cos(wd) sin(wd);-sin(wd) cos(wd)]*[x_r x_c;y_r y_c];
%                    x_rr = xy_rot(1,1:4);
%                    y_rr = xy_rot(2,1:4);
                   x_cr = xy_rot(1,5);
                   y_cr = xy_rot(2,5);
                   
%                    ifov = xtrack_inp(il2);
                   vcd = vcd_inp(il2);
                   vcd_unc = vcde_inp(il2);
                   
                   dy = y_r(1)-y_r(3);
                   dx = x_r(1)-x_r(3);
                   if dx >= 0
                       t_km_local = atan(dy/dx);
                   else
                       t_km_local = pi+atan(dy/dx);
                   end
                   v_km_local = sqrt(dx^2+dy^2)/2;
                   u_km_local = sqrt((y_r(4)-y_r(2))^2+(x_r(4)-x_r(2))^2)/2;
                   
                   pixel_edge = max([v_km_local,u_km_local]);
                   pixel_edge_inflate = max([xmargin ymargin]);
                   
                   local_left = x_c-pixel_edge_inflate*pixel_edge;
                   local_right = x_c+pixel_edge_inflate*pixel_edge;
                   
                   local_bottom = y_c-pixel_edge_inflate*pixel_edge;
                   local_top = y_c+pixel_edge_inflate*pixel_edge;
                   
                   x_index = xgrid >= local_left & xgrid <= local_right;
                   y_index = ygrid >= local_bottom & ygrid <= local_top;
                   
                   x_local_mesh = xmesh(y_index,x_index);
                   y_local_mesh = ymesh(y_index,x_index);
                   
                   SG = F_2D_SG_rotate(x_local_mesh,y_local_mesh,x_c,y_c,2*v_km_local,2*u_km_local,k,-t_km_local);
                   
%                    close
%                    hold on
%                    h = pcolor(x_local_mesh,y_local_mesh,double(SG));
%                    set(h,'edgecolor','none')
%                    plot(x_r,y_r,'*',x_c,y_c,'o',x_rr,y_rr,'*',x_cr,y_cr,'o')
%                    quiver(x_c,y_c,ws*cos(wd),ws*sin(wd))
                   
                   sum_above_local = zeros(nrows,ncols,'single');
                   sum_below_local = zeros(nrows,ncols,'single');
                   D_local = zeros(nrows,ncols,'single');
                   
                   D_local(y_index,x_index) = SG;
                   sum_above_local(y_index,x_index) = SG/(u_km_local*v_km_local)/vcd_unc*vcd;
                  sum_below_local(y_index,x_index) = SG/(u_km_local*v_km_local)/vcd_unc;
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
                   
                   SG = F_2D_SG_rotate(x_local_mesh,y_local_mesh,x_cr,y_cr,2*v_km_local,2*u_km_local,k,-t_km_local+wd);
                   
%                    close
%                    hold on
%                    h = pcolor(x_local_mesh,y_local_mesh,double(SG));
%                    set(h,'edgecolor','none')
%                    plot(x_r,y_r,'*',x_c,y_c,'o',x_rr,y_rr,'*',x_cr,y_cr,'o')
%                    quiver(x_c,y_c,ws*cos(wd),ws*sin(wd))
                   
                   sum_above_local = zeros(nrows,ncols,'single');
                   sum_below_local = zeros(nrows,ncols,'single');
                   D_local = zeros(nrows,ncols,'single');
                   
                   D_local(y_index,x_index) = SG;
                   sum_above_local(y_index,x_index) = SG/(u_km_local*v_km_local)/vcd_unc*vcd;
                   sum_below_local(y_index,x_index) = SG/(u_km_local*v_km_local)/vcd_unc;
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
