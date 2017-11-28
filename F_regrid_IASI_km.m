function output_regrid = F_regrid_IASI_km(inp,output_subset)
% function to take in the output from F_subset_IASI.m and regrid these L2
% data to a L3 grid, centered at clon and clat and bounded by max_x and 
% max_y with resolution res in km.

% Modified from F_regrid_IASI.m by Kang Sun on 2017/09/24

% Significant update on 2017/11/21 to correct rotation angle in km
% projection

output_regrid = [];

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
% Ifov = output_subset.ifov(validmask);
ColNH3 = output_subset.colnh3(validmask);
ColNH3e = output_subset.colnh3error(validmask);

% "e" stands for "ellipse"
Ue = output_subset.u(validmask);
Ve = output_subset.v(validmask);
Te = output_subset.t(validmask);

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

disp('Calculating spatial response functions pixel by pixel...')
Sum_Above = zeros(nrows,ncols,'single');
Sum_Below = zeros(nrows,ncols,'single');
D = zeros(nrows,ncols,'single');
count = 1;
for i = 1:nL2
    x = outp_xy.x(i,5);
    y = outp_xy.y(i,5);
    x_r = outp_xy.x(i,1:4);
    y_r = outp_xy.y(i,1:4);

%     ifov = Ifov(i);
    colnh3 = ColNH3(i);
    colnh3e = ColNH3e(i);
    
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
    
    local_left = x-pixel_edge_inflate*pixel_edge;
                   local_right = x+pixel_edge_inflate*pixel_edge;
                   
                   local_bottom = y-pixel_edge_inflate*pixel_edge;
                   local_top = y+pixel_edge_inflate*pixel_edge;
    
    x_local_index = xgrid >= local_left & xgrid <= local_right;
    y_local_index = ygrid >= local_bottom & ygrid <= local_top;
    
    x_local_mesh = xmesh(y_local_index,x_local_index);
    y_local_mesh = ymesh(y_local_index,x_local_index);
    SG = F_2D_SG(x_local_mesh,y_local_mesh,x,y,2*v_km_local,2*u_km_local,2,2,-t_km_local);
    
    Sum_Above(y_local_index,x_local_index) = Sum_Above(y_local_index,x_local_index)+...
        SG/(v_km_local*u_km_local)/colnh3e*colnh3;
    Sum_Below(y_local_index,x_local_index) = Sum_Below(y_local_index,x_local_index)+...
        SG/(v_km_local*u_km_local)/colnh3e;
    D(y_local_index,x_local_index) = D(y_local_index,x_local_index)+SG;
    
    if i == count*round(nL2/10)
        disp([num2str(count*10),' % finished'])
        count = count+1;
    end
end
output_regrid.A = Sum_Above;
output_regrid.B = Sum_Below;
output_regrid.C = Sum_Above./Sum_Below;
output_regrid.D = D;

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
