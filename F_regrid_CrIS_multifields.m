function output_regrid = F_regrid_CrIS_multifields(inp,output_subset)
% updated from F_regrid_IASI.m by Kang Sun on 2019/02/21 to accomadate
% CrIS/AMoN comparison


output_regrid = [];
Res = inp.Res;
MinLon = inp.MinLon;
MaxLon = inp.MaxLon;
MinLat = inp.MinLat;
MaxLat = inp.MaxLat;
if isfield(inp,'MarginLon')
    MarginLon = inp.MarginLon;
else
    MarginLon = 0.5;
end
if isfield(inp,'MarginLat')
    MarginLat = inp.MarginLat;
else
    MarginLat = 0.5;
end
if ~isfield(inp,'k')
    k = 2;
else
    k = inp.k;
end

if ~isfield(inp,'min_dofs')
    min_dofs = 0;
else
    min_dofs = inp.min_dofs;
end

if isfield(inp,'vcdname')
    output_subset.colnh3_true = output_subset.colnh3;
    output_subset.colnh3 = output_subset.(inp.vcdname);
end
Startdate = inp.Startdate;
Enddate = inp.Enddate;

% define x y grids
xgrid = (MinLon+0.5*Res):Res:MaxLon;
ygrid = (MinLat+0.5*Res):Res:MaxLat;
nrows = length(ygrid);
ncols = length(xgrid);

% define x y mesh
[xmesh, ymesh] = meshgrid(single(xgrid),single(ygrid));

% construct a rectangle envelopes the orginal pixel
xmargin = 3; % how many times to extend zonally
ymargin = 2; % how many times to extend meridonally

% how many points to define an IASI ellipse?
npoint_ellipse = 10;


if length(Startdate) == 3
    start_utc = datenum([Startdate 0 0 0]);
    end_utc = datenum([Enddate 23 59 59]);
else
    start_utc = datenum(Startdate);
    end_utc = datenum(Enddate);
end
f1 = output_subset.utc >= start_utc ...
    & output_subset.utc <= end_utc;

f2 = output_subset.lat >= MinLat-MarginLat & output_subset.lat <= MaxLat+MarginLat...
    & output_subset.lon >= MinLon-MarginLon & output_subset.lon <= MaxLon+MarginLon;
f3 = ~isnan(output_subset.colnh3);
f4 = ~isnan(output_subset.colnh3error);
f5 = output_subset.dofall >= min_dofs;
disp(['Regriding pixels from ',datestr(start_utc),' to ',...
    datestr(end_utc)])
disp(['In time & in box: ',num2str(sum(f1&f2))...
    '; NH3 not NaN: ',num2str(sum(f1&f2&f3)),...
    '; NH3 error not NaN: ',num2str(sum(f1&f2&f3&f4))]);
validmask = f1&f2&f3&f4&f5;
nL2 = sum(validmask);disp([num2str(nL2),' pixels to be regridded...'])
if nL2 == 0
    warning('No valid L2 pixels were found!')
    return;
end
Lon = output_subset.lon(validmask);
Lat = output_subset.lat(validmask);
U = output_subset.u(validmask);
V = output_subset.v(validmask);
T = output_subset.t(validmask);
ColNH3 = output_subset.colnh3(validmask);
ColNH3e = output_subset.colnh3error(validmask);
sfcvmr = output_subset.sfcvmrall(validmask)/1e6;
dofs = output_subset.dofall(validmask);

Tsurf = output_subset.T_sfc_NARR(validmask);
psurf = output_subset.psurf(validmask);
pblp = output_subset.pblp(validmask);
ugperm3_a = sfcvmr.*psurf*0.017031/8.31446./Tsurf*1e9;
Gamma = 1;
ugperm3_r = ColNH3*1e4*9.8*0.02896/6.022140857e23./Gamma./pblp ...
    .*psurf*0.017031/8.31446./Tsurf*1e9;

ColNH3 = log10(ColNH3);
ColNH3e = log10(ColNH3e);
ugperm3_a = log10(ugperm3_a);
ugperm3_r = log10(ugperm3_r);
ofield = {'ColNH3','Tsurf','ugperm3_a','ugperm3_r','dofs'};
Sum_Above = [];
for ifield = 1:length(ofield)
    Sum_Above.(ofield{ifield}) = zeros(nrows,ncols);
end
Sum_Below = zeros(nrows,ncols);
D = zeros(nrows,ncols,'single');

count = 1;
for i = 1:nL2
    lon = Lon(i);
    lat = Lat(i);
    u = U(i);
    v = V(i);
    t = T(i);
    %     colnh3 = ColNH3(i);
    colnh3e = ColNH3e(i);
    % minlon_e is the minimum lon of the elliptical pixel, does not have to
    % be super accurate; minlat_e is the minmum lat; X is the polygon
    [~, minlon_e, minlat_e] =...
        F_construct_ellipse([lon;lat],v,u,t,npoint_ellipse,0);
    
    local_left = lon-xmargin*(lon-minlon_e);
    local_right = lon+xmargin*(lon-minlon_e);
    
    local_bottom = lat-ymargin*(lat-minlat_e);
    local_top = lat+ymargin*(lat-minlat_e);
    
    x_local_index = xgrid >= local_left & xgrid <= local_right;
    y_local_index = ygrid >= local_bottom & ygrid <= local_top;
    
    x_local_mesh = xmesh(y_local_index,x_local_index);
    y_local_mesh = ymesh(y_local_index,x_local_index);
    SG = F_2D_SG_rotate(x_local_mesh,y_local_mesh,lon,lat,2*v,2*u,k,-t);
    SG(isnan(SG)|isinf(SG)) = 0;
    for ifield = 1:length(ofield)
        tmpvar = eval([ofield{ifield},'(i)']);
        Sum_Above.(ofield{ifield})(y_local_index,x_local_index) ...
            = Sum_Above.(ofield{ifield})(y_local_index,x_local_index)+...
            SG/(v*u)/colnh3e*tmpvar;
    end
    Sum_Below(y_local_index,x_local_index) = Sum_Below(y_local_index,x_local_index)+...
        SG/(v*u)/colnh3e;
    D(y_local_index,x_local_index) = D(y_local_index,x_local_index)+SG;
    if i == count*round(nL2/10)
        disp([num2str(count*10),' % finished'])
        count = count+1;
    end
end
output_regrid.A = Sum_Above;
output_regrid.B = Sum_Below;
% output_regrid.C = Sum_Above./Sum_Below;
output_regrid.D = D;

output_regrid.xgrid = xgrid;
output_regrid.ygrid = ygrid;
