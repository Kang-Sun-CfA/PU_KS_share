function output_tessellate = F_tessellate_IASI(inp,output_subset)
% function to take in the output from F_subset_IASI.m and regrid these L2
% data to a L3 grid, centered at clon and clat and bounded by max_x and 
% max_y with resolution res in km. 

% Use tessellation method, originated from Lei Zhu and Kai Yang

% Modified by Kang Sun from F_tessellate_IASI_km.m on 2018/06/28

output_tessellate = [];

Startdate = inp.Startdate;
Enddate = inp.Enddate;

Res = inp.Res;
MinLon = inp.MinLon;
MaxLon = inp.MaxLon;
MinLat = inp.MinLat;
MaxLat = inp.MaxLat;

min_lat = MinLat;
max_lat = MaxLat;
min_lon = MinLon;
max_lon = MaxLon;

npoints = inp.npoints;


% define x y grids
xgrid = (MinLon+0.5*Res):Res:MaxLon;
ygrid = (MinLat+0.5*Res):Res:MaxLat;
nrows = length(ygrid);
ncols = length(xgrid);

% define x y mesh
[xmesh, ymesh] = meshgrid(single(xgrid),single(ygrid));
Lon_mesh = xmesh;
Lat_mesh = ymesh;

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

% "e" stands for "ellipse"
Ue = output_subset.u(validmask);
Ve = output_subset.v(validmask);
Te = output_subset.t(validmask);

% calculate points for each oval
Lon_r = nan(nL2,npoints,'single');
Lat_r = nan(nL2,npoints,'single');
for il2 = 1:nL2
    X = F_construct_ellipse([Lon(il2);Lat(il2)],Ve(il2),Ue(il2),Te(il2),npoints+1,0);
    Lon_r(il2,:) = X(1,1:end-1);
    Lat_r(il2,:) = X(2,1:end-1);
end

disp('Tessellating pixel by pixel...')
Sum_Above = zeros(nrows,ncols,'single');
Sum_Below = zeros(nrows,ncols,'single');
D = zeros(nrows,ncols,'single');

count = 1;

for iL2 = 1:nL2
    lat_r = Lat_r(iL2,:);
    lon_r = Lon_r(iL2,:);
    lat = Lat(iL2);
    lon = Lon(iL2);
    u = Ue(iL2);
    v = Ve(iL2);
    
    colnh3 = ColNH3(iL2);
    colnh3e = ColNH3e(iL2);
    
    pixel = [];
    pixel.nv = npoints;
    pixel.vList = [(lon_r(:)-MinLon)/Res,(lat_r(:)-MinLat)/Res];
    pixel.center = [(lon-MinLon)/Res,(lat-MinLat)/Res];
    % Used for the next step
%     id_all = 0;
    pixel_area = 0;
    
    % Perform Horizontal cut first at the integer grid lines
    [sub_pixels,n_sub_pixels] = F_HcakeCut( pixel );

    Sub_Area = zeros(nrows,ncols,'single');
    Pixels_count = Sub_Area;
    % Then perform Vertical cut for each sub pixel obtainted
    % from the Horizontal cut at the integer grid lines
    for id_sub = 1: n_sub_pixels
        [final_pixels, n_final_pixels] = F_VcakeCut( sub_pixels(id_sub) );
        for id_final = 1: n_final_pixels
            row = floor(min(final_pixels(id_final).vList(1:final_pixels(id_final).nv,2))) + 1;
            col = floor(min(final_pixels(id_final).vList(1:final_pixels(id_final).nv,1))) + 1;
            
            if row >= 1 && row <= nrows && col >= 1 && col <= ncols
%             id_all = id_all + 1;
            % temp_area is the area of each sub polygon, at this stage
            [ifsquare,edges] = F_if_square(final_pixels(id_final));
            if ifsquare
                temp_area = edges(1)*edges(2)*Res^2;
            else
                temp_area = F_polyarea(final_pixels(id_final))*Res^2;
            end
            
            pixel_area = pixel_area + temp_area;
            
            % Get the overlaped area between the pixel and each cell
            Sub_Area(row,col) = temp_area;
            
            Pixels_count(row,col) = Pixels_count(row,col) + 1;
            end
        end
    end
%     if abs((pixel_area-A_IASI(Ifov(iL2))*pi)/(A_IASI(Ifov(iL2))*pi)) > 0.1
%         error('Pixel area is more inaccurate than 10%. Increase your npoints!')
%     end
    % Sum weighted value and weights
    % Here use temp_area/A/VCD_Unc(p) as averaging weight, meaning that
    % averaging weight is assumed to be proportional to the ratio of the overlap area (temp_area) to the
    % pixel size (A) and inversely proportional to the error standard deviation (VCD_Unc(p)).
    % If you just want fraction of overlap area as averaging weight, use: temp_area/A
    % If you just want area weighted average, use: temp_area
    if ~isnan(colnh3) && ~isnan(colnh3e)
    Sum_Above = Sum_Above+Sub_Area/(u*v)/colnh3e*colnh3;
    Sum_Below = Sum_Below+Sub_Area/(u*v)/colnh3e;
    D = D+Pixels_count;
    end
    if iL2 == count*round(nL2/10)
        disp([num2str(count*10),' % finished'])
        count = count+1;
    end
end

output_tessellate.A = Sum_Above;
output_tessellate.B = Sum_Below;
output_tessellate.C = Sum_Above./Sum_Below;
output_tessellate.D = D;

output_tessellate.xgrid = xgrid;
output_tessellate.ygrid = ygrid;

output_tessellate.xmesh = single(xmesh);
output_tessellate.ymesh = single(ymesh);

output_tessellate.max_lon = max_lon;
output_tessellate.min_lon = min_lon;

output_tessellate.max_lat = max_lat;
output_tessellate.min_lat = min_lat;
