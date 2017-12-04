function output_tessellate = F_tessellate_IASI_km(inp,output_subset)
% function to take in the output from F_subset_IASI.m and regrid these L2
% data to a L3 grid, centered at clon and clat and bounded by max_x and 
% max_y with resolution res in km. 

% Use tessellation method, originated from Lei Zhu and Kai Yang

% Modified by Kang Sun on 2017/11/29

output_tessellate = [];
% IASI pixel area in km
A_IASI = [176.112319946289,175.997787475586,197.295578002930,197.296386718750,134.232391357422,134.307479858398,147.597915649414,147.511016845703,106.910308837891,106.910293579102,115.835662841797,115.835952758789,88.1102142333984,88.1099700927734,94.3334960937500,94.3336944580078,74.6136627197266,74.6135177612305,79.0981597900391,79.0983047485352,64.6809158325195,64.6807861328125,68.0252914428711,68.0253448486328,57.2129020690918,57.1962089538574,59.7423744201660,59.7241744995117,51.4984931945801,51.5112380981445,53.4395408630371,53.4396286010742,47.1394691467285,47.1294631958008,48.6127166748047,48.6127738952637,43.7793922424316,43.7793731689453,44.9149169921875,44.9149055480957,41.2354736328125,41.2354354858398,42.0903854370117,42.0904350280762,39.3517951965332,39.3517684936523,39.9780311584473,39.9780654907227,38.0385055541992,38.0384712219238,38.4638633728027,38.4605789184570,37.2202949523926,37.2202720642090,37.4690208435059,37.4670410156250,36.8619995117188,36.8632621765137,36.9433822631836,36.9425735473633,36.9458961486816,36.9475517272949,36.8649864196777,36.8641242980957,37.4771308898926,37.4791526794434,37.2297134399414,37.2281494140625,38.4804687500000,38.4804649353027,38.0526008605957,38.0526084899902,40.0019073486328,40.0066566467285,39.3770217895508,39.3770751953125,42.1226539611816,42.1291275024414,41.2644309997559,41.2644348144531,44.9572257995606,44.9657287597656,43.8254928588867,43.8177566528320,48.6675682067871,48.6785697937012,47.1893501281738,47.1893539428711,53.5103225708008,53.5103416442871,51.5756797790527,51.5885810852051,59.8340911865234,59.8340873718262,57.2961006164551,57.3128280639648,68.1210937500000,68.1691513061523,64.7891693115234,64.7891540527344,79.2890853881836,79.2891082763672,74.7851409912109,74.7851104736328,94.5494537353516,94.5494995117188,88.3023147583008,88.3407821655273,116.077011108398,116.198059082031,107.175537109375,107.228660583496,147.949081420898,148.125152587891,134.687866210938,134.687728881836,198.105438232422,198.105773925781,176.683227539063,176.683029174805];

Startdate = inp.Startdate;
Enddate = inp.Enddate;

res = inp.res;
max_x = inp.max_x;
max_y = inp.max_y;
clon = inp.clon;
clat = inp.clat;
npoints = inp.npoints;

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

% "e" stands for "ellipse"
Ue = output_subset.u(validmask);
Ve = output_subset.v(validmask);
Te = output_subset.t(validmask);

disp('Converting pixel lat lon to x y in km...')
% calculate points for each oval
Lon_r = nan(nL2,npoints,'single');
Lat_r = nan(nL2,npoints,'single');
for il2 = 1:nL2
    X = F_construct_ellipse([Lon(il2);Lat(il2)],Ve(il2),Ue(il2),Te(il2),npoints+1,0);
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

disp('Tessellating pixel by pixel...')
Sum_Above = zeros(nrows,ncols,'single');
Sum_Below = zeros(nrows,ncols,'single');
D = zeros(nrows,ncols,'single');

count = 1;

for iL2 = 1:nL2
    x = outp_xy.x(iL2,end);
    y = outp_xy.y(iL2,end);
    x_r = outp_xy.x(iL2,1:end-1);
    y_r = outp_xy.y(iL2,1:end-1);
    
    colnh3 = ColNH3(iL2);
    colnh3e = ColNH3e(iL2);
    
    pixel = [];
    pixel.nv = npoints;
    pixel.vList = [x_r(:)/res,y_r(:)/res];
    pixel.center = [x/res,y/res];
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
            row = max_x/res+floor(min(final_pixels(id_final).vList(1:final_pixels(id_final).nv,2))) + 1;
            col = max_y/res+floor(min(final_pixels(id_final).vList(1:final_pixels(id_final).nv,1))) + 1;
            
            if row >= 1 && row <= nrows && col >= 1 && col <= ncols
%             id_all = id_all + 1;
            % temp_area is the area of each sub polygon, at this stage
            [ifsquare,edges] = F_if_square(final_pixels(id_final));
            if ifsquare
                temp_area = edges(1)*edges(2)*res^2;
            else
                temp_area = F_polyarea(final_pixels(id_final))*res^2;
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
    Sum_Above = Sum_Above+Sub_Area/(A_IASI(Ifov(iL2))*pi)/colnh3e*colnh3;
    Sum_Below = Sum_Below+Sub_Area/(A_IASI(Ifov(iL2))*pi)/colnh3e;
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

disp('Transfer x y coordinates back to lat lon...')
[tmplat,tmplon] = minvtran(outp_xy.mstruct,xmesh(:),ymesh(:));
output_tessellate.latmesh = single(reshape(tmplat,size(xmesh)));
output_tessellate.lonmesh = single(reshape(tmplon,size(xmesh)));

output_tessellate.mstruct = outp_xy.mstruct;
output_tessellate.max_lon = max_lon;
output_tessellate.min_lon = min_lon;

output_tessellate.max_lat = max_lat;
output_tessellate.min_lat = min_lat;
