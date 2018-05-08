% script to read in L2g file and interpolate weather data at each L2 pixel
% location. This way no need to always interpolate during regridding.
% Written by Kang Sun on 2018/05/07

addpath('C:\Users\Kang Sun\Documents\GitHub\PU_KS_share\')
NARR_download_dir = 'C:\data_ks\NARR\';
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
%%
clc
L2gdir = 'C:\data_ks\IASIa\L2g\';
for iyear = 2008
    if exist([L2gdir,'CONUS_v2_',num2str(iyear),'.mat'],'file')
        l2gfn = [L2gdir,'CONUS_v2_',num2str(iyear),'.mat'];
        load(l2gfn);
    else
        l2gfn = [L2gdir,'CONUS_',num2str(iyear),'.mat'];
        load(l2gfn);
    end
    if isfield(output_subset,'lonc')
        Lat = output_subset.latc;Lon = output_subset.lonc;
    else
        Lat = output_subset.lat;Lon = output_subset.lon;
    end
    
    [UTC, I] = sort(output_subset.utc);
    Lat = Lat(I);
    Lon = Lon(I);
    % allocate space for wind, associated with each satellite pixel
    U_sat_pbl = nan(size(Lat),'single');
    V_sat_pbl = U_sat_pbl;
    T_sat_pbl = U_sat_pbl;
    
    day_1 = datevec(UTC(1));
    day_end = datevec(UTC(end));
    
    nL2 = length(Lat);
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
                    
                    inp_interp_narr.mstruct_narr = mstruct_narr;
                    
                    inp_interp_narr.P_pblmax = 150;% maximum pbl thickness in hPa
                    outp_interp_narr = F_interp_narr_v2(inp_interp_narr);
                    
                    idx_interp_hour = idx_interp_day(bin == ibin);
                    U_sat_pbl(idx_interp_hour) = outp_interp_narr.u_sat_pbl;
                    V_sat_pbl(idx_interp_hour) = outp_interp_narr.v_sat_pbl;
                    T_sat_pbl(idx_interp_hour) = outp_interp_narr.T_sat_pbl;
                end
            end
        end
    end
    fn = fieldnames(output_subset);
    for ifn = 1:length(fn)
        output_subset.(fn{ifn}) = output_subset.(fn{ifn})(I,:);
    end
    output_subset.U_sat_pbl = U_sat_pbl;
    output_subset.V_sat_pbl = V_sat_pbl;
    output_subset.T_sat_pbl = T_sat_pbl;
    save(l2gfn,'output_subset','inp_subset');
end
%% test regrid
clc
inp_regrid = [];
% Resolution of oversampled L3 data?
inp_regrid.Res = 0.02; % in degree
% Colorado
inp_regrid.MinLon = -128;
inp_regrid.MaxLon = -65;
inp_regrid.MinLat = 25;
inp_regrid.MaxLat = 50;

inp_regrid.vcdname = 'T_sat_pbl';
inp_regrid.Startdate = [2008 5 1];
inp_regrid.Enddate = [2008 9 1];
output_regrid = F_regrid_IASI(inp_regrid,output_subset);
%%
h = pcolor(output_regrid.xgrid,output_regrid.ygrid,output_regrid.C);
set(h,'edgecolor','none')