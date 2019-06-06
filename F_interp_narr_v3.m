function outp_interp_narr = F_interp_narr_v3(inp_interp_narr)
% interpolate NARR PBL weighted-mean data to lat lon and time
% written on 2019/02/04

NARR_download_dir = inp_interp_narr.NARR_download_dir;

nn = length(inp_interp_narr.hour);

% maximum pbl thickness in Pa
if ~isfield(inp_interp_narr,'P_pblmax')
    P_pblmax = 200*100;
else
    P_pblmax = inp_interp_narr.P_pblmax*100;% from hPa to Pa
end

lon_sat = inp_interp_narr.lon;
lat_sat = inp_interp_narr.lat;
mstruct_narr = inp_interp_narr.mstruct_narr;
[x_sat,y_sat] = mfwdtran(mstruct_narr,lat_sat,lon_sat);

scaleH = 7500; % scale height in km

tt = nan(nn,1);

for ihour = 1:nn
    narr_datenum = inp_interp_narr.datenum+inp_interp_narr.hour(ihour)/24;
    tt(ihour) = narr_datenum;
    tmp = datevec(narr_datenum);
    narr_year = num2str(tmp(1));
    narr_month = num2str(tmp(2),'%02d');
    narr_day = num2str(tmp(3),'%02d');
    narr_hour = num2str(tmp(4),'%02d');
    try
    narr_data_3d = load([NARR_download_dir,'/',narr_year,'/',narr_month,'/',...
        'subset_3D_',narr_day,'_',narr_hour,'.mat']);
    narr_data_sfc = load([NARR_download_dir,'/',narr_year,'/',narr_month,'/',...
        'subset_sfc_',narr_day,'_',narr_hour,'.mat']);
    catch s
        disp(s)
        continue
    end
    
    x = narr_data_3d.x;
    y = narr_data_3d.y;
    
    u = narr_data_3d.u;
    v = narr_data_3d.v;
    T = narr_data_3d.T;
    
    P_surf = narr_data_sfc.P_surf;
    T_surf = narr_data_sfc.T_surf;
    PBLH = narr_data_sfc.PBLH;
    if isfield(narr_data_3d,'P')
    P = narr_data_3d.P*100;
    else
        P = [550 600 650 700 725 750 775 800 825 850 875 900 925 950 975 1000]'*100;
    end
    if ~exist('nx','var')
        nx = length(x);
        ny = length(y);
        u_pbl = nan(ny,nx,nn,'single');
        u_sfc = nan(ny,nx,nn,'single');
        v_pbl = nan(ny,nx,nn,'single');
        v_sfc = nan(ny,nx,nn,'single');
        T_pbl = nan(ny,nx,nn,'single');
        T_sfc = nan(ny,nx,nn,'single');
    end
    P_pbl = P_surf.*exp(-PBLH/scaleH);
    P_pbl(P_pbl < P_surf-P_pblmax) = P_surf(P_pbl < P_surf-P_pblmax)-P_pblmax;
    for ix = 1:nx
        for iy = 1:ny
            int = P >= P_pbl(iy,ix) & P <= P_surf(iy,ix) & P >= P_surf(iy,ix)-P_pblmax;
            % sometimes the pbl is so thin...
            if sum(int) == 0
                int = find(P <= P_surf(iy,ix),1,'last');
            end
            localu = squeeze(u(int,iy,ix));
            localv = squeeze(v(int,iy,ix));
            localT = squeeze(T(int,iy,ix));
            u_pbl(iy,ix,ihour) = sum(localu.*P(int)/sum(P(int)));
            v_pbl(iy,ix,ihour) = sum(localv.*P(int)/sum(P(int)));
            
            u_sfc(iy,ix,ihour) = localu(end);
            v_sfc(iy,ix,ihour) = localv(end);
            T_sfc(iy,ix,ihour) = localT(end);
            
            T_pbl(iy,ix,ihour) = sum(squeeze(T(int,iy,ix)).*P(int)/sum(P(int)));
        end
    end
    
end
%%
outp_interp_narr.u_pbl = interpn(double(y),double(x),double(tt),double(u_pbl), ...
    double(y_sat),double(x_sat),inp_interp_narr.utc);
outp_interp_narr.v_pbl = interpn(double(y),double(x),double(tt),double(v_pbl), ...
    double(y_sat),double(x_sat),inp_interp_narr.utc);
outp_interp_narr.T_pbl = interpn(double(y),double(x),double(tt),double(T_pbl), ...
    double(y_sat),double(x_sat),inp_interp_narr.utc);
outp_interp_narr.u_sfc = interpn(double(y),double(x),double(tt),double(u_sfc), ...
    double(y_sat),double(x_sat),inp_interp_narr.utc);
outp_interp_narr.v_sfc = interpn(double(y),double(x),double(tt),double(v_sfc), ...
    double(y_sat),double(x_sat),inp_interp_narr.utc);
outp_interp_narr.T_sfc = interpn(double(y),double(x),double(tt),double(T_sfc), ...
    double(y_sat),double(x_sat),inp_interp_narr.utc);
%%
% close
% hold on
% quiver(inp_interp_narr.lon,inp_interp_narr.lat,outp_interp_narr.u_pbl,outp_interp_narr.v_pbl)
% quiver(inp_interp_narr.lon,inp_interp_narr.lat,outp_interp_narr.u_sfc,outp_interp_narr.v_sfc)
% 
% %%
% sfcH = interpn(narr_data_3d.P,narr_data_3d.y,narr_data_3d.x,narr_data_3d.H,...
%     narr_data_sfc.P_surf/100,ymesh,xmesh);
% sfcH(isnan(sfcH)) = 0;
% 
% pbl_topH = sfcH + narr_data_sfc.PBLH;
% h = pcolor(narr_data_sfc.x,narr_data_sfc.y,P_pbl);set(h,'edgecolor','none');colorbar