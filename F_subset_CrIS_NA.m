function output = F_subset_CrIS_NA(inp)
% Matlab function to subset CrIS data from Shephard and Cady-Pereira.
% Saved from F_subset_CrIS.m on 2019/01/07 for Karen's AER data, and
% save different a priori shape information

olddir = pwd;
Startdate = inp.Startdate;
Enddate = inp.Enddate;
MinLat = inp.MinLat;
MinLon = inp.MinLon;
MaxLat = inp.MaxLat;
MaxLon = inp.MaxLon;

% grey area, L2 pixel center cannot be there, but pixel boundaries can
MarginLat = 0.5;
MarginLon = 0.5;

L2dir = inp.L2dir;

if ~isfield(inp,'errormodel')
    errormodel = 'total';
else
    errormodel = inp.errormodel;
end

day_array = (datenum(Startdate):1:datenum(Enddate))';
datevec_array = datevec(day_array);
year_array = datevec_array(:,1);
month_array = datevec_array(:,2);
date_array = datevec_array(:,3);
nday = length(day_array);
varname = {'DOF','Day_Night_Flag','LandFraction','Latitude','Longitude'...
    'Quality_Flag','Run_ID','mdate','rvmr','rvmr_error','tot_col','xretv',...
    'total_covariance_error','noise_error_covariance','pressure','xa','avg_kernel','xa_Type'};
latall = single([]);
lonall = latall;
ifovall = latall;
colnh3all = latall;
utcall = []; % time needs to be double
totErrall = latall;
sfcvmrall = latall;
dofall = latall;
xa_type_all = cell(0);
for iday = 1:nday
    day_dir = [L2dir,num2str(year_array(iday)),'/',...
        num2str(month_array(iday),'%02d'),'/',num2str(date_array(iday),'%02d'),'/'];
    cd(day_dir)
    flist = dir([day_dir,'*.nc']);
    nfile = length(flist);
    for ifile = 1:nfile
        fn = flist(ifile).name;
        datavar = [];
        try
            % load netcdf data as hdf-5
            file_id = H5F.open (fn, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
            for i = 1:length(varname)
                % Open the dataset.
                DATAFIELD_NAME = ['/',varname{i}];
                data_id = H5D.open (file_id, DATAFIELD_NAME);
                datavar.(varname{i}).data=H5D.read (data_id);
            end
            H5D.close (data_id);
            H5F.close (file_id);
        catch ME
            warning(['Error occurs for ',datestr(day_array(iday)),', orbit ',...
                num2str(ifile),', the message is "',ME.message,'"']);
            continue;
        end
        % extract useful information
        f1 = datavar.Latitude.data >= MinLat+MarginLat & datavar.Latitude.data <= MaxLat-MarginLat & ...
            datavar.Longitude.data >= MinLon+MarginLon & datavar.Longitude.data <= MaxLon-MarginLon;
        f2 = datavar.Quality_Flag.data == 3 | datavar.Quality_Flag.data == 4;
        f3 = datavar.Day_Night_Flag.data == 1;
        f4 = ~(datavar.tot_col.data > -999.501 & datavar.tot_col.data < -999.499);
        validmask = f1 & f2 & f3 & f4;
        
        if sum(validmask) > 0
            disp(['You have ',sprintf('%5d',sum(validmask(:))),...
                ' valid L2 pixels on ',datestr(day_array(iday)),', orbit ',...
                num2str(ifile)]);
            utcall = cat(1,utcall,datavar.mdate.data(validmask)+366);
            
            latall = cat(1,latall,single(datavar.Latitude.data(validmask)));
            lonall = cat(1,lonall,single(datavar.Longitude.data(validmask)));
            tmprunID = datavar.Run_ID.data(validmask);
            tmpfov = nan(length(tmprunID),1);
            tmpfor = tmpfov;
            for irun = 1:length(tmprunID)
                tmpfov(irun) = str2double(tmprunID{irun}(end-2:end));
                tmpfor(irun) = str2double(tmprunID{irun}(end-7:end-4));
            end
            ifovall = cat(1,ifovall,single((tmpfor-1)*9+tmpfov));
            colnh3all = cat(1,colnh3all,single(datavar.tot_col.data(validmask)));
            
            dofall = cat(1,dofall,single(datavar.DOF.data(validmask)));
            xa_type_all = cat(1,xa_type_all,datavar.xa_Type.data(validmask));
            
            switch errormodel
                case 'rvmr'
                    error0 = single(datavar.rvmr_error.data(1,validmask)'./datavar.rvmr.data(1,validmask)'.*datavar.tot_col.data(validmask));
                    sfcvmr0 = single(datavar.xretv.data(2,validmask))';
                case 'total'
                    nobs = sum(validmask);
                    pressure0 = datavar.pressure.data(:,validmask);
                    xretv0 = datavar.xretv.data(:,validmask);
                    noise_error0 = datavar.noise_error_covariance.data(:,:,validmask);
                    total_error0 = datavar.total_covariance_error.data(:,:,validmask);
                    ak0 = datavar.avg_kernel.data(:,:,validmask);
                    % xa0 = datavar.xa.data(:,validmask);
                    ak_colm = 0*xretv0;
                    tot_col_test = nan(nobs,1);
                    sfcvmr0 = nan(nobs,1);
                    noise_error_colm = datavar.Latitude.data(validmask)*0;
                    total_error_colm = noise_error_colm;
                    % loop over observations
                    for io = 1:nobs
                        index = pressure0(:,io) > 0;
                        pressure = pressure0(index,io);
                        nlev = length(pressure);
                        dp = zeros(nlev,1);
                        dp(1) = (pressure(1)-pressure(2))/2;
                        for ip = 2:nlev-1
                            dp(ip) = (pressure(ip-1)-pressure(ip))/2+(pressure(ip)-pressure(ip+1))/2;
                        end
                        dp(nlev) = pressure(nlev-1)-pressure(nlev);
                        trans = 2.12e16*dp;
                        % calculate column AK
                        xretv = xretv0(index,io);
                        %     xa = xa0(index,io);
                        ak = squeeze(ak0(index,index,io));
                        noise_error = squeeze(noise_error0(index,index,io));
                        total_error = squeeze(total_error0(index,index,io));
                        ak_colm(index,io) = (trans.*xretv)'*ak;
                        % calculate errors
                        xarr = diag(1e-6*ones(nlev,1).*xretv);
                        trans_error = trans*1e6;
                        sx = xarr*noise_error*xarr;
                        noise_error_colm(io) = sqrt(trans_error'*sx*trans_error);
                        st = xarr*total_error*xarr;
                        total_error_colm(io) = sqrt(trans_error'*st*trans_error);
                        tot_col_test(io) = sum(trans.*xretv);
                        sfcvmr0(io) = xretv(1);
                    end
                    error0 = total_error_colm;
                case 'noise'
                    nobs = sum(validmask);
                    pressure0 = datavar.pressure.data(:,validmask);
                    xretv0 = datavar.xretv.data(:,validmask);
                    noise_error0 = datavar.noise_error_covariance.data(:,:,validmask);
                    total_error0 = datavar.total_covariance_error.data(:,:,validmask);
                    ak0 = datavar.avg_kernel.data(:,:,validmask);
                    % xa0 = datavar.xa.data(:,validmask);
                    ak_colm = 0*xretv0;
                    tot_col_test = nan(nobs,1);
                    sfcvmr0 = nan(nobs,1);
                    noise_error_colm = datavar.Latitude.data(validmask)*0;
                    total_error_colm = noise_error_colm;
                    % loop over observations
                    for io = 1:nobs
                        index = pressure0(:,io) > 0;
                        pressure = pressure0(index,io);
                        nlev = length(pressure);
                        dp = zeros(nlev,1);
                        dp(1) = (pressure(1)-pressure(2))/2;
                        for ip = 2:nlev-1
                            dp(ip) = (pressure(ip-1)-pressure(ip))/2+(pressure(ip)-pressure(ip+1))/2;
                        end
                        dp(nlev) = pressure(nlev-1)-pressure(nlev);
                        trans = 2.12e16*dp;
                        % calculate column AK
                        xretv = xretv0(index,io);
                        %     xa = xa0(index,io);
                        ak = squeeze(ak0(index,index,io));
                        noise_error = squeeze(noise_error0(index,index,io));
                        total_error = squeeze(total_error0(index,index,io));
                        ak_colm(index,io) = (trans.*xretv)'*ak;
                        % calculate errors
                        xarr = diag(1e-6*ones(nlev,1).*xretv);
                        trans_error = trans*1e6;
                        sx = xarr*noise_error*xarr;
                        noise_error_colm(io) = sqrt(trans_error'*sx*trans_error);
                        st = xarr*total_error*xarr;
                        total_error_colm(io) = sqrt(trans_error'*st*trans_error);
                        tot_col_test(io) = sum(trans.*xretv);
                        sfcvmr0(io) = xretv(1);
                    end
                    error0 = noise_error_colm;
            end
            
            totErrall = cat(1,totErrall,...
                error0);
            sfcvmrall = cat(1,sfcvmrall,sfcvmr0);
            
        end
    end
end

output.colnh3 = colnh3all;
output.colnh3error = totErrall;
output.lat = latall;
output.lon = lonall;
output.ifov = ifovall;
output.utc = utcall;
output.sfcvmrall = sfcvmrall;
output.dofall = dofall;
output.xa_type_all = xa_type_all;
cd(olddir);