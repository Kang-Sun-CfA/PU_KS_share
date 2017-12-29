function output = F_subset_IASI_v2(inp)
% extract useful IASI data from L2 files in L2dir. Designed for the v2.1
% data downloaded from http://cds-espri.ipsl.fr/etherTypo/index.php?id=1700&L=1

% Re-written from F_subset_IASI.m by Kang Sun on 2017/12/28

Startdate = inp.Startdate;
Enddate = inp.Enddate;
MinLat = inp.MinLat;
MinLon = inp.MinLon;
MaxLat = inp.MaxLat;
MaxLon = inp.MaxLon;
MaxCF = inp.MaxCF*100;

% grey area, L2 pixel center cannot be there, but pixel boundaries can
MarginLat = 0.5;
MarginLon = 0.5;

L2dir = inp.L2dir;
if ~isfield(inp,'a_or_b')
    a_or_b = 'a';
else
    a_or_b = inp.a_or_b;
end

if isfield(inp,'update_ifov')
    update_ifov = inp.update_ifov;
else
    update_ifov = true;
end

day_array = (datenum(Startdate):1:datenum(Enddate))';
datevec_array = datevec(day_array);
year_array = datevec_array(:,1);
month_array = datevec_array(:,2);
date_array = datevec_array(:,3);
nday = length(day_array);

latall = single([]);
lonall = latall;
ifovall = latall;
colnh3all = latall;
utcall = []; % time needs to be double
totErrall = latall;
cfall = latall;

for iday = 1:nday
    ystr = num2str(year_array(iday));
    mstr = num2str(month_array(iday),'%02d');
    dstr = num2str(date_array(iday),'%02d');
    fn = [L2dir,ystr,'/',num2str(month_array(iday)),'/',...
        'IASI_metop',a_or_b,'_L2_NH3_',ystr,mstr,dstr,'_V2.1_AM_GLOBAL.nc'];
    ifovfn = [L2dir,ystr,'/',num2str(month_array(iday)),'/',...
        'IASI_metop',a_or_b,'_L2_NH3_',ystr,mstr,dstr,'_ifov.mat'];
    ifov_save = [];
    ifov_save.ifovfn = ifovfn;
    ifov_save.update = update_ifov;
    if ~exist(fn,'file')
        disp([datestr(day_array(iday)),' L2 file is missing!'])
        continue
    end
    data = F_manipulate_iasi_nc(fn,ifov_save);
    f1 = data.lat >= MinLat+MarginLat & data.lat <= MaxLat-MarginLat & ...
        data.lon >= MinLon+MarginLon & data.lon <= MaxLon-MarginLon;
    f2 = ~isnan(data.colnh3) & ~isnan(data.error);
    f3 = data.CLcov <= MaxCF;
    validmask = f1 & f2 & f3;
    
    if sum(validmask) > 0
        disp(['You have ',sprintf('%5d',sum(validmask(:))),' valid L2 pixels on ',datestr(day_array(iday))]);
        utcall = cat(1,utcall,data.fractional_day(validmask)+day_array(iday));
        
        latall = cat(1,latall,single(data.lat(validmask)));
        lonall = cat(1,lonall,single(data.lon(validmask)));
        ifovall = cat(1,ifovall,single(data.ifov(validmask)));
        colnh3all = cat(1,colnh3all,single(data.colnh3(validmask)));
        totErrall = cat(1,totErrall,single(data.error(validmask)));
        cfall = cat(1,cfall,single(data.CLcov(validmask)));
    end
    
end

% NH3 column error, absolute
colnh3errorall = abs(colnh3all.*totErrall/100);

output.colnh3 = colnh3all;
output.colnh3error = colnh3errorall;
output.lat = latall;
output.lon = lonall;
output.ifov = ifovall;
output.utc = utcall;
output.cloudfrac = cfall;