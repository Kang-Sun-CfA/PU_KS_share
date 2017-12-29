function data = F_manipulate_iasi_nc(fn,ifov_save)

% not really 'manipulate' because the iasi nc files are read-only. Have to
% calculate ifov and save separately.

% written by Kang Sun on 2017/12/28

% addpath('c:\Users\Kang Sun\Dropbox\matlab functions\')
% addpath('C:\Users\Kang Sun\Documents\GitHub\PU_KS_share\')
% fn = 'C:\data_ks\IASIa\L2\IASI_metopa_L2_NH3_20071001_V2.1_AM_GLOBAL.nc';
% % fn = 'C:\data_ks\IASIa\L2\IASI_metopa_L2_NH3_20160701_V2.1_AM_GLOBAL.nc';
% clc

ncid = netcdf.open(fn);

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
variable = [];
for i = 0:nvars-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,i);
    variable.(varname).data = netcdf.getVar(ncid,i);
    for k = 0:natts-1
        attname = netcdf.inqAttName(ncid,i,k);
        fieldattname = attname;
        variable.(varname).(fieldattname) = netcdf.getAtt(ncid,i,attname);
    end
end
netcdf.close(ncid)

inp = [];
inp.lat = variable.latitude.data;
inp.lon = variable.longitude.data;
inp.time = variable.time.data;
inp.angle = variable.angle.data;

if ifov_save.update
    
    outp = F_guess_ifov_60or120(inp);
    inpez = [];
    inpez.ifov = outp.ifov;
    inpez.angle = inp.angle;
    outpez = F_ez_ifov(inpez);
    
    data = [];
    data.lon = inp.lon(:);
    data.lat = inp.lat(:);
    data.fractional_day = outp.UTC(:);
    data.ifov = outpez.ifov(:);
    data.colnh3 = variable.column.data(:);
    data.error = variable.error.data(:);
    data.CLcov = variable.CLcov.data(:);
    data.angle = variable.angle.data(:);
    UTC = outp.UTC(:);
    ifov = outpez.ifov(:);
    save(ifov_save.ifovfn,'UTC','ifov')
else
    load(ifov_save.ifovfn)
    data = [];
    data.lon = inp.lon(:);
    data.lat = inp.lat(:);
    data.fractional_day = UTC;
    data.ifov = ifov;
    data.colnh3 = variable.column.data(:);
    data.error = variable.error.data(:);
    data.CLcov = variable.CLcov.data(:);
    data.angle = variable.angle.data(:);
end
%     ifovid = netcdf.defVar(ncid,'ifov','NC_FLOAT',dimids(1));
%     netcdf.endDef(ncid);
%     netcdf.putVar(ncid,ifovid,outpez.ifov);


