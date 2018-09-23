function data = F_manipulate_iasi_nc(fn,ifov_save)

% not really 'manipulate' because the iasi nc files are read-only. Have to
% calculate ifov and save separately.

% written by Kang Sun on 2017/12/28
% updated by Kang Sun on 2018/09/23 to accomodate the new v2.2 and v2.2 r
% versions, which is messy and lack of consistency

% addpath('c:\Users\Kang Sun\Dropbox\matlab functions\')
% addpath('C:\Users\Kang Sun\Documents\GitHub\PU_KS_share\')
% fn = 'C:\data_ks\IASIa\L2\IASI_metopa_L2_NH3_20071001_V2.1_AM_GLOBAL.nc';
% % fn = 'C:\data_ks\IASIa\L2\IASI_metopa_L2_NH3_20160701_V2.1_AM_GLOBAL.nc';
% clc

if ~isstruct(fn)
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
else
    ncid = netcdf.open(fn.fn_v22);
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
    
    ncid = netcdf.open(fn.fn_v22r);
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    variabler = [];
    for i = 0:nvars-1
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,i);
        variabler.(varname).data = netcdf.getVar(ncid,i);
        for k = 0:natts-1
            attname = netcdf.inqAttName(ncid,i,k);
            fieldattname = attname;
            variabler.(varname).(fieldattname) = netcdf.getAtt(ncid,i,attname);
        end
    end
    netcdf.close(ncid)
end

if ~isstruct(fn)
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
    
else
    
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
        data.lon = variabler.longitude.data(:);
        data.lat = variabler.latitude.data(:);
        f_all = ismember([variable.longitude.data,variable.latitude.data],...
            [variabler.longitude.data,variabler.latitude.data],'rows');
        if sum(f_all) ~= length(variabler.time.data)
            error('I cannot perfectly map v2.2r pixels to v2.2 pixels!')
        end
        data.fractional_day = outp.UTC(f_all);
        data.ifov = outpez.ifov(f_all);
        data.colnh3 = variabler.column.data(:);
        data.error = variabler.error.data(:);
        data.CLcov = variabler.CLcov.data(:);
        data.angle = variabler.angle.data(:);
        UTC = outp.UTC(f_all);
        ifov = outpez.ifov(f_all);
        save(ifov_save.ifovfn,'UTC','ifov')
    else
        load(ifov_save.ifovfn)
        data = [];
        data.lon = inp.lon(:);
        data.lat = inp.lat(:);
        data.fractional_day = UTC;
        data.ifov = ifov;
        data.colnh3 = variabler.column.data(:);
        data.error = variabler.error.data(:);
        data.CLcov = variabler.CLcov.data(:);
        data.angle = variabler.angle.data(:);
    end
end

