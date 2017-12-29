% matlab script to subset v2.1 iasi data, downloaded using read_iasi_v2.m.
% Turn on 'update_ifov' to generate ifov using F_guess_ifov(incomplete).m
% Written by Kang Sun on 2017/12/29
clear;clc
if ispc
    L2dir = 'C:\data_ks\IASIa\L2\';
    L2gdir = 'C:\data_ks\IASIa\L2g\';
    cd('C:\Users\Kang Sun\Documents\GitHub\PU_KS_share\')
    pixel_size_file = '.\daysss.mat';
else
    L2dir = '/data/tempo1/Shared/kangsun/IASI/IASIa/L2/';
    L2gdir = '/data/tempo1/Shared/kangsun/IASI/IASIa/L2g/';
    cd('/home/kangsun/IASI/PU_KS_share/')
    pixel_size_file = './daysss.mat';
end
%% subsetting L2 data, using F_subset_IASI.m. Save the results every year
pixel_shape = load(pixel_size_file);
inp_subset = [];
% CONUS
inp_subset.MinLat = 25;
inp_subset.MaxLat = 50;
inp_subset.MinLon = -130;
inp_subset.MaxLon = -63;
inp_subset.MaxCF = 0.25;
inp_subset.L2dir = L2dir;

inp_subset.a_or_b = 'a';
inp_subset.update_ifov = true;

for iyear = 2007:2016
    inp_subset.Startdate = [iyear 1 1];
    inp_subset.Enddate = [iyear 12 31];
%     inp_subset.Startdate = [iyear 11 18];
%     inp_subset.Enddate = [iyear 11 18];
    output_subset = F_subset_IASI_v2(inp_subset);
    % calculate pixel shape
    [output_subset.u, output_subset.v, output_subset.t] =...
        F_define_IASI_pixel(output_subset.lat,output_subset.ifov,...
        pixel_shape.uuu4,pixel_shape.vvv4,pixel_shape.ttt4);
    L2g_fn = ['CONUS_v2_',num2str(iyear),'.mat'];
    save([L2gdir,L2g_fn],'inp_subset','output_subset')
end
%%
% scatter(output_subset.lon,output_subset.lat,[],output_subset.utc);colorbar