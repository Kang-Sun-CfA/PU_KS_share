function [u, v, t] = F_define_CrIS_pixel(lat,ifov,uuu4,vvv4,ttt4)
% calculate ellipse pixel with look up table
ifov = round(ifov);
u = single(interp2(-90:90,1:270,uuu4',lat,ifov));
v = single(interp2(-90:90,1:270,vvv4',lat,ifov));
t = single(interp2(-90:90,1:270,ttt4',lat,ifov));