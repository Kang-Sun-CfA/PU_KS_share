function outp = F_guess_ifov(inp)

% tried really hard to get ifov if it is not provided (WHY DOING
% THIS!!@@!!!!!!!). Test code is in read_v21.m, try_ifov_from_angle.m

% when a single sounding is missing in a scan position, the whole group of
% four is disgarded. Too hard to save them. Try lat lon information if you
% are really into this. <- Upgraded on 12/21. Now they are saved.

% written by Kang Sun from 12/18 to 12/20/2017.
% added option to estimate ifov with missing data (error = +/- 1)
% Improved on 12/21/2017 to nail down each ifov at mid-low latitude

% fixed a few bugs on 2017/12/28

angle = inp.angle(:);
lat = inp.lat(:);
lon = inp.lon(:);
time = inp.time(:);

ttbdr = [-3.7509395e-05,-3.5051035e-05,-3.2507902e-05,-3.0009985e-05,-2.7506412e-05,-2.5008492e-05,-2.2504919e-05,-2.0006999e-05,-1.7503429e-05,-1.5005508e-05,-1.2462377e-05,-1.0004016e-05,-7.5004441e-06,-4.9573127e-06,-2.4989522e-06,-1.0319354e-09,2.5420995e-06,5.0004601e-06,7.5435919e-06,1.0047163e-05,1.2505524e-05,1.5003444e-05,1.7546576e-05,2.0004936e-05,2.2548069e-05,2.5051639e-05,2.7549559e-05,3.0053132e-05,3.2551052e-05,3.5054625e-05,3.7558195e-05];

anglebdr = [0,2.0300002,3.8375001,5.6925001,7.5725002,9.4624996,11.352500,13.242500,15.142500,17.042500,18.932499,20.834999,22.762499,24.689999,26.609999,28.537498,30.485001,32.439999,34.415001,36.397499,38.395000,40.407501,42.442497,44.494999,46.587502,48.702499,50.857502,53.044998,55.292500,57.582500,Inf];
[~,anglen] = histc(angle,anglebdr);
ref_anglen = [29,29,30,30,27,27,28,28,25,25,26,26,23,23,24,24,21,21,22,22,19,19,20,20,17,17,18,18,15,15,16,16,13,13,14,14,11,11,12,12,9,9,10,10,7,7,8,8,5,5,6,6,3,3,4,4,1,1,2,2,2,2,1,1,4,4,3,3,6,6,5,5,8,8,7,7,10,10,9,9,12,12,11,11,14,14,13,13,16,16,15,15,18,18,17,17,20,20,19,19,22,22,21,21,24,24,23,23,26,26,25,25,28,28,27,27,30,30,29,29];
if isempty(strfind(anglen(:)',ref_anglen))
    error('120 fovs are not complete!')
end
ref_anglen_pm = ref_anglen;
ref_anglen_pm(61:end) = -ref_anglen_pm(61:end);
ifov = nan*anglen;
HH = double(floor(time/10000));
MM = double(floor((time-10000*HH)/100));
SS = double(mod(time,100));
UTC = HH(:)/24+MM(:)/60/24+SS(:)/3600/24;

difflat = diff(lat);

ind_orbit = find(difflat > 90);

if isempty(ind_orbit) % what if there is only one (or less) orbit of this day?!
    int_orbit_cell = [];
    int_orbit_cell{1} = single(1:length(lat));
elseif length(ind_orbit) == 1 % two orbits?
    int_orbit_cell = [];
    int_orbit_cell{1} = single(1:ind_orbit);
    int_orbit_cell{2} = single(ind_orbit+1:length(lat));
else
    if length(lat)-ind_orbit(end) > 5000
        ind_orbit = [ind_orbit;length(lat)]; % disgard incomplete orbit
    end
    if ind_orbit(1) > 5000
        ind_orbit = [0;ind_orbit];
    end
    
    int_orbit_cell = cell(length(ind_orbit)-1,1);
    for i = 1:length(ind_orbit)-1
        int_orbit_cell{i} = single((ind_orbit(i)+1):ind_orbit(i+1));
    end
end

int_orbit_cell = cell(length(ind_orbit)-1,1);
for i = 1:length(ind_orbit)-1
    int_orbit_cell{i} = single((ind_orbit(i)+1):ind_orbit(i+1));
end

position_array = 1:30;
array1 = 3:4:120;
array2 = 1:4:120;
array3 = 4:4:120;
array4 = 2:4:120;
upper_ind = nan(60,1);
lower_ind = nan(60,1);
upper_ind(1:2:end) = array3;
upper_ind(2:2:end) = array2;
lower_ind(1:2:end) = array1;
lower_ind(2:2:end) = array4;

for iorbit = 1:length(int_orbit_cell)
    
    orbit_anglen = anglen(int_orbit_cell{iorbit});
    orbit_anglen = orbit_anglen(:)';
    pattern_str = [1 1 2 2 2 2 1 1];
    center_ind = strfind(orbit_anglen,pattern_str)-1+int_orbit_cell{iorbit}(1);
    center_time = center_ind*nan;
    for icenter_ind = 1:length(center_ind)
        center_time(icenter_ind) = mean(UTC(center_ind(icenter_ind):center_ind(icenter_ind)+length(pattern_str)-1));
    end
    scan_time = median(diff(center_time)); % precisely, 8 seconds per full scan
    orbit_start = min(UTC(int_orbit_cell{iorbit}));
    orbit_end = max(UTC(int_orbit_cell{iorbit}));
    if iorbit ~= length(int_orbit_cell)
    % add missing center time at the beginning of orbit
    count = 0;
    while center_time(1) > orbit_start+scan_time
        center_time = cat(2,center_time(1)-scan_time,center_time);
        count = count+1;
    end
    % add missing center time at the end of orbit
    count = 0;
    while center_time(end) < orbit_end-scan_time
        center_time = cat(2,center_time,center_time(end)+scan_time);
        count = count+1;
    end
    end
    miss_ind = find(diff(center_time) > scan_time*1.5);
    
    center_time1 = center_time;
    total_shift = 0;
    for imiss = miss_ind
        nmiss = round((center_time(imiss+1)-center_time(imiss))/scan_time)-1;
        miss_str = center_time(imiss)+(1:nmiss)*scan_time;
        center_time1 = cat(2,center_time1(1:imiss+total_shift),miss_str,center_time1(imiss+1+total_shift:end));
        total_shift = total_shift+nmiss;
    end
    
    full_scan_ind = strfind(orbit_anglen,ref_anglen)-1+int_orbit_cell{iorbit}(1);
    center_time2 = center_time1;
    
    for ifull_scan = 1:length(full_scan_ind)
        full_scan_int = full_scan_ind(ifull_scan):full_scan_ind(ifull_scan)+119;
        if mean(UTC(full_scan_int(60:61)))-UTC(full_scan_int(1)) < abs(ttbdr(1)) ...
            && UTC(full_scan_int(end))-mean(UTC(full_scan_int(60:61))) < abs(ttbdr(end))
        ifov(full_scan_int) = 1:120;
        
        [mintdiff,I] = min(abs(mean(UTC(full_scan_int(60:61)))-center_time1));
        if mintdiff > 0.1/86400
            error('time matching went wrong')
        end
        center_time2(I) = nan;
        end
    end
    center_time2 = center_time2(~isnan(center_time2));
    
    % iscan = 200;
    % scan_center_time = center_time1(iscan);
    % scan_UTC = UTC(UTC > scan_center_time-3.5/86400 & UTC < scan_center_time+3.5/86400);
    % position_time_all = mean(reshape(scan_UTC,[4,30]))-scan_center_time;
    % position_time = median(diff(position_time_all)); % position time, 0.2163 s
    % ttbdr = [position_time_all-0.5*position_time,position_time_all(end)+0.5*position_time];
    % orbit_UTC = UTC(int_orbit_cell{iorbit});
    for iscan = 1:length(center_time2)
        scan_center_time = center_time2(iscan);
        [binc,C] = histc(UTC,ttbdr+scan_center_time);
        full_position = find(binc == 4);
        full_fov = [full_position*4-3;full_position*4-2;full_position*4-1;full_position*4];
        full_fov = sort(full_fov);
        full_int = ismember(C,full_position);
        ifov(full_int) = full_fov;
        lon_full = lon(full_int);
        lat_full = lat(full_int);
        lon_up = lon_full(upper_ind(upper_ind <= length(lon_full)));
        lat_up = lat_full(upper_ind(upper_ind <= length(lon_full)));
        lon_lo = lon_full(lower_ind(lower_ind <= length(lon_full)));
        lat_lo = lat_full(lower_ind(lower_ind <= length(lon_full)));
        lon_ct = 0.5*(lon_up+lon_lo);
        lat_ct = 0.5*(lat_up+lat_lo);

        for iposition = position_array(binc(1:30)' > 0 & binc(1:30)' < 4)
            int_bad_position = find(C == iposition);
                local_lat = lat(int_bad_position);
                          local_lon = lon(int_bad_position);
                local_anglen = anglen(int_bad_position);
                if iposition > 15
                    local_anglen = -local_anglen;
                end
            if max(abs(local_lat)) > 60 % do not use lat info for high latitude 
                guess_ind = strfind(ref_anglen_pm(:)',local_anglen(:)');
                if isempty(guess_ind)
                    warning('The angle pattern should not be like this, not sure why')
                else
                ifov(int_bad_position) = guess_ind(1):(guess_ind(1)+length(local_anglen)-1);
                end 
            else
                for ip = 1:binc(iposition)
                    p_anglen = local_anglen(ip);
                    p_anglen_s = p_anglen;
                    if p_anglen_s < 0
                        p_anglen_s = p_anglen_s+1;
                    end
                    guess_ind = strfind(ref_anglen_pm(:)',p_anglen(:)');
                    if length(guess_ind) == 1
                        ifov(int_bad_position(ip)) = guess_ind;
                    else
                    p_lat = local_lat(ip);
                    p_lon = local_lon(ip);
                    c_lat = interp1(lon_ct,lat_ct,p_lon,'linear','extrap');
                    above_or_below = p_lat > c_lat;
                    if (mod(p_anglen_s,2) == 0 && above_or_below) || ...
                            (mod(p_anglen_s,2) == 1 && ~above_or_below)
                    ifov(int_bad_position(ip)) = guess_ind(2);
                    else
                        ifov(int_bad_position(ip)) = guess_ind(1);
                    end
                    end
                end

            end
        end
    end
    disp(['Orbit # ',num2str(iorbit),' finished...'])
end
outp.UTC = UTC;
outp.ifov = ifov;