clear all

load geometries.mat
load final_times.mat
load dissociation_marker.mat
load dissociated_multiplicities.mat
load diss_times.mat
load cslabels.mat

geometries = permute(geometries,[1,2,4,3]);
[n_atoms, xyz, nts, ntraj] = size(geometries);
dt = 0.5;
time_grid = linspace(0,nts*dt,nts);
number_diss_traj = sum(diss_marker == 1);
diss_cutoff = 3.4;

distance_mat = zeros(size(geometries)); % calculate distance matrix 
for traj=1:ntraj
    for ts=1:nts
        for a=1:xyz
            for b=a+1:xyz
                dd = norm(geometries(a,1:3,ts,traj) - geometries(b,1:3,ts,traj));
            distance_mat(a,b,ts,traj) = dd;
            distance_mat(b,a,ts,traj) = dd;
            end
        end
    end
end

extended_distance_mat = distance_mat;
figure
for i=1:ntraj% loop over trajs
    cslabels(20) = 0;
    cs2_dist = distance_mat(2,1,:,i);
    cs3_dist = distance_mat(3,1,:,i);
    ss_dist = distance_mat(2,3,:,i); 
    cs2_dist = squeeze(cs2_dist)';
    cs3_dist = squeeze(cs3_dist)';
    ss_dist = squeeze(ss_dist)';
    last_idx = find(cs2_dist==0,1,'first'); % index at which no geom
   if cslabels(i) == 2
    diss_cutoff_idx = find(cs2_dist > 3.4,1,'first');
    [peaks, peaks_loc] = findpeaks(cs3_dist);  % find peaks/troughs 
    [troughs, troughs_loc] = findpeaks(-cs3_dist); % values and time index
    troughs = troughs*-1;
    num_peaks = length(peaks);
    num_troughs = length(troughs);
    peak_diffs = abs(diff(peaks)); % difference between consecutive
    trough_diffs = abs(diff(troughs)); % peaks/ troughs
    
    if cs3_dist(peaks_loc(end)+1) == 0 % remove last peak
        peaks(end) = [];
        peaks_loc(end) = [];
        num_peaks = length(peaks);
        peak_diffs = abs(diff(peaks));
    elseif cs3_dist(troughs_loc(end)+1) == 0 % remve last trough
        troughs(end) = [];
        troughs_loc(end) = [];
        num_troughs = length(troughs);
        trough_diffs = abs(diff(troughs));
    end
    
    A = (peaks(end)-troughs(end))/2; % ampltiude from final peak/trough
    
    if peak_diffs(end) < trough_diffs(end) % use peaks if smaller diff between two last peaks 
       T = (time_grid(peaks_loc(end))-time_grid(peaks_loc(end-1))); % Period
       wave_center = peaks(end)-A; % centre of wave
       omega = 2*pi/T; % angular frequencey
        freq = 1/T;
    elseif peak_diffs(end) > trough_diffs(end) % use troughs if smaller diff between two troughs
        T =  (time_grid(troughs_loc(end))-time_grid(troughs_loc(end-1)));
        wave_center = troughs(end)+A;
        omega = 2*pi/T;
        freq = 1/T;
    end
    
    deriv = abs((cs2_dist(last_idx)-cs2_dist(diss_cutoff_idx))/... % gradient between last value and 
        (time_grid(last_idx)-time_grid(diss_cutoff_idx))); % dissociated value for C-S1
    deriv_ss = abs((ss_dist(last_idx)-ss_dist(diss_cutoff_idx))/... % gradient for S-S 
        (time_grid(last_idx)-time_grid(diss_cutoff_idx)));
    c = cs2_dist(last_idx-1) - deriv*time_grid(last_idx-1); % intercept for C-S1
    c_ss = (ss_dist(last_idx-1) - deriv_ss*time_grid(last_idx-1)); % intercept for S-S
  
    for t=1:last_idx
        cs3_ext(t) = wave_center + A * sin(omega*time_grid(t));
    end
    [peaks_ext, peaks_loc_ext] = findpeaks(cs3_ext);
    dtau = (time_grid(peaks_loc(end)) - time_grid(peaks_loc_ext(end)));
    phi = omega*dtau;
    
    for t=1:last_idx
        cs3_extm(t) = wave_center + A * sin(pi-omega*time_grid(t)+phi);
        cs3_extp(t) = wave_center + A * sin(pi+omega*time_grid(t)+phi);
    end
    [peaks_extm, peaks_loc_extm] = findpeaks(cs3_extm);
    [peaks_extp, peaks_loc_extp] = findpeaks(cs3_extp);
    diff_p = abs(time_grid(peaks_loc(end))-time_grid(peaks_loc_extp(end)));
    diff_m = abs(time_grid(peaks_loc(end))-time_grid(peaks_loc_extm(end)));
    
    if diff_m < diff_p
        for t=last_idx:nts % loop over time from last complete value to end
        cs2_dist(t) = (time_grid(t)*deriv) + c; % linear fit for C-S1 distance
        cs3_dist(t) = wave_center + A * sin(pi-omega*time_grid(t) +phi); % harmonic fit for C-S2 distance
        ss_dist(t) = (time_grid(t)*deriv_ss) + c_ss; % linear fit for S-S distance
        end
    elseif diff_m > diff_p
        for t=last_idx:nts % loop over time from last complete value to end
        cs2_dist(t) = (time_grid(t)*deriv) + c; % linear fit for C-S1 distance
        cs3_dist(t) = wave_center + A * sin(pi+omega*time_grid(t) +phi); % harmonic fit for C-S2 distance
        ss_dist(t) = (time_grid(t)*deriv_ss) + c_ss; % linear fit for S-S distance
        end
    end
    
    hold on
    plot(time_grid(1:last_idx), cs3_dist(1:last_idx), 'b')
    plot(time_grid(last_idx:end), cs3_dist(last_idx:end), 'r')
    plot(time_grid(1:end), cs2_dist(1:end))
    plot(time_grid(1:end), ss_dist(1:end))
    xlabel('Time (fs)')
    xlim([0 1000])
    ylim([1 10])
    ylabel('C-S Bond Length')
    % Extend distance matrix with extrapolated values 
    extended_distance_mat(2,1,[last_idx:end],i) = cs2_dist(last_idx:end);
    extended_distance_mat(1,2,[last_idx:end],i) = cs2_dist(last_idx:end);
    extended_distance_mat(3,1,[last_idx:end],i) = cs3_dist(last_idx:end);
    extended_distance_mat(1,3,[last_idx:end],i) = cs3_dist(last_idx:end);
    extended_distance_mat(2,3,[last_idx:end],i) = ss_dist(last_idx:end);
    extended_distance_mat(3,2,[last_idx:end],i) = ss_dist(last_idx:end);
    
    for x = last_idx:ts
        cs_2 = cs2_dist(x);
        cs_3 = cs3_dist(x);
        ss_1 = ss_dist(x);
        rad = acos((cs_2^2+cs_3^2-ss_1^2)/(2*cs_2*cs_3));
        theta = real(rad * 180/pi);
        x_1 = real(cs_3*cos(rad));
        y_1 = sqrt(cs_3^2-x_1^2);
        geometries(2,1,i,x) = cs2_dist(x);
        geometries(3,1,i,x) = x_1;
        geometries(3,2,i,x) = y_1;
    end
    
   elseif cslabels(i) == 3
    diss_cutoff_idx = find(cs3_dist > 3.4,1,'first');
    [peaks, peaks_loc] = findpeaks(cs2_dist);  % find peaks/troughs 
    [troughs, troughs_loc] = findpeaks(-cs2_dist); % values and time index
    troughs = troughs*-1;
    num_peaks = length(peaks);
    num_troughs = length(troughs);
    peak_diffs = abs(diff(peaks)); % difference between consecutive
    trough_diffs = abs(diff(troughs)); % peaks/ troughs
    
    if cs2_dist(peaks_loc(end)+1) == 0 % remove last peak
        peaks(end) = [];
        peaks_loc(end) = [];
        num_peaks = length(peaks);
        peak_diffs = abs(diff(peaks));
    elseif cs2_dist(troughs_loc(end)+1) == 0 % remve last trough
        troughs(end) = [];
        troughs_loc(end) = [];
        num_troughs = length(troughs);
        trough_diffs = abs(diff(troughs));
    end
    
    A = (peaks(end)-troughs(end))/2; % ampltiude from final peak/trough
    
    if peak_diffs(end) < trough_diffs(end) % use peaks if smaller diff between two last peaks 
       T = (time_grid(peaks_loc(end))-time_grid(peaks_loc(end-1))); % Period
       wave_center = peaks(end)-A; % centre of wave
       omega = 2*pi/T; % angular frequencey
       freq = 1/T;
    elseif peak_diffs(end) > trough_diffs(end) % use troughs if smaller diff between two troughs
        T =  (time_grid(troughs_loc(end))-time_grid(troughs_loc(end-1)));
        wave_center = troughs(end)+A;
        omega = 2*pi/T;
        freq = 1/T;
    end
    
    deriv = abs((cs3_dist(last_idx)-cs3_dist(diss_cutoff_idx))/... % gradient between last value and 
        (time_grid(last_idx)-time_grid(diss_cutoff_idx))); % dissociated value for C-S1
    deriv_ss = abs((ss_dist(last_idx)-ss_dist(diss_cutoff_idx))/... % gradient for S-S 
        (time_grid(last_idx)-time_grid(diss_cutoff_idx)));
    c = cs3_dist(last_idx-1) - deriv*time_grid(last_idx-1); % intercept for C-S1
    c_ss = (ss_dist(last_idx-1) - deriv_ss*time_grid(last_idx-1)); % intercept for S-S
    
    for t=1:last_idx
       cs2_ext(t) = wave_center + A * sin(omega*time_grid(t));
    end
    [peaks_ext, peaks_loc_ext] = findpeaks(cs2_ext);
    dtau = (time_grid(peaks_loc(end)) - time_grid(peaks_loc_ext(end)));
    phi = omega*dtau;
    
    for t=1:last_idx
        cs2_extm(t) = wave_center + A * sin(pi-omega*time_grid(t)+phi);
        cs2_extp(t) = wave_center + A * sin(pi+omega*time_grid(t)+phi);
    end
    [peaks_extm, peaks_loc_extm] = findpeaks(cs2_extm);
    [peaks_extp, peaks_loc_extp] = findpeaks(cs2_extp);
    diff_p = abs(time_grid(peaks_loc(end))-time_grid(peaks_loc_extp(end)));
    diff_m = abs(time_grid(peaks_loc(end))-time_grid(peaks_loc_extm(end)));
    
    if diff_m < diff_p
        for t=last_idx:nts % loop over time from last complete value to end
        cs3_dist(t) = (time_grid(t)*deriv) + c; % linear fit for C-S1 distance
        cs2_dist(t) = wave_center + A * sin(pi-omega*time_grid(t) +phi); % harmonic fit for C-S2 distance
        ss_dist(t) = (time_grid(t)*deriv_ss) + c_ss; % linear fit for S-S distance
        end
    elseif diff_m > diff_p
        for t=last_idx:nts % loop over time from last complete value to end
        cs3_dist(t) = (time_grid(t)*deriv) + c; % linear fit for C-S1 distance
        cs2_dist(t) = wave_center + A * sin(pi+omega*time_grid(t) +phi); % harmonic fit for C-S2 distance
        ss_dist(t) = (time_grid(t)*deriv_ss) + c_ss; % linear fit for S-S distance
        end
    end
    
    hold on
    plot(time_grid(1:last_idx), cs2_dist(1:last_idx), 'b')
    plot(time_grid(last_idx:end), cs2_dist(last_idx:end), 'r')
    %plot(time_grid(1:end), cs3_dist(1:end))
    %plot(time_grid(1:end), ss_dist(1:end))
    xlabel('Time (fs)')
    xlim([0 1000])
    ylim([1 3])
    ylabel('C-S Bond Length')
    % Extend distance matrix with extrapolated values 
    extended_distance_mat(2,1,[last_idx:end],i) = cs2_dist(last_idx:end);
    extended_distance_mat(1,2,[last_idx:end],i) = cs2_dist(last_idx:end);
    extended_distance_mat(3,1,[last_idx:end],i) = cs3_dist(last_idx:end);
    extended_distance_mat(1,3,[last_idx:end],i) = cs3_dist(last_idx:end);
    extended_distance_mat(2,3,[last_idx:end],i) = ss_dist(last_idx:end);
    extended_distance_mat(3,2,[last_idx:end],i) = ss_dist(last_idx:end);
    
    for x = last_idx:ts
        cs_2 = cs2_dist(x);
        cs_3 = cs3_dist(x);
        ss_1 = ss_dist(x);
        rad = acos((cs_2^2+cs_3^2-ss_1^2)/(2*cs_2*cs_3));
        theta = real(rad * 180/pi);
        x_1 = real(cs_3*cos(rad));
        y_1 = sqrt(cs_3^2-x_1^2);
        geometries(2,1,i,x) = cs_2;
        geometries(3,1,i,x) = x_1;
        geometries(3,2,i,x) = y_1;
    end
    
   end
end

%geometries_extended = permute(geom_array,[1,2,4,3]);
save('geometries_extended', 'geometries')
