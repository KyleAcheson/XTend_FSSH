% KYLE ACHESON - Script to extend SHARC CS2 trajectories that have crashed
% after the dissociation cutoff has been reached, but before the simulation
% time limit. 
% REQUIRES - geometries.mat array with dimensions (n_atoms, 3, ntraj, nts) 
% cslabels.mat vector (1, ntraj) containing labels of which atom
% dissociated. final_time.mat containing times at which calculations
% crashed. USE trajectory_collector.py to compile these.
% OUTPUT - geometries_extended.mat - contains extended trajectories (xyz).
% distances_extended.mat - contains extended distance matrix.
% NOTE - If using for generating non-rotationally averaged scattering
% patterns, may be a source of error stemming setting C1 to (0,0,0) and
% assuming molecule is in XY plane.
% EDIT trajectory loop to visually inspect individual trajectories.

clear all

load geometries_2.mat
load final_t_2.mat
load cslabels_2.mat

%EDIT
dt = 0.5; % time step 
diss_cutoff = 3.4; % limit of C-S bond length at which dissociation occurs 
FLAGplot = 1; % 1 = plots n trajectories, 0 = does not plot 
%%%%%%%%%%%%
geometries = geom_array;
geometries = permute(geometries,[1,2,4,3]);
[n_atoms, xyz, ntraj, nts] = size(geometries);
time_grid = linspace(0,nts*dt,nts);
end_time = (nts-1)*dt;

distance_mat = zeros(size(geometries)); % calculate distance matrix 
for traj=1:ntraj
    for ts=1:nts
        for a=1:xyz
            for b=a+1:xyz
                dd = norm(geometries(a,1:3,traj,ts) - geometries(b,1:3,traj,ts));
            distance_mat(a,b,traj,ts) = dd;
            distance_mat(b,a,traj,ts) = dd;
            end
        end
    end
end
geometries_old = geometries;
D_old = distance_mat;

if FLAGplot == 1
figure
end
for i=1:ntraj% loop over trajs - EDIT for visual inspection

    % Initialise bond length variables
    cs2_distance = distance_mat(2,1,i,:);
    cs3_distance = distance_mat(3,1,i,:);
    ss_distance = distance_mat(2,3,i,:); 
    cs2_distance = squeeze(cs2_distance)';
    cs3_distance = squeeze(cs3_distance)';
    ss_distance = squeeze(ss_distance)';
    last_idx = find(cs2_distance==0,1,'first'); % Index at which no geom

   % If atom 2 (first sulphur) dissociated
   if cslabels(i) == 2 && final_t(i) ~= end_time
       dissociated_bl = cs2_distance;
       fragment_bl = cs3_distance;
       % Calculate extrapolated values to t=nts using a harmonic oscillator
       % approximation for the remaining fragment and a linear fit for the
       % dissociated atom
       [fragment_bl, dissociated_bl, ss_distance, diss_cutoff_idx] = harmonic_extrapolation...
           (nts, time_grid, diss_cutoff, dissociated_bl, fragment_bl, ss_distance, last_idx);
       
       if FLAGplot == 1
           plot_trajectories(time_grid, fragment_bl, dissociated_bl, ss_distance, diss_cutoff_idx, last_idx, nts, dt);
       end
       
       % Update distance matrix with extrapolated values
       distance_mat(2,1,i,[last_idx:end]) = dissociated_bl(last_idx:end);
       distance_mat(1,2,i,[last_idx:end]) = dissociated_bl(last_idx:end);
       distance_mat(3,1,i,[last_idx:end]) = fragment_bl(last_idx:end);
       distance_mat(1,3,i,[last_idx:end]) = fragment_bl(last_idx:end);
       distance_mat(2,3,i,[last_idx:end]) = ss_distance(last_idx:end);
       distance_mat(3,2,i,[last_idx:end]) = ss_distance(last_idx:end);
       
       % Convert from distance matrix to cartesian coordinates, Carbon is
       % set to the origin (0,0,0) - molecule in XY plane
       
       for ts = last_idx:last_idx+1
           cs_2 = dissociated_bl(ts);
           cs_3 = fragment_bl(ts);
           ss_1 = ss_distance(ts);
           rad = real(acos((cs_2^2+cs_3^2-ss_1^2)/(2*cs_2*cs_3)));
           theta = real(rad * 180/pi);
           x_1 = real(cs_3*cos(rad));
           y_1 = sqrt(cs_2^2-x_1^2);
           geometries(2,1,i,ts) = x_1;
           geometries(3,1,i,ts) = cs_3;
           geometries(2,2,i,ts) = y_1;
       end
       unit_vector = geometries(2,[1:2],i,last_idx+1)-geometries(2,[1:2],i,last_idx);
       unit_vector = [unit_vector 0];
       for ts=last_idx+2:nts
           geometries(3,1,i,ts) = fragment_bl(ts);
           geometries(2,:,i,ts) = geometries(2,:,i,ts-1)+unit_vector;
       end
       
       
   % If atom 3 (second sulphur) dissociated
   elseif cslabels(i) == 3 && final_t(i) ~= end_time
       dissociated_bl = cs3_distance;
       fragment_bl = cs2_distance;
       % Do extrapolation
       [fragment_bl, dissociated_bl, ss_distance, diss_cutoff_idx] = harmonic_extrapolation...
           (nts, time_grid, diss_cutoff, dissociated_bl, fragment_bl, ss_distance, last_idx);
       
       if FLAGplot == 1
           plot_trajectories(time_grid, fragment_bl, dissociated_bl, ss_distance, diss_cutoff_idx, last_idx, nts, dt);
       end
       
       % Update matrices
       distance_mat(2,1,i,[last_idx:end]) = dissociated_bl(last_idx:end);
       distance_mat(1,2,i,[last_idx:end]) = dissociated_bl(last_idx:end);
       distance_mat(3,1,i,[last_idx:end]) = fragment_bl(last_idx:end);
       distance_mat(1,3,i,[last_idx:end]) = fragment_bl(last_idx:end);
       distance_mat(2,3,i,[last_idx:end]) = ss_distance(last_idx:end);
       distance_mat(2,3,i,[last_idx:end]) = ss_distance(last_idx:end);
       
       % Convert to cartesian - same as above
       for ts = last_idx:last_idx+1
           cs_3 = dissociated_bl(ts);
           cs_2 = fragment_bl(ts);
           ss_1 = ss_distance(ts);
           rad = real(acos((cs_2^2+cs_3^2-ss_1^2)/(2*cs_2*cs_3)));
           theta = real(rad * 180/pi);
           x_1 = real(cs_2*cos(rad));
           y_1 = sqrt(cs_3^2-x_1^2);
           geometries(2,1,i,ts) = cs_2;
           geometries(3,1,i,ts) = x_1;
           geometries(3,2,i,ts) = y_1;
       end
       unit_vector = geometries(3,[1:2],i,last_idx+1)-geometries(3,[1:2],i,last_idx);
       unit_vector = [unit_vector 0];
       for ts=last_idx+2:nts
           geometries(3,:,i,ts) = geometries(3,:,i,ts-1)+unit_vector;
           geometries(2,1,i,ts) = fragment_bl(ts);
       end
       
   elseif cslabels(i) ==0
       geometries(:,:,i,:) = geometries(:,:,i,:);
   end
end

for traj=1:ntraj
    for ts=1:nts
        for a=1:xyz
            for b=a+1:xyz
                dd = norm(geometries(a,1:3,traj,ts) - geometries(b,1:3,traj,ts));
            distance_mat(a,b,traj,ts) = dd;
            distance_mat(b,a,traj,ts) = dd;
            end
        end
    end
end
figure
fragment_bl_total= zeros(1,nts);
dissociated_bl_total = zeros(1,nts);
for i=1:ntraj
    cs2_distance = distance_mat(2,1,i,:);
    cs3_distance = distance_mat(3,1,i,:);
    ss_distance = distance_mat(2,3,i,:); 
    cs2_distance = squeeze(cs2_distance)';
    cs3_distance = squeeze(cs3_distance)';
    ss_distance = squeeze(ss_distance)'; 
    
    if cslabels(i) == 2 && final_t(i) ~= end_time
       dissociated_bl = cs2_distance;
       fragment_bl = cs3_distance;
       fragment_bl_total(1,:) = fragment_bl_total(1,:) + fragment_bl(1,:);
       dissociated_bl_total(1,:) = dissociated_bl_total(1,:) + dissociated_bl(1,:);
       hold on
       plot(time_grid(1:last_idx), fragment_bl(1:last_idx), 'b')
        plot(time_grid(last_idx:end), fragment_bl(last_idx:end), 'r')
       plot(time_grid(1:end), dissociated_bl(1:end))
       plot(time_grid(1:end), ss_distance(1:end))
       ylim([1 10])
       ylabel('Bond Length')
       xlabel('Time (fs)')
       xlim([0 1000])
    elseif cslabels(i) == 3 && final_t(i) ~= end_time
       dissociated_bl = cs3_distance;
       fragment_bl = cs2_distance;
       fragment_bl_total(1,:) = fragment_bl_total(1,:) + fragment_bl(1,:);
       dissociated_bl_total(1,:) = dissociated_bl_total(1,:) + dissociated_bl(1,:);
       hold on
       plot(time_grid(1:last_idx), fragment_bl(1:last_idx), 'b')
       plot(time_grid(last_idx:end), fragment_bl(last_idx:end), 'r')
       plot(time_grid(1:end), dissociated_bl(1:end))
       plot(time_grid(1:end), ss_distance(1:end), 'g')
    end
end
fragment_bl_avg = fragment_bl_total./ntraj;
dissociated_bl_avg = dissociated_bl_total./ntraj;

%geometries = real(geometries);
%save('geometries_extended_231', 'geometries')
%save('distances_extended', 'distance_mat')
