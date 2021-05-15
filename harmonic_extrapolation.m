% OUPUTS:
% fragment_bl - harmonic extrapolated values for the remaining C-S fragment
% dissociated_bl - linear extrapolated values for dissociated C-S bond
% ss_distance - linear extrapolated values for S-S distance

function [fragment_bl, dissociated_bl, ss_distance,diss_cutoff_idx] = harmonic_extrapolation(nts, time_grid, diss_cutoff, dissociated_bl, fragment_bl, ss_distance, last_idx)

    diss_cutoff_idx = find(dissociated_bl > diss_cutoff,1,'first'); % Index at which dissociation reached
    [peaks, peaks_loc] = findpeaks(fragment_bl);  % Find peaks/troughs of CS fragment stretching motion
    [troughs, troughs_loc] = findpeaks(-fragment_bl);
    troughs = troughs*-1;
    num_peaks = length(peaks);
    num_troughs = length(troughs);
    peak_diffs = abs(diff(peaks)); % Difference between consecutive peaks/troughs
    trough_diffs = abs(diff(troughs));
    
    if fragment_bl(peaks_loc(end)+1) == 0 % If last know geometry is at a peak value 
        peaks(end) = [];                  % remove last peak from list
        peaks_loc(end) = [];              % as this peak is not an actual peak
        num_peaks = length(peaks);
        peak_diffs = abs(diff(peaks));    % Update difference between peaks
    elseif fragment_bl(troughs_loc(end)+1) == 0 % If at trough - remove last
        troughs(end) = [];
        troughs_loc(end) = [];
        num_troughs = length(troughs);
        trough_diffs = abs(diff(troughs));      % Update differences
    end
    
    A = (peaks(end)-troughs(end))/2; % Ampltiude from final peak/trough
    
    % Check if peaks or troughs have a larger difference between last two
    % known values - use whichever has a smaller difference for following
    % calculations
    if peak_diffs(end) < trough_diffs(end) % Use peaks
        T = (time_grid(peaks_loc(end))-time_grid(peaks_loc(end-1))); % Calculate period
        wave_center = peaks(end)-A; % Calculate centre of wave
        omega = 2*pi/T; % Calculate angular frequencey
        freq = 1/T; % Calculate frequencey 
    elseif peak_diffs(end) > trough_diffs(end) % Use troughs
        T =  (time_grid(troughs_loc(end))-time_grid(troughs_loc(end-1)));
        wave_center = troughs(end)+A;
        omega = 2*pi/T;
        freq = 1/T;
    end
    
    % Calculate derivatives and intercepts for the dissociated and S-S
    % length - for use in a linear extrapolation (approximate atoms are
    % pulled apart with a speed that is constant between the dissociation
    % time and the simulation time limit)
    deriv = abs((dissociated_bl(last_idx-1)-dissociated_bl(last_idx-2))/... % Gradient for C-S 
        (time_grid(last_idx-1)-time_grid(last_idx-2)));
    deriv_ss = abs((ss_distance(last_idx-1)-ss_distance(last_idx-2))/... % Gradient for S-S 
        (time_grid(last_idx-1)-time_grid(last_idx-2)));
    c = dissociated_bl(last_idx-1) - deriv*time_grid(last_idx-1); % Intercept for dissociated C-S
    c_ss = (ss_distance(last_idx-1) - deriv_ss*time_grid(last_idx-1)); % Intercept for S-S
  
    for t=1:last_idx
        inverse_fragment(t) = wave_center + A * sin(omega*time_grid(t));
    end % Inverse harmonic extrapolation from last known geometry to t=0 geomtry 
    
    [~, peaks_loc_inverse] = findpeaks(inverse_fragment);
    dtau = (time_grid(peaks_loc(end)) - time_grid(peaks_loc_inverse(end)));
    phi = omega*dtau; % Calculate phase shift from time shift between known
                      % sinusoidal wave and the harmonic fit to the wave
    
    for t=1:last_idx % Check if sin(theta) = theta OR (pi - theta) is the best solution
        inverse_fragment_pi_minus(t) = wave_center + A * sin(pi-omega*time_grid(t)+phi);
        inverse_fragment_pi_plus(t) = wave_center + A * sin(pi+omega*time_grid(t)+phi);
    end 
    [~, peaks_loc_inverse_pi_minus] = findpeaks(inverse_fragment_pi_minus);
    [~, peaks_loc_inverse_pi_plus] = findpeaks(inverse_fragment_pi_plus);
    diff_pi_plus = abs(time_grid(peaks_loc(end))-time_grid(peaks_loc_inverse_pi_plus(end)));
    diff_pi_minus = abs(time_grid(peaks_loc(end))-time_grid(peaks_loc_inverse_pi_minus(end)));
    
    % Extrapolate to final time using the solution of sin(theta) that give the best
    % match to the known sinusoidal wave
    if diff_pi_minus < diff_pi_plus 
        for t=last_idx:nts
            dissociated_bl(t) = (time_grid(t)*deriv) + c; % Linear fit
            fragment_bl(t) = wave_center + A * sin(pi-omega*time_grid(t) +phi); % Harmonic fit
            ss_distance(t) = (time_grid(t)*deriv_ss) + c_ss; % Linear fit
        end
    elseif diff_pi_minus > diff_pi_plus
        for t=last_idx:nts 
            dissociated_bl(t) = (time_grid(t)*deriv) + c; % Linear fit
            fragment_bl(t) = wave_center + A * sin(pi+omega*time_grid(t) +phi); % Harmonic fit
            ss_distance(t) = (time_grid(t)*deriv_ss) + c_ss; % Linear fit     
        end
    elseif diff_pi_minus == diff_pi_plus
        for t=last_idx:nts
            dissociated_bl(t) = (time_grid(t)*deriv) + c; % Linear fit
            fragment_bl(t) = wave_center + A * sin(pi+omega*time_grid(t) +phi); % Harmonic fit
            ss_distance(t) = (time_grid(t)*deriv_ss) + c_ss; % Linear fit
        end
    end

end

