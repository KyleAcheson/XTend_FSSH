function [] = plot_trajectories(time_grid, fragment_bl, dissociated_bl, ss_distance, diss_cutoff_idx, last_idx, nts, dt)
    hold on 
    plot(time_grid(1:last_idx), fragment_bl(1:last_idx), 'b')
    plot(time_grid(last_idx:end), fragment_bl(last_idx:end), 'r')
    %plot(time_grid(diss_cutoff_idx), dissociated_bl(diss_cutoff_idx), 'r*')
    %plot(time_grid(diss_cutoff_idx), ss_distance(diss_cutoff_idx), 'r*')
    %plot(time_grid(last_idx), dissociated_bl(last_idx), 'r*')
    %plot(time_grid(last_idx), ss_distance(last_idx), 'r*')
    plot(time_grid(1:end), dissociated_bl(1:end))
    plot(time_grid(1:end), ss_distance(1:end))
    xlabel('Time (fs)')
    xlim([0 (nts-1)*dt])
    ylim([1 10])
    ylabel('Bond Length')
end

