function generate_plots_variance_topk_coordinates(x, ...
                                                  true_variances, ...
                                                  subsampling_mean, ...
                                                  subsampling_std, ...
                                                  bootstrap_mean, ...
                                                  bootstrap_std)
    %% Create a new figure with a white background and set its size
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    hold on;
    
    %% Enhance the current axes properties for a nicer look.
    ax = gca;
    ax.FontSize   = 16;        % Increase tick label font size
    ax.LineWidth  = 1.5;       % Thicken the axis lines
    ax.TickDir    = 'out';     % Ticks on the outside for a cleaner look
    ax.Box        = 'off';     % Remove the top and right box lines
    grid on;
    grid minor;
    
    %% Plot the true variances as a line plot.
    plot(x, true_variances(x), 'k-', 'LineWidth', 3, 'DisplayName', 'True Variance');
    
    %% Plot the subsampling mean and its uncertainty as a shaded region.
    % Construct the x-values for the patch.
    x_patch = [x, fliplr(x)];
    % Construct the y-values for the shaded region (mean +/- std),
    % ensuring the lower bound is at least 0.
    lower_sub = max(subsampling_mean(x) - subsampling_std(x), 0);
    upper_sub = subsampling_mean(x) + subsampling_std(x);
    y_patch_sub = [lower_sub, fliplr(upper_sub)];
    % Plot the shaded uncertainty region.
    fill(x_patch, y_patch_sub, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
         'DisplayName', 'OjaVarEst Uncertainty');
    % Plot the subsampling mean over the shaded region.
    plot(x, subsampling_mean(x), 'b-', 'LineWidth', 3, 'DisplayName', 'OjaVarEst Mean');
    
    %% Plot the bootstrap mean and its uncertainty as a shaded region.
    % Construct the y-values for the shaded region.
    lower_boot = max(bootstrap_mean(x) - bootstrap_std(x), 0);
    upper_boot = bootstrap_mean(x) + bootstrap_std(x);
    y_patch_boot = [lower_boot, fliplr(upper_boot)];
    % Plot the shaded uncertainty region.
    fill(x_patch, y_patch_boot, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
         'DisplayName', 'Bootstrap Uncertainty');
    % Plot the bootstrap mean over the shaded region.
    plot(x, bootstrap_mean(x), 'r-', 'LineWidth', 3, 'DisplayName', 'Bootstrap Mean');
    
    %% Add labels and title with increased font sizes and bold styling.
    xlabel('Coordinate', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Variance', 'FontSize', 18, 'FontWeight', 'bold');
    title('Variance Estimates with Uncertainty Regions', 'FontSize', 20, 'FontWeight', 'bold');
    
    %% Display the legend with larger font size and improved location.
    legend('show', 'FontSize', 14, 'Location', 'best');
    
    hold off;
end