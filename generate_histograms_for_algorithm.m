function generate_histograms_for_algorithm(algorithmName, topCoords, ojaVecResults, varianceResults, data_params, commonXLimits)
% generate_histograms_for_algorithm plots, for each coordinate in topCoords, an empirical
% density (displayed as a shaded, staircase‐shaped patch with an overlaid stairs line)
% of the oja_vec estimates along with an overlaid Gaussian density. The Gaussian density 
% has mean equal to the true coordinate value (data_params.trueV(i)) and variance taken 
% either from a randomly chosen experiment or from the mean variance (if commonXLimits is provided).
%
% INPUTS:
%   algorithmName   - A string ('ojavarest' or 'bootstrap') used for titles.
%   topCoords       - A vector of coordinate indices to plot (e.g., 1:3).
%   ojaVecResults   - A [num_experiments x d] matrix of estimated oja_vec values.
%   varianceResults - A [num_experiments x d] matrix of estimated variances.
%   data_params     - A structure containing at least the field trueV.
%   commonXLimits   - (Optional) A cell array with one element per coordinate.
%                     Each cell is a two-element vector [xmin, xmax] that sets the
%                     x-axis limits for that coordinate. When provided, both algorithms
%                     will be plotted on the same scale.
%
% The function creates one figure per coordinate, formatted for publication.

    % Number of experiments (assumed to be the number of rows)
    num_experiments = size(ojaVecResults, 1);
    
    % Choose a face color based on the algorithm type.
    switch lower(algorithmName)
        case 'ojavarest'
            faceColor = [0.2, 0.6, 0.8];  % bluish
        case 'bootstrap'
            faceColor = [0.8, 0.4, 0.4];  % reddish
        otherwise
            faceColor = [0.5, 0.5, 0.5];  % gray (default)
    end
    
    % Loop over each coordinate to create separate figures.
    for idx = 1:length(topCoords)
        i = topCoords(idx);
        
        % Create a new figure with a white background.
        figure('Color','w','Position',[100 100 700 500]);
        ax = gca;
        ax.FontSize = 16;
        ax.LineWidth = 1.5;
        hold on;
        
        % Extract the data for coordinate i.
        data_i = ojaVecResults(:, i);
        
        % Determine the x-range and variance for the Gaussian overlay.
        if nargin < 6 || isempty(commonXLimits)
            % Use a randomly chosen experiment's variance.
            r = randi(num_experiments);
            var_i = varianceResults(r, i);
            x_min = data_params.trueV(i) - 4*sqrt(var_i);
            x_max = data_params.trueV(i) + 4*sqrt(var_i);
        else
            % Use the provided common x-limits.
            x_min = commonXLimits{idx}(1);
            x_max = commonXLimits{idx}(2);
            var_i = mean(varianceResults(:, i));
            xlim([x_min, x_max]);
        end
        
        % Compute the empirical density using histogram counts.
        [counts, edges] = histcounts(data_i, 'Normalization', 'pdf');
        nBins = length(counts);
        
        % Build patch coordinates for a staircase shape.
        % The patch will start at (edges(1), 0), then for each bin it will
        % go from (edges(j), counts(j)) to (edges(j+1), counts(j)), and finally back to (edges(end),0).
        x_patch = zeros(1, 2*nBins + 2);
        y_patch = zeros(1, 2*nBins + 2);
        x_patch(1) = edges(1);
        y_patch(1) = 0;
        for j = 1:nBins
            x_patch(2*j)   = edges(j);
            y_patch(2*j)   = counts(j);
            x_patch(2*j+1) = edges(j+1);
            y_patch(2*j+1) = counts(j);
        end
        x_patch(end) = edges(end);
        y_patch(end) = 0;
        
        % Plot the empirical density as a filled patch (shaded area).
        % Hide its handle from the legend.
        fill(x_patch, y_patch, faceColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility','off');
        
        % Plot the staircase line on top.
        stairs(edges, [counts, counts(end)], 'Color', faceColor, 'LineWidth', 2, ...
               'DisplayName', 'Empirical Density');
        
        % Generate the x-range for the Gaussian overlay.
        x_range = linspace(x_min, x_max, 100);
        % Compute the Gaussian density with mean = trueV(i) and std = sqrt(var_i).
        y_pdf = normpdf(x_range, data_params.trueV(i), sqrt(var_i));
        
        % Overlay the Gaussian density.
        plot(x_range, y_pdf, 'r-', 'LineWidth', 2, 'DisplayName', 'Gaussian Density');
        
        % Mark the true value with a solid red dot at y = 0 on the x-axis.
        plot(data_params.trueV(i), 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
             'LineWidth', 2, 'DisplayName', 'True Value');
        
        % Add titles and labels.
        title(sprintf('%s: Empirical Density for Coordinate %d', algorithmName, i), ...
              'FontSize', 18, 'FontWeight', 'bold');
        xlabel(sprintf('Coordinate %d Value', i), 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Density', 'FontSize', 16, 'FontWeight', 'bold');
        grid on;
        legend('show', 'FontSize', 14, 'Location', 'best');
        hold off;
    end
end