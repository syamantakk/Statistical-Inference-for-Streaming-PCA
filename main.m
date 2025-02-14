%% Data Parameters
n = 2000;  % Number of data samples
d = 2000;  % Data dimension
c = 0.001; % Parameter of covariance matrix
b = 0.2;     % Parameter of covariance matrix
data_params = generate_data(n, d, c, b);

%% Experiment Parameters
alpha = 5; % Learning rate parameter
num_experiments = 200; % Number of experiments
verbose = 0;

% Parameters for subsampling estimator (Our algorithm)
m1 = 3;  % Number of batches for median
m2 = floor(max(log(n), log(d)));  % Number of batches for mean
B  = floor(n/(m1*m2)); % Number of elements in each batch

% Parameters for bootstrap estimator (Lunde, Sarkar, Ward (NeurIPS 2021))
num_bootstrap_samples = 1;

% Recording experiment results (each row corresponds to one experiment)
variance_subsampling_estimate_results = zeros(num_experiments, d);
variance_bootstrap_estimate_results   = zeros(num_experiments, d);

oja_vec_sampling_estimate_results   = zeros(num_experiments, d);
oja_vec_subsampling_estimate_results = zeros(num_experiments, d);
oja_vec_bootstrap_estimate_results   = zeros(num_experiments, d);

% Arrays to record execution times
subsampling_times = zeros(num_experiments, 1);
bootstrap_times   = zeros(num_experiments, 1);

%% Variance Estimation Experiments
for exp_num = 1:num_experiments
    fprintf("Experiment %d\n", exp_num);

    % Pre-process the IID samples
    Z = sqrt(3) * (2 * unifrnd(0, 1, [n, d]) - 1);
    data = Z * data_params.Sigma_true_sqrtm';

    %%% Sampling estimator (baseline)
    sampling_oja_vec = get_oja_vec(data, n, d, alpha, data_params);
    % Align the sign with the true vector.
    sampling_oja_vec = sampling_oja_vec * sign(sampling_oja_vec' * data_params.trueV);
    oja_vec_sampling_estimate_results(exp_num, :) = sampling_oja_vec;

    %%% Subsampling estimator (Our Algorithm)
    tic;  % Start timing subsampling estimator
    var_estimates_subsampling = get_var_estimates_subsampling(data, n, d, ...
                                            alpha, data_params, ...
                                            m1, m2, B, verbose);
    subsampling_times(exp_num) = toc;  % Record elapsed time
    
    variance_subsampling_estimate_results(exp_num, :) = var_estimates_subsampling.variance;
    subsampling_oja_vec = var_estimates_subsampling.oja_vec;
    oja_vec_subsampling_estimate_results(exp_num, :) = subsampling_oja_vec*sign(subsampling_oja_vec'*data_params.trueV);

    %%% Bootstrap estimator
    tic;  % Start timing bootstrap estimator
    var_estimates_bootstrap = get_var_estimates_bootstrap(data, n, d, ...
                                            alpha, data_params, num_bootstrap_samples, verbose);
    bootstrap_times(exp_num) = toc;  % Record elapsed time
    
    variance_bootstrap_estimate_results(exp_num, :) = var_estimates_bootstrap.variance;
    bootstrap_oja_vec = var_estimates_bootstrap.oja_vec;
    oja_vec_bootstrap_estimate_results(exp_num, :) = bootstrap_oja_vec*sign(bootstrap_oja_vec'*data_params.trueV);
    
    fprintf("-------------\n");
end

%% Compute and print average times and standard deviations

% Subsampling estimator timings
avg_subsampling_time = mean(subsampling_times);
std_subsampling_time = std(subsampling_times);

% Bootstrap estimator timings
avg_bootstrap_time = mean(bootstrap_times);
std_bootstrap_time = std(bootstrap_times);

fprintf('Subsampling estimator average time: %.4f +/- %.4f seconds\n', ...
        avg_subsampling_time, std_subsampling_time);
fprintf('Bootstrap estimator average time:   %.4f +/- %.4f seconds\n', ...
        avg_bootstrap_time, std_bootstrap_time);

%% Compute "true" variances from the sampling estimator estimates
% (This is our plug-in variance for the sampling estimator, available only
%  as an overall (across experiments) value.)
true_variances = var(oja_vec_sampling_estimate_results, 1);

%% Coverage Calculation for the Top 10 Coordinates
% For each coordinate i, we construct a 95% confidence interval using
% the formula:
%
%   CI_i = [ oja_vec(i) - 1.96 * sqrt(variance(i)),
%            oja_vec(i) + 1.96 * sqrt(variance(i)) ]
%
% Then, we check whether the true coordinate (data_params.trueV(i)) is contained
% in that interval for each experiment.

z_val = 1.96;
top_coords = 1:10;  % Define the top 10 coordinates

% Preallocate coverage arrays (one value per coordinate)
coverage_sampling   = zeros(1, length(top_coords));
coverage_subsampling = zeros(1, length(top_coords));
coverage_bootstrap   = zeros(1, length(top_coords));

for idx = 1:length(top_coords)
    i = top_coords(idx);
    
    % --- Sampling estimator ---
    % For sampling, we do not have an individual variance estimate per experiment.
    % We use the overall empirical variance computed across experiments.
    lower_bounds_sampling = oja_vec_sampling_estimate_results(:, i) - z_val * sqrt(true_variances(i));
    upper_bounds_sampling = oja_vec_sampling_estimate_results(:, i) + z_val * sqrt(true_variances(i));
    coverage_sampling(idx) = mean((data_params.trueV(i) >= lower_bounds_sampling) & ...
                                  (data_params.trueV(i) <= upper_bounds_sampling));
    
    % --- Subsampling estimator ---
    lower_bounds_sub = oja_vec_subsampling_estimate_results(:, i) - z_val * sqrt(variance_subsampling_estimate_results(:, i));
    upper_bounds_sub = oja_vec_subsampling_estimate_results(:, i) + z_val * sqrt(variance_subsampling_estimate_results(:, i));
    coverage_subsampling(idx) = mean((data_params.trueV(i) >= lower_bounds_sub) & ...
                                     (data_params.trueV(i) <= upper_bounds_sub));
    
    % --- Bootstrap estimator ---
    lower_bounds_boot = oja_vec_bootstrap_estimate_results(:, i) - z_val * sqrt(variance_bootstrap_estimate_results(:, i));
    upper_bounds_boot = oja_vec_bootstrap_estimate_results(:, i) + z_val * sqrt(variance_bootstrap_estimate_results(:, i));
    coverage_bootstrap(idx) = mean((data_params.trueV(i) >= lower_bounds_boot) & ...
                                   (data_params.trueV(i) <= upper_bounds_boot));
end

%% Print coverage results
fprintf('\nCoverage for top 10 coordinates (Sampling estimator):\n');
for i = 1:length(top_coords)
    fprintf('  Coordinate %d: %.2f%%\n', top_coords(i), coverage_sampling(i)*100);
end

fprintf('\nCoverage for top 10 coordinates (Subsampling estimator):\n');
for i = 1:length(top_coords)
    fprintf('  Coordinate %d: %.2f%%\n', top_coords(i), coverage_subsampling(i)*100);
end

fprintf('\nCoverage for top 10 coordinates (Bootstrap estimator):\n');
for i = 1:length(top_coords)
    fprintf('  Coordinate %d: %.2f%%\n', top_coords(i), coverage_bootstrap(i)*100);
end

%% Plot results for all experiments

% Compute the subsampling mean and standard deviation.
subsampling_var_mean = mean(variance_subsampling_estimate_results, 1);
subsampling_var_std  = std(variance_subsampling_estimate_results, 1);

% Compute the bootstrap mean and standard deviation.
bootstrap_var_mean = mean(variance_bootstrap_estimate_results, 1);
bootstrap_var_std  = std(variance_bootstrap_estimate_results, 1);

%% Plot 1: Variances of interesting coordinates
% Define the coordinates to plot (e.g., the first 10 coordinates)
interesting_coordinates = 1:10;  % x-axis representing coordinate indices
generate_plots_variance_topk_coordinates(interesting_coordinates, true_variances, ...
                                         subsampling_var_mean, subsampling_var_std, ...
                                         bootstrap_var_mean, bootstrap_var_std);

%% Plot 2 : Histograms of interesting coordinates
% Plot additional histograms with Gaussian overlays for the top 3 coordinates
% Define the top 3 coordinates you wish to plot.
top_coords_for_hist = 1:3;

% Plot for the Subsampling estimator.
generate_histograms_for_algorithm('OjaVarEst', top_coords_for_hist, ...
    oja_vec_subsampling_estimate_results, variance_subsampling_estimate_results, data_params);

% Plot for the Bootstrap estimator.
generate_histograms_for_algorithm('Bootstrap', top_coords_for_hist, ...
    oja_vec_bootstrap_estimate_results, variance_bootstrap_estimate_results, data_params);