function var_estimates_subsampling = get_var_estimates_subsampling(data, n, d, ...
                                                                 alpha, ...
                                                                 data_params, ...
                                                                 m1, m2, B, ...
                                                                 verbose)
    % Compute the overall oja_vec using all n samples.
    oja_vec = get_oja_vec(data, n, d, alpha, data_params);
  
    % Preallocate an array (m1 x d) to store the mean squared residuals for each group.
    mean_residual_squared_all_batches = zeros(m1, d);
    
    % Loop over groups (i = 1,...,m1) and within each group over batches (j = 1,...,m2)
    for i = 1:m1
        for j = 1:m2
            % Compute the overall batch index.
            batch_idx = m2*(i-1) + (j-1);
            % Determine the start and end indices for this batch within the variance estimation data.
            batch_start = B * batch_idx + 1;
            batch_end = min(batch_start + B - 1, n);
            
            % Get the oja_vec from this batch.
            vhat_t = get_oja_vec(data(batch_start:batch_end, :), ...
                                 batch_end - batch_start + 1, d, alpha, data_params);
                             
            % Compute the residual vector (projecting out the component along central_vec).
            residual_t = vhat_t - (oja_vec' * vhat_t) * oja_vec;
            
            % Add the elementwise square of the residual.
            % If residual_t is a column vector (d x 1), then (residual_t.^2)' is a row vector (1 x d).
            mean_residual_squared_all_batches(i, :) = mean_residual_squared_all_batches(i, :) + (residual_t.^2)';
        end
        % Average over the m2 batches in this group.
        mean_residual_squared_all_batches(i, :) = mean_residual_squared_all_batches(i, :) / m2;
    end
    
    % For each coordinate, take the median over the m1 groups.
    if m1 == 1
        mean_residual_squared = mean_residual_squared_all_batches(1, :);
    else
        mean_residual_squared = median(mean_residual_squared_all_batches, 1);
    end
    
    % Adjust the mean residual squared error using the learning rates.
    eta_n = get_learning_rate(n, alpha, data_params.eigengap);
    eta_B = get_learning_rate(B, alpha, data_params.eigengap);
    mean_residual_squared = mean_residual_squared * (eta_n / eta_B);
    
    % Return the structure containing the variance estimates and the oja_vec.
    var_estimates_subsampling = struct('variance', mean_residual_squared, ...
                                       'oja_vec', oja_vec);
    
    if verbose
        fprintf("Subsampling Estimate : sin^2 error : %.5f\n", 1 - (oja_vec' * data_params.trueV)^2);
    end
    % keyboard;
end