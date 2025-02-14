function var_estimates_bootstrap = get_var_estimates_bootstrap(data, n, d, ...
                                                             alpha, ...
                                                             data_params, ...
                                                             b, verbose)
    init = randn(d,1);
    oja_vec = init;
    bootstrap_vectors = zeros(d, b);
    for j=1:b
        bootstrap_vectors(:,j) = init;
    end
    eta = get_learning_rate(n, alpha, data_params.eigengap);
    W = randn(n, b)/sqrt(2);
    for t=2:n
        x_t = data(t, :)';
        x_t_1 = data(t-1, :)';
        oja_vec = oja_vec + eta * (x_t' * oja_vec) * x_t;
        oja_vec = oja_vec / norm(oja_vec);
        for k=1:b
            h_k = (x_t' * bootstrap_vectors(:,k))*x_t;
            g_k = (x_t_1' * bootstrap_vectors(:,k))*x_t_1;
            bootstrap_vectors(:,k) = bootstrap_vectors(:,k) + eta*(h_k + W(t, k)*(h_k-g_k));
            bootstrap_vectors(:,k) = bootstrap_vectors(:,k)/norm(bootstrap_vectors(:,k));
        end
    end
    mean_residual_squared = zeros(d,1);
    for k=1:b
        mean_residual_squared = mean_residual_squared + (bootstrap_vectors(:,k) - oja_vec*sign(bootstrap_vectors(:,k)'*oja_vec)).^2;
    end
    mean_residual_squared = mean_residual_squared/b;
    var_estimates_bootstrap = struct('variance', mean_residual_squared, ...
                                     'oja_vec', oja_vec);
    if(verbose)
        fprintf("Bootstrap Estimate : sin^2 error : %.5f\n", 1 - (oja_vec'*data_params.trueV)^2);
    end
end