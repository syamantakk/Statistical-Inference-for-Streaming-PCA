function oja_vec = get_oja_vec(data, n, d, alpha, data_params)
    oja_vec = randn(d,1);
    eta = get_learning_rate(n, alpha, data_params.eigengap);
    % Run Oja's algorithm
    for t = 1:n
        x_t = data(t, :)';
        oja_vec = oja_vec + eta * (x_t' * oja_vec) * x_t;
        oja_vec = oja_vec / norm(oja_vec);
    end
end