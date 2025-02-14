function data_params = generate_data(n, d, c, b)
%GENERATE_DATA Generates the covariance matrix and its spectral properties.
%
%   data_cov = GENERATE_DATA(n, d, rho, iid_gap, coord) returns a structure
%   containing the following fields:
%
%       n                - Number of data points (default: 5000)
%       d                - Dimensionality (default: 10)

%       c                - Parameter for the covariance computation (default: set to 1)
%       b                - Exponent parameter (default: set to 1)
%       Sigma_true       - The true covariance matrix
%       Sigma_true_sqrtm - The square root of Sigma_true (via sqrtm)
%       trueV            - The leading eigenvector of Sigma_true
%       v2               - The second eigenvector (from eigs)
%       eigengap         - The eigenvalue gap between the first two eigenvalues
%       Vp               - All remaining eigenvectors (sorted in descending order)
%
%   Example:
%       data_cov = generate_data();
%
%   See also: sqrtm, eig, eigs, meshgrid.

    % Set default parameters if not provided
    if nargin < 1 || isempty(n)
        n = 5000;
    end
    if nargin < 2 || isempty(d)
        d = 10;
    end
    if nargin < 3 || isempty(c)
        c = 1;
    end
    if nargin < 4 || isempty(b)
        b = 1;
    end

    %% Compute the average covariance matrix
    sigma = (5 * (1:d).^(-b))';
    [p, q] = meshgrid(1:d, 1:d);
    Sigma_true = exp(-abs(p - q) * c) .* (sigma * sigma');
    Sigma_true_sqrtm = real(sqrtm(Sigma_true));

    %% Compute the leading eigenvector of the covariance matrix
    % Use eigs to compute the top two eigenvectors
    [trueV, trueE] = eigs(Sigma_true, 2);
    v2    = trueV(:, 2);    % Second eigenvector (if needed)
    trueV = trueV(:, 1);     % Leading eigenvector
    if(trueV(1) < 0)
        trueV = -trueV;
    end
    eigengap = trueE(1,1) - trueE(2,2);  % Difference between the first two eigenvalues

    % Compute full eigendecomposition for additional eigen information
    [V, D] = eig(Sigma_true);
    eigenvalues = diag(D);
    [~, idx] = sort(eigenvalues, 'descend');
    V_sorted = V(:, idx);
    Vp = V_sorted(:, 2:end);  % All eigenvectors except the first

    %% Pack all the data into a structure
    data_params = struct('n', n, ...
                      'd', d, ...
                      'c', c, ...
                      'b', b, ...
                      'Sigma_true', Sigma_true, ...
                      'Sigma_true_sqrtm', Sigma_true_sqrtm, ...
                      'trueV', trueV, ...
                      'v2', v2, ...
                      'eigengap', eigengap, ...
                      'Vp', Vp);
end
