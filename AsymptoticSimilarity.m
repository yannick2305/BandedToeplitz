%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2025]
    Description:  [Asymptotic Similarity Transform for Toeplitz matrices]
    ----------------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 16;             % Truncation size for a_k
    p = 6.5;            % Decay rate upwards
    q = 6.5;            % Decay rate downwards
    DimT = 30;          % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 300;   % Number of plotting points (keep it large)
    fs = 18;            % Fontsize for annotation
    
% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    
    % --- Add noise to the coeffiients ---
    ai = 1;
    bi = 1.5;

    col(1) = 1; % Diagonal element
    row(1) = 1;

    % Downward decay (sub-diagonals)
    for k = 1:n
        r = ai + (bi-ai)*rand;
        col(k) = r / ((k+1)^q);

        r = ai + (bi-ai)*rand;
        row(k+1) = r / ((k+1)^p);
    end
   
    % --- Coefficients of the symbol function ---
    a = [ col(end:-1:1)', row];
       
    % --- Generate finite Toeplitz matrix ---
    T = fourier_to_toeplitz(a, DimT);
    eigT = sort(eig(T));

% --- Get approximation of asymptoric spectrum ---
    lambda_range = [min(real(eigT)) - 0.01, max(real(eigT)) + 0.01];

% --- Approximate the open limit ---
    % --- Generate frequency range ---
    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);

    % --- Preallocate for all polynomial coefficients ---
    % P(z) = Q(z) - lambda*z^m
    % We need to subtract lambda from the coefficient of z^m
    P_coeffs_matrix = repmat(a(:), 1, num_lambda);
    
    % The coefficient of z^m is at position (m+1) in the array
    % Subtract lambda values from that position
    P_coeffs_matrix(n + 1, :) = P_coeffs_matrix(n + 1, :) - lambda_vals;
    
    % --- Compute all roots at once ---
    all_roots = zeros(2*n, num_lambda);
    for k = 1:num_lambda
        all_roots(:, k) = roots(P_coeffs_matrix(:, k));
    end
    
    % --- Sort roots by magnitude for each lambda ---
    [~, sort_idx] = sort(abs(all_roots), 1);
    % Convert linear indices for sorted access
    linear_idx = sort_idx + (0:num_lambda-1) * (2*n);
    all_roots_sorted = all_roots(linear_idx);
    
    % --- Extract m-th and (m+1)-th roots ---
    candidate_roots = all_roots_sorted(n:n+1, :);
    
    % --- Check if they have approximately the same modulus ---
    tol = 1e-6;
    mag1 = abs(candidate_roots(1, :));
    mag2 = abs(candidate_roots(2, :));
    valid_mask = abs(mag1 - mag2) < tol;
    
    % --- Collect valid roots ---
    openLimit = candidate_roots(1, valid_mask);

    % --- Compute polar curve in polar coordinates ---
    open = [openLimit, conj(openLimit(end:-1:1))] ;
    openLimit = open;

% --- Plot the set Λ(f) ---
    Lambda_f = [openLimit, openLimit(1)];
    figure;
    plot(real(Lambda_f), imag(Lambda_f), 'k-', 'LineWidth', 2.5)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    grid on;
    box on;
    hold off;

% --- compute Fourier coefficients of f(p(z)) ---

    phase = zeros(length(openLimit),1);

    for i = 1:length(openLimit)
        phase(i) = angle(openLimit(i));
    end
    
    % --- Evaluate f(p(z)) on the torus ---
    N = length(openLimit);
    fp_values = zeros(N, 1);
    
    for j = 1:N
        for k = -n:n
            k_idx = k + n + 1; 
            fp_values(j) = fp_values(j) + a(k_idx) * openLimit(j)^(-k);
        end
    end

    % --- Compute the Fourier Transform of f(p(z)) ---
    FourierFP = fourier_nonuniform(phase, fp_values, DimT);

    % --- Plot the decay of the Fouier Coefficients of f(p(z)) ---
    figure;
    loglog(1:DimT, abs(FourierFP(DimT+1: 2*DimT)), 'k.', 'LineWidth', 1.5)
    hold on
    plot(row)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    grid on;
    hold off;

    % --- Toeplitz matrix for deformed path ---
    T_b = fourier_to_toeplitz(FourierFP, DimT);
    eigT_b = sort(eig(T_b));
    similar = norm(eigT - eigT_b);
    fprintf('Similarity Coefficient: %f\n', similar);

% --- Plot the eigenvalues before and after asymptotic Similarity transform ---
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    plot(real(eigT_b), imag(eigT_b), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    box on;
    xlabel('$\mathrm{Re}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 250]); 
    axis equal;
    hold off;


%% --- Defining functions ---


function [c, nrange] = fourier_nonuniform(alpha, g, M)
    % Compute Fourier coefficients using trapezoidal rule

    % --- Ensure column vectors ---
    alpha = alpha(:);
    g     = g(:);

    % --- Sort by angle to construct quadrature weights ---
    [alpha, idx] = sort(alpha);
    g = g(idx);

    % Periodic trapezoidal weights for nonuniform samples
    dalpha = diff([alpha; alpha(1)+2*pi]);
    dalpha_prev = [dalpha(end); dalpha(1:end-1)];
    w = 0.5*(dalpha + dalpha_prev);    % weights w_j

    % Fourier index range
    nrange = (-M:M).';
    L = length(nrange);

    % Build weighted design matrix: E(j,n) = exp(i n alpha_j)
    E = exp(1i * alpha * nrange.');

    % Weighted least squares setup
    Wsqrt = sqrt(w);               % sqrt of quadrature weights
    EW = E .* Wsqrt;               % apply weights row-wise
    gw = g .* Wsqrt;

    % Tikhonov regularization for stability
    lambda = 1e-10;

    % Solve LS system: (EW'*EW + λI)c = EW'*gw
    c = (EW' * EW + lambda*eye(L)) \ (EW' * gw);

end

function T = fourier_to_toeplitz(a, dimT)
    K = (length(a) - 1) / 2;
    a_0 = a(K+1);          
    
    col = zeros(dimT, 1);
    col(1) = a_0;
    row = zeros(1, dimT);
    row(1) = a_0;

    for k = 1:min(K, dimT-1)
        col(k+1) = a(K+1-k);
        row(k+1) = a(K+1+k);
    end
   
    T = toeplitz(col, row);
end
