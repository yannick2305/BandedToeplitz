%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2025]
    Description:  [Complex Band Structure for banded Toeplitz]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    n = 4;                      % Truncation size for a_k
    p = 3.8;                    % Decay rate upwards
    q = 5.8;                    % Decay rate downwards
    lambda_range = [0.7, 1.2];  % Frequency range
    algebraic = 1;              % 1 for alg. decay and 0 for exp. decay
    num_lambda = 500;           % Discretisation for frequency range
    DimT = 10;                  % Dimension of finite Toeplitz matrix
    fs = 18;                    % Fontsize of plot annotation

% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    col(1) = 1; 
    row(1) = 1;
    
    % --- Add noise to the coefficients ---
    ai = 1.0;
    bi = 1.3;
    
    % --- Downward decay (sub-diagonals) ---
    for k = 1:n
        r = ai + (bi-ai)*rand;
        if algebraic == 1
            col(k) = 1 / ((k+1)^q);
        else 
            col(k) = r * exp(-q*k);
        end
    end
    
    % --- Upward decay (super-diagonals) ---
    for k = 1:n
        r = ai + (bi-ai)*rand;
        if algebraic == 1
            row(k+1) = 1 / ((k+1)^p);
        else
            row(k+1) = r * exp(-p*k);
        end
    end

    % --- Retrieve coefficient vector ---
    a = [ col(end:-1:1)', row];

    % --- Generate the symbol function ---
    m = (length(a) - 1) / 2;
    syms z
    Qz = sym(0);
    for j = -m:m
        Qz = Qz + a(j + m + 1) * z^(m - j);
    end

    % --- Precompute the Polynomials in vectorised fashion ---
    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);
    Pz_all = repmat(Qz, num_lambda, 1) - lambda_vals(:) * z^m;

    alpha_all  = zeros(2*n*num_lambda, 1);
    beta_all   = zeros(2*n*num_lambda, 1);
    lambda_all = zeros(2*n*num_lambda, 1);

    idx = 1; 
    
    for k = 1:length(lambda_vals)
        lambda = lambda_vals(k);
        Pz = Qz - lambda * z^m;
        roots_z = roots(sym2poly(Pz));
        roots_z = roots_z(~(abs(roots_z) < 1e-12)); % Avoid log(0)
        alpha = angle(roots_z);
        beta = -log(abs(roots_z));
        
        n_roots = length(alpha);
        alpha_all(idx:idx+n_roots-1)  = alpha(:);
        beta_all(idx:idx+n_roots-1)   = beta(:);
        lambda_all(idx:idx+n_roots-1) = lambda;
        idx = idx + n_roots;
    end
    
    % --- Trim to actual size ---
    alpha_all  = alpha_all(1:idx-1);
    beta_all   = beta_all(1:idx-1);
    lambda_all = lambda_all(1:idx-1);

% --- Plot the complex band structure and compare with numerical beta ---
    figure;
    scatter(alpha_all, lambda_all, 10, 'k', 'filled');
    hold on;
    scatter(beta_all,  lambda_all, 10, 'r', 'filled');
    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', fs);
    grid on;
    box on;
    xlim([-pi - 0.01, pi + 0.01]);
    ylim([lambda_range(1), lambda_range(2)]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);

% --- Add decay length of eigenvectors to CBS ---
    % --- Generate finite Toeplitz matrix ---
    colT = zeros(DimT, 1);
    colT(1) = 1;
    for ind = 1:length(col)
        
        colT(ind + 1) = col(ind);
    end
    
    rowT = zeros(1,DimT);
    for ind = 1:length(row)
        rowT(ind) = row(ind);
    end
    
    T = toeplitz(colT, rowT);
    [V, D] = eig(T);
    eigenvalues = diag(D);

    % --- Initialize decay rate storage ---
    decay_rates = zeros(DimT, 1);

    % --- Loop through eigenvectors and measure exponential decay rate ---
    for i = 1: DimT
        v = V(:, i);
        v = v / norm(v); 
        
        x = (1: DimT)';
        y = log(abs(v)); 
        
        valid = isfinite(y);
        x = x(valid);
        y = y(valid);
        
        % Fit line: log(|v_j|) ≈ -beta * j + const
        p = polyfit(x, y, 1);
        beta = p(1);
        
        decay_rates(i) = beta;
    end

    % --- Add decay lengths to CBS ---
    scatter(decay_rates, real(eigenvalues), 'x', 'SizeData', 100, 'MarkerEdgeColor', 'b','LineWidth', 2);   
    hold off;



