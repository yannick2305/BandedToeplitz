%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [July 2025]
    Description:  [Pseudospectrum convergence]
    --------------------------------------------------------------
%}

close all;
clear all;
clc;

fs = 18;

% --- Parameters ---
    DimT = 50;     % Dimension Toeplitz "Operator" 
    N = 20;        
    p = 1.4;        % Decay rate upwards
    q = 3.6;        % Decay rate downwards
    algebraic = 0;  

% --- Initialize vector to store residuals ---
    steps = 2: 1: N-1;
    residuals = zeros(length(steps), 1);

% --- Loop over different truncation sizes ---
    for idx = 1:length(steps)
        n = steps(idx);
    
        % --- Generate sequence a_k for n-banded ---
        col = zeros(n,1);
        row = zeros(1,n+1);
        
        col(1) = 1; 
    
        % Inverval for coefficients [ai, bi]
        ai = 1;
        bi = 1;
        
        % Downward decay (sub-diagonals)
        for k = 1:n
            r = ai + (bi-ai)*rand;
            if algebraic == 1
                col(k) = 1 / ((k+1)^q);
            else 
                col(k) = r * exp(-q*k);
            end
        end
        
        % Upward decay (super-diagonals)
        for k = 1:n
            r = ai + (bi-ai)*rand;
            if algebraic == 1
                row(k+1) = 1 / ((k+1)^p);
            else
                row(k+1) = r * exp(-p*k);
            end
        end
        row(1) = 1;
    
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

        % Compute all eigenvalues and eigenvectors
        [V, D] = eig(T); 
    
        selectEigval = 1;
    
        % Extract and sort eigenvalues
        [lambda_sorted, sort_idx] = sort(diag(D), 'ascend');
    
        % Sort the eigenvectors accordingly
        V = V(:, sort_idx);
    
        lambda = lambda_sorted(selectEigval);
        v = V(:, selectEigval);
        
        % Truncate matrix and eigenvector
        T_n = T(1:n+1, 1:n+1);
        v_n = v(1:n+1);
        v_n = v_n / norm(v_n); 
        
        % Compute eigenvalues of truncated matrix
        eigenvalues_T_n = eig(T_n);
        
        % Find closest eigenvalue in T_n to the full lambda
        [~, idx_closest] = min(abs(eigenvalues_T_n - lambda));
        lambda_n = eigenvalues_T_n(idx_closest);
    
        res = (T_n - lambda * eye(n+1)) * v_n;
        residual = norm(res);
        residuals(idx) = residual;
    end


% --- Overlay the residual with the trendline for the pseudospectrum ---
    n_vals = (1:N-1)';
    
    % --- Pseudospectrum estimate for algebraic off-diagonal decay ---
    
    if algebraic == 1
        % --- trendline for algebraic off-diagonal decay ---
        p = min(p,q);
        trend = zeta(p)  .* ( ( n_vals.^(-2*p) ) .* (n_vals.^p - 1) ) ./ (n_vals.^(p ./ n_vals) - 1);

        trend = 0.5 * trend;
    else 
        % --- trendline for exponential off-diagonal decay ---
        nonreciprocit = 0.5 * log ( row(2) / col(1) );
        trend = exp(- n_vals .* nonreciprocit);

        constant = residuals(end) / trend(end);
        trend = constant * trend;
    end
    
% --- Plot residuals ---
    figure;
    if algebraic == 1
        loglog(steps, residuals, '-o', 'LineWidth', 2);
    else
        semilogy(steps, residuals, '-o', 'LineWidth', 2);
    end
    hold on;
    plot(n_vals, trend, '--r', 'LineWidth', 2);
    
    xlabel('Matrix size $N$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\varepsilon_N$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);
    set(gcf, 'Position', [100, 100, 500, 300]);
    grid on;

