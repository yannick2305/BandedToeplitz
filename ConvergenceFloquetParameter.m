%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [July 2025]
    Description:  [Convergence of complex Floquet parameter]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    n = 10;         % Truncation size for a_k (Bandsize)
    p = 1.6;        % Decay rate upwards
    q = 6.8;        % Decay rate downwards (take q > p)
    lambda = 1;     % Set frequency (in the Bulk)
    algebraic = 0;  % Set 1 for algebraic off-diag decay or 0 for exponential
    fs = 16;

% --- Loop over truncation size --- 
    col = zeros(n,1);
    row = zeros(1,n+1);
        
    col(1) = 1;
        
    % --- Add noise to coefficients ---
    ai = 1;
    bi = 1.8;
        
    % Downward decay (sub-diagonals)
    for k = 1:n
        r = ai + (bi-ai)*rand;
        if algebraic == 1
            col(k) = r / ((k+1)^q);    % Algebraic off diagonal decay
        else
            col(k) = r * exp(-q*k);    % Exponential off diagonal decay
        end
    end
        
    % Upward decay (super-diagonals)
    for k = 1:n
        r = ai + (bi-ai)*rand;
        if algebraic == 1
            row(k+1) = r / ((k+1)^p);   % Algebraic off diagonal decay
        else
            row(k+1) = r * exp(-p*k);   % Exponential off diagonal decay
        end
    end
    row(1) = 1;

    median_moduli = zeros(1, n-1); 

% --- Compute the complex band structure --- 
    for kt = 1:n
    
        col_band = col(1:kt);
        row_band = row(1:kt+1);
        
        a = [ col_band(end:-1:1)', row_band];
        coeffs = a;
    
        % Subtract lambda * z^k
        coeffs(kt + 1) = coeffs(kt + 1) - lambda;
    
        % Compute roots
        z_roots = roots(coeffs);
        
        % Sort by modulus and pick median (closest to origin)
        moduli = abs(z_roots);
        moduli = sort(moduli);
        median_moduli(kt) = moduli(ceil(length(moduli)/2)+1);
    end
    
    n_vals = (1:n)'; 

% --- Compute the expected convergence ---
    if algebraic == 1
        trend = (log(n_vals) .* p) ./ n_vals;
        Constant = log(abs(median_moduli(n))) / trend(n);
    else 
        beta = log(abs(median_moduli(n)));
        trend = exp( -(p + beta) .* n_vals);
        Constant = abs(median_moduli(n-2) - median_moduli(n)) / trend(n-2);
    end

    trend = Constant * trend;

% --- Plotting the expected vs the numerical convergence ---
    figure;
    if algebraic == 1
        plot(1:n, log(abs(median_moduli)) , 'o-', 'LineWidth', 2);
    else
        semilogy(1:n, (abs(median_moduli - median_moduli(n)) ), 'o-', 'LineWidth', 2);
    end
    hold on;
    plot(n_vals, trend, '--r', 'LineWidth', 2);
    xlabel('Bandwidth $m$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$|z_0- z_{0,m}|$',    'Interpreter', 'latex', 'FontSize', fs);
    grid on;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);
    set(gcf, 'Position', [100, 100, 500, 300]); 
   
