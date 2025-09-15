%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [September 2025]
    Description:  [Comparison: Jaffard vs Damko vs CBS]
    --------------------------------------------------------------
%}

close all;
clear all;

% --- Parameters ---
    n = 30;      % Matrix size
    p = 3.8;     % Decay rate upwards
    q = 3.8;     % Decay rate downwards
    bandw = 20;  % Bandwidth of the Toeplitz matrix

    fs = 18;        % Fontsize for plot annotations
    lambda = 2.4;   % Defect eigenfrequency

% --- Generate the m-banded Toeplitz Matrix ---

    col = zeros(n,1);
    row = zeros(1,n);
    
    col(1) = 1; % Diagonal element
    
    for k = 2:n
        col(k) = 1 / ((k)^q); % Downward decay (sub-diagonals)
        row(k) = 1 / ((k)^p); % Upward decay (super-diagonals)
    end

    col(bandw+1 : n) = 0;
    col_band = col(1 : bandw);

    row(bandw+1 : n) = 0;
    row_band = row(1 : bandw);
    
    T = toeplitz(col, row);
    
% --- Generate Green's function and defect eigenmode ---
    Green = (T - lambda * eye(n)) \ eye(n);

    indxv = 1;
    middle_column = Green(:, indxv);

% --- Compute Floquet Parameter ---
    a = [ col_band(end:-1:1)', row_band(2:end)];
    coeffs = a;
    
    % Subtract lambda * z^k
    coeffs(bandw) = coeffs(bandw) - lambda;
    
    % Compute roots
    z_roots = roots(coeffs);
        
    % Sort by modulus and pick median (closest to origin)
    moduli = abs(z_roots);
    moduli = sort(abs(moduli));
    median_moduli = min(abs(moduli-1));
    median_moduli = abs(median_moduli -1);

    beta = log(median_moduli);

% --- Add the predicted trendlines ---
    n_vals = (1:n)';

    kappa = cond(T); 
    qc = (sqrt(kappa) - 1) / (sqrt(kappa) + 1);
    
    % --- Complex band structure estimate ---
    trend_CBS = exp(beta .* (n_vals-1));
    Constant_CBS = abs(middle_column(bandw))/ trend_CBS(bandw);

    trend_CBS = Constant_CBS * trend_CBS;

    % --- Demko estimate
    trend_Demko = exp( 2/bandw .* log(qc) .* (n_vals-1) ) ;


    % --- Jaffard estimate ---
    trend_Jaffard = abs(middle_column(1)) .* (n_vals).^( -max(p, q) );
    Constant_Jaffard = abs(middle_column(bandw))/ trend_Jaffard(bandw);

    trend_Jaffard = Constant_Jaffard * trend_Jaffard;

% --- Plot the results ---
    figure;
    semilogy( 1:n, abs(middle_column)+eps, 'k', 'LineWidth', 2.5);
    hold on;
    plot(n_vals, trend_Jaffard,  ':b',  'LineWidth', 2);
    plot(n_vals, trend_CBS,      '--r', 'LineWidth', 2);
    plot(n_vals, trend_Demko,    'xg',  'LineWidth', 2);
    xlabel('$|i-j|$',                 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\mathbf{A}^{-1}(i,j)$',  'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);
    legend('Eigenmode', 'Jaffard Estimate', 'Complex band structure', 'Demko Estimate', 'Interpreter', 'latex', 'Location', 'southwest');
    grid on;
    set(gcf, 'Position', [100, 100, 500, 300]);  
