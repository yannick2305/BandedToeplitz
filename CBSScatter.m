%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [June 2025]
    Description:  [Complex Band Structure for banded Toeplitz]
    --------------------------------------------------------------
%}

clear all;
close all;

fs = 18;

% --- Parameters ---
    n = 4;      % Truncation size for a_k
    p = 3.8;    % Decay rate upwards
    q = 5.8;    % Decay rate downwards
    lambda_range = [0.5, 1.3]; % Frequency range
    algebraic = 1;

% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    
    col(1) = 1; % Diagonal element

    % Inverval for coefficients [a, b]
    ai = 0.85;
    b = 1.1;
    
    % Downward decay (sub-diagonals)
    for k = 1:n
        r = ai + (b-ai)*rand;
        if algebraic == 1
            col(k) = 1 / ((k+1)^q);
        else 
            col(k) = r * exp(-q*k);
        end
    end
    
    % Upward decay (super-diagonals)
    for k = 1:n
        r = ai + (b-ai)*rand;
        if algebraic == 1
            row(k+1) = 1 / ((k+1)^p);
        else
            row(k+1) = r * exp(-p*k);
        end
    end
    row(1) = 1;
    
    a = [ col(end:-1:1)', row];

    num_lambda = 500; 

    % Check input
    m = (length(a) - 1) / 2;
    if mod(m, 1) ~= 0
        error('Input vector a must have odd length (for degrees -m to m)');
    end

    % Symbolic setup for Q(z) = sum_{j = -m}^m a_j z^{m - j}
    syms z
    Qz = sym(0);
    for j = -m:m
        Qz = Qz + a(j + m + 1) * z^(m - j);
    end

    % Sweep over lambda values
    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);

    alpha_all = [];
    beta_all = [];
    lambda_all = [];

    for k = 1:length(lambda_vals)
        lambda = lambda_vals(k);
        Pz = Qz - lambda * z^m;
        roots_z = roots(sym2poly(Pz));
        roots_z = roots_z(~(abs(roots_z) < 1e-12)); % Avoid log(0)

        alpha = angle(roots_z);
        beta = -log(abs(roots_z));

        alpha_all = [alpha_all; alpha(:)];
        beta_all  = [beta_all; beta(:)];
        lambda_all = [lambda_all; lambda * ones(length(alpha), 1)];
    end

% --- Plot the complex band structure and compare with numerical beta ---
    figure;
    scatter(alpha_all, lambda_all, 10, 'k', 'filled'); % lambda vs alpha
    hold on;
    scatter(beta_all,  lambda_all, 10, 'r', 'filled'); % lambda vs beta

    % Plot formatting

    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', fs);
    grid on;

    xlim([-pi - 0.01, pi + 0.01]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);

% --- Add decay length to CBS ---

DimT = 10; 

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

% Eigendecomposition
[V, D] = eig(T);
eigenvalues = diag(D);
n = size(T, 1);

% Initialize decay rate storage
decay_rates = zeros(n, 1);

% Loop through eigenvectors
for i = 1:n
    v = V(:, i);
    v = v / norm(v);  % normalize
    
    % Compute log(abs(v)) over index
    x = (1:n)';
    y = log(abs(v));  % component-wise decay
    
    % Only use components where v ≠ 0 to avoid -Inf
    valid = isfinite(y);
    x = x(valid);
    y = y(valid);
    
    % Fit line: log(|v_j|) ≈ -alpha * j + const
    p = polyfit(x, y, 1);
    alpha = p(1);  % decay rate
    
    decay_rates(i) = alpha;
end

% Plot decay rate vs eigenvalue
scatter(decay_rates, real(eigenvalues), 'x', 'SizeData', 100, 'MarkerEdgeColor', 'b','LineWidth', 2);   
hold off;

