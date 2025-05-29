
close all;
clear all;
clc;


n = 1000;
a = 2;              % Main diagonal
b = -exp(-0.2);      % 1st off-diagonal
c = -exp(-0.9);
d = 0;      % 4th off-diagonal
e = 0;        % 5th off-diagonal

T = banded_toeplitz(n, a, b, c);

% Compute eigenvalues
eigvals = eig(T);

% Limit of the Spectrum

Lower_Spectrum = min(real(eigvals));
Upper_Spectrum = max(real(eigvals));



% Define the symbol function

%f_z =  a + 2*b*cos(alpha)*cosh(beta) + 2*c * cos(2*alpha)*cosh(2*beta);

alpha = linspace(-pi-0.1, pi+0.1, 1000);  
beta  = linspace(-5,      5,      1000);  
[Alpha, Beta] = meshgrid(alpha, beta);

% Imaginary part
Im_z = 2*b*sin(Alpha).*sinh(Beta) + 2*c *sin(2*Alpha).* sinh(2*Beta) + 2 * d *sin(3*Alpha).* sinh(3*Beta) + 2 * e *sin(4*Alpha).* sinh(4*Beta) ;


% Plot only the contour where Im_z = 0
figure;
contour(Alpha, Beta, Im_z, [0 0], 'LineWidth', 2, 'LineColor', 'r')
xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\beta$',  'Interpreter', 'latex', 'FontSize', 12);
title('Contour where Im(f(z)) = 0')
grid on


%%

% Range of lambda values
lambda_vals = linspace(-5, 10, 1000);

% Preallocate storage
alpha_all = [];
beta_all = [];
lambda_all = [];

% Loop over lambda (cannot avoid this because each gives different coefficients)
for k = 1:length(lambda_vals)
    lambda = lambda_vals(k);
    
    % Polynomial coefficients: c*z^4 + b*z^3 + (a - lambda)*z^2 + b*z + c = 0
    coeffs = [c, b, (a - lambda), b, c];
    
    % Find roots of the quartic polynomial
    z_roots = roots(coeffs);
    
    % Filter valid roots (non-zero, finite)
    z_roots = z_roots(~isnan(z_roots) & abs(z_roots) > 1e-10);
    
    % Compute alpha and beta
    alpha = angle(z_roots);         % arg(z)
    beta = -log(abs(z_roots));      % -log|z|
    
    % Store results
    alpha_all = [alpha_all; alpha];
    beta_all = [beta_all; beta];
    lambda_all = [lambda_all; repmat(lambda, size(alpha))];
end

% Get current x-limits (after plotting, or set manually)
x_min = min([alpha_all(:); beta_all(:)]);
x_max = max([alpha_all(:); beta_all(:)]);

% Start plot
figure;

% Shade background region (blue, semi-transparent)
fill([x_min, x_max, x_max, x_min], ...
     [Lower_Spectrum, Lower_Spectrum, Upper_Spectrum, Upper_Spectrum], ...
     [0.8 0.9 1], 'EdgeColor', 'none'); % Light blue
hold on;

% Plot alpha and beta vs lambda
scatter(alpha_all, lambda_all, 10, 'b', 'filled'); 
scatter(beta_all, lambda_all, 10, 'r', 'filled');

% Labels and formatting
xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('\lambda');
xlim([-pi, pi]);
legend('','\alpha', '\beta');
grid on;
uistack(findobj(gca,'Type','scatter'),'top'); % Ensure scatter is above shading


%%


function T = banded_toeplitz(n, a, b, c)

    col = [a; b; c; zeros(n-3, 1)];
    row = [a, b, c, zeros(1, n-3)];

    T = toeplitz(col, row);
end


function D = exponential_diagonal_matrix(n, gamma)
% Generates an n x n diagonal matrix where D_ii = exp(gamma * i)

    i = (1:n)';               % Column vector of indices
    d = exp(gamma * i);       % Compute diagonal entries
    D = diag(d);              % Create diagonal matrix
end


function decay_rate = analyze_column_decay(A, i)
% Analyzes and plots the decay of the i-th column of matrix A
% Returns the estimated exponential decay rate

    % Get the i-th column
    col = abs(A(:, i));
    n = length(col);
    x = (1:n)';

    % Remove zeros to avoid log issues
    nonzero_idx = col > 0;
    x_nz = x(nonzero_idx);
    col_nz = col(nonzero_idx);

    % Semilog plot
    figure;
    semilogy(x, col, 'bo-', 'LineWidth', 1.5);
    xlabel('Index');
    ylabel('|A_{ji}|');
    title(['Semilog plot of column ', num2str(i)]);
    grid on;

    % Fit exponential: log(|A_{ji}|) = log(C) - alpha * j
    p = polyfit(x_nz, log(col_nz), 1);
    alpha = -p(1); % Decay rate is the negative slope
    C = exp(p(2)); % Intercept gives C

    % Plot fitted curve
    hold on;
    fitted_curve = C * exp(-alpha * x);
    semilogy(x, fitted_curve, 'r--', 'LineWidth', 2);
    legend('Column Data', sprintf('Fit: C * exp(-%.4f * j)', alpha));

    % Output decay rate
    decay_rate = alpha;
    fprintf('Estimated decay rate (alpha): %.4f\n', alpha);
end
