clear all;
close all;

fs = 18;

% --- Parameters ---
n = 4;      % Truncation size for a_k
p = 1.8;    % Decay rate upwards
q = 2.4;    % Decay rate downwards
lambda_range = [0.0, 1.5]; % Frequency range
num_lambda = 1000;
algebraic = 0;

% --- Generate sequence a_k ---
col = zeros(n,1);
row = zeros(1,n+1);

col(1) = 1; % Diagonal element

ai = 0.85;
b = 1.1;

for k = 1:n
    r = ai + (b - ai) * rand;
    if algebraic == 1
        col(k) = 1 / ((k + 1)^q);
    else
        col(k) = r * exp(-q * k);
    end
end

for k = 1:n
    r = ai + (b - ai) * rand;
    if algebraic == 1
        row(k + 1) = 1 / ((k + 1)^p);
    else
        row(k + 1) = r * exp(-p * k);
    end
end
row(1) = 1;

a = [col(end:-1:1)', row];
m = (length(a) - 1) / 2;
if mod(m, 1) ~= 0
    error('Input vector a must have odd length (for degrees -m to m)');
end

% --- Prepare polynomial coefficients vector ---
coeff_Q = zeros(1, 2*m + 1);
for j = -m:m
    coeff_Q(m - j + 1) = a(j + m + 1);
end


lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);
idx_lambda_m = m + 1;

% Preallocate for alpha and beta
coeffs_first = coeff_Q;
coeffs_first(idx_lambda_m) = coeffs_first(idx_lambda_m) - lambda_vals(1);
roots_first = roots(coeffs_first);
roots_first = roots_first(~(abs(roots_first) < 1e-12));
num_roots = length(roots_first);

alpha_mat = NaN(num_roots, num_lambda);
beta_mat = NaN(num_roots, num_lambda);

for k = 1:num_lambda
    lambda = lambda_vals(k);
    coeffs = coeff_Q;
    coeffs(idx_lambda_m) = coeffs(idx_lambda_m) - lambda;
    
    roots_z = roots(coeffs);
    roots_z = roots_z(~(abs(roots_z) < 1e-12));
    
    [~, order] = sort(angle(roots_z));
    roots_z = roots_z(order);
    
    if length(roots_z) == num_roots
        alpha_mat(:, k) = angle(roots_z);
        beta_mat(:, k) = -log(abs(roots_z));
    end
end

% --- Augment alpha_mat with -pi where alpha ≈ pi ---
tol_pi = 1e-8;

alpha_cols = cell(1, num_lambda);
beta_cols  = cell(1, num_lambda);
max_len = 0;

for k = 1:num_lambda
    alpha_col = alpha_mat(:, k);
    beta_col = beta_mat(:, k);

    idx_pi = abs(alpha_col - pi) < tol_pi;

    alpha_new = alpha_col;
    beta_new  = beta_col;

    if any(idx_pi)
        alpha_new = [alpha_new; -pi * ones(sum(idx_pi), 1)];
        beta_new  = [beta_new;  beta_col(idx_pi)];
    end

    alpha_cols{k} = alpha_new;
    beta_cols{k}  = beta_new;
    max_len = max(max_len, length(alpha_new));
end

% Pad all columns to same length with NaNs
alpha_mat_ext = NaN(max_len, num_lambda);
beta_mat_ext  = NaN(max_len, num_lambda);

for k = 1:num_lambda
    len = length(alpha_cols{k});
    alpha_mat_ext(1:len, k) = alpha_cols{k};
    beta_mat_ext(1:len,  k) = beta_cols{k};
end

% Replace originals
alpha_mat = alpha_mat_ext;
beta_mat = beta_mat_ext;
num_roots = size(alpha_mat, 1);


% --- Toeplitz matrix and decay rate computation ---
DimT = 15;

colT = zeros(DimT, 1);
colT(1) = 1;
colT(2:length(col)+1) = col;

rowT = zeros(1, DimT);
rowT(1:length(row)) = row;

T = toeplitz(colT, rowT);

[V, D] = eig(T);
eigenvalues = diag(D);
n = size(T, 1);
decay_rates = zeros(n, 1);

for i = 1:n
    v = V(:, i);
    v = v / norm(v);
    x = (1:n)';
    y = log(abs(v));
    valid = isfinite(y);
    p = polyfit(x(valid), y(valid), 1);
    decay_rates(i) = p(1);
end

% --- Plot: continuous lines with breakpoints on jumps ---
tol_x = 0.13;  % jump threshold in alpha/beta
tol_y = 1.9 * abs(lambda_vals(2) - lambda_vals(1));

tol_x_alpha = 0.09;  % jump threshold in alpha/beta
tol_y_alpha = 1.9 * abs(lambda_vals(2) - lambda_vals(1));


figure;
hold on;

for r = 1:num_roots
    % ---- Alpha (black solid) ----
    x = alpha_mat(r, :);
    y = lambda_vals;
    x_clean = x(1); y_clean = y(1);
    
    for i = 2:length(x)
        dx = abs(x(i) - x(i-1));
        dy = abs(y(i) - y(i-1));
        if dx < tol_x_alpha && dy < tol_y_alpha
            x_clean(end+1) = x(i); %#ok<SAGROW>
            y_clean(end+1) = y(i); %#ok<SAGROW>
        else
            x_clean(end+1:end+2) = [NaN, x(i)];
            y_clean(end+1:end+2) = [NaN, y(i)];
        end
    end
    plot(x_clean, y_clean, '-', 'Color', 'k', 'LineWidth', 3.5);
end

for r = 1:num_roots
    % ---- Beta (red dashed) ----
    x = beta_mat(r, :);
    y = lambda_vals;
    x_clean = x(1); y_clean = y(1);
    
    for i = 2:length(x)
        dx = abs(x(i) - x(i-1));
        dy = abs(y(i) - y(i-1));
        if dx < tol_x && dy < tol_y
            x_clean(end+1) = x(i); %#ok<SAGROW>
            y_clean(end+1) = y(i); %#ok<SAGROW>
        else
            x_clean(end+1:end+2) = [NaN, x(i)];
            y_clean(end+1:end+2) = [NaN, y(i)];
        end
    end
    plot(-x_clean, y_clean, '-', 'Color', 'r', 'LineWidth', 3.5);
end

% --- Decay rates (blue Xs) ---
scatter(decay_rates, real(eigenvalues), 100, ...
    'MarkerEdgeColor', 'b', 'Marker', 'x', 'LineWidth', 2);

xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', fs);
grid on;
xlim([-pi - 0.01, pi + 0.01]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs + 2);
hold off;
box on
set(gcf, 'Position', [100, 100, 500, 300]);  % [left, bottom, width, height]

