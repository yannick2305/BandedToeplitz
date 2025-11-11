%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2025]
    Description:  [Reality of open limit]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    n = 8;      % Truncation size for a_k
    p = 2.3;    % Decay rate upwards
    q = 6.5;    % Decay rate downwards
    DimT = 80;  % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 100; 
    
% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    
    % --- Add noise to the coeffiients ---
    ai = 1;
    bi = 1;

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
    col = [1, col'];
    col = col';

    col(n+1 : DimT) = 0;
    col_band = col(1 : n);

    row(n+1 : DimT) = 0;
    row_band = row(1 : n);
    
    T = toeplitz(col, row);
    eigT = sort(eig(T));

    % --- Define interval for the CBS ---
    lambda_range = [min(real(eigT)) - 0.2, max(real(eigT)) + 0.2];

    % --- Plot the open spectrum ---
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    grid on; 
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 250]); 
    ylim([-1.2*max(imag(eigT)), 1.2*max(imag(eigT))]);
    hold off;


%% --- Generate the Complex band structure ---

    m = (length(a) - 1) / 2;

    % --- Generate the symbol function ---
    syms z
    Qz = sym(0);
    for j = -m:m
        Qz = Qz + a(j + m + 1) * z^(m - j);
    end

    % --- Generate frequency range ---
    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);

    % --- Initiallise the CBS ---
    alpha_all  = [];
    beta_all   = [];
    lambda_all = [];

    % --- Compute derivative of symbol function for confluent roots ---
    Qz_derivative = diff(Qz*z^(-m), z);
    
    % --- Confluent roots z satisfy f'(z) = 0 ---
    roots_Qz_derivative = double(solve(Qz_derivative, z)); 
   
    % --- Plot the set Λ(f) ---
    figure;
    plot(real(roots_Qz_derivative), imag(roots_Qz_derivative), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    hold on;

    % Initialize arrays to collect roots across all lambda values
    all_positive_imag = [];
    all_negative_imag = [];

    % --- Loop through lambda values and collect roots ---
    for k = 1:length(lambda_vals)
        lambda = lambda_vals(k);
        Pz = Qz - lambda * z^m;
        roots_z = roots(sym2poly(Pz));
        
        % Sort roots by magnitude
        [~, sort_idx] = sort(abs(roots_z));
        roots_z_sorted = roots_z(sort_idx);
        
        % Keep the m-th and (m+1)-th roots
        candidate_roots = roots_z_sorted(m:m+1);
        
        % Check if they have approximately the same modulus
        tol = 1e-6;
        mag1 = abs(candidate_roots(1));
        mag2 = abs(candidate_roots(2));
        
        if abs(mag1 - mag2) < tol
            selected_roots = candidate_roots;
            
            % Collect positive and negative imaginary parts
            positive_imag = selected_roots(imag(selected_roots) > 0);
            negative_imag = selected_roots(imag(selected_roots) < 0);
            
            all_positive_imag = [all_positive_imag; positive_imag(:)];
            all_negative_imag = [all_negative_imag; negative_imag(:)];
        end
        
        roots_z = roots_z(~(abs(roots_z) < 1e-12));
        
        % Compute the derivative Q
        alpha = angle(roots_z);
        beta = -log(abs(roots_z));
        
        alpha_all = [alpha_all; alpha(:)];
        beta_all  = [beta_all; beta(:)];
        lambda_all = [lambda_all; lambda * ones(length(alpha), 1)];
    end

    % Now plot the collected roots with lines
    if ~isempty(all_positive_imag)
        % Sort by real part for smooth curve
        [~, idx_pos] = sort(real(all_positive_imag));
        all_positive_imag = all_positive_imag(idx_pos);
        plot(real(all_positive_imag), imag(all_positive_imag), 'k.', 'LineWidth', 1.5)
    end
    
    if ~isempty(all_negative_imag)
        [~, idx_neg] = sort(real(all_negative_imag));
        all_negative_imag = all_negative_imag(idx_neg);
        plot(real(all_negative_imag), imag(all_negative_imag), 'k.', 'LineWidth', 1.5)
    end

    ylim([-7, 7])
    xlim([-4.3, 4])
    legend('Confluent roots', '$f_m^{-1}(R)$', 'Interpreter', 'latex', 'FontSize', 16);
    %legend('$f_m^{-1}(R)$', 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('$\Re$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\Im$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    
    hold off;

%%  --- Plot the corresponding complex band structure ---
    figure; 
    hold on;
    scatter(alpha_all, lambda_all, 10, 'k', 'filled'); % lambda vs alpha
    scatter(beta_all,  lambda_all, 10, 'r', 'filled'); % lambda vs beta

    % Plot formatting
    xlabel('\alpha (blue), \beta (green)');
    ylabel('\lambda');
    title('\lambda as a function of \alpha and \beta');
    grid on;


%%

   
analyze_laurent_roots(a, lambda_range, 1000);
%
%analyze_laurent_roots_hp(a, lambda_range, 1000, 80);  % 80-digit precision

lambda_guess = 1.2;
tol = 1e-12;         % convergence tolerance
maxIter = 50;        % max Newton iterations

%% Convert Laurent coefficients to standard polynomial coefficients
m = (length(a) - 1)/2;
if mod(m,1) ~= 0
    error('Length of a must be odd');
end
Q_coeffs = a;  % assume a is ordered a_{-m}..a_m

%% Define beta function numerically
beta_func = @(lambda) compute_beta(lambda, Q_coeffs, m);

%% Newton-Raphson iteration
lambda = lambda_guess;
for iter = 1:maxIter
    F = beta_func(lambda);
    h = 1e-8;  % finite difference step
    dF = (beta_func(lambda+h) - beta_func(lambda-h)) / (2*h);
    lambda_new = lambda - F/dF;

    if abs(lambda_new - lambda) < tol
        fprintf('Converged to lambda = %.12f after %d iterations\n', lambda_new, iter);
        lambda_zero_beta = lambda_new;
        break
    end
    lambda = lambda_new;
end

if iter == maxIter
    warning('Newton-Raphson did not converge within maxIter iterations');
    lambda_zero_beta = lambda;
end

disp('Lambda for which beta ≈ 0:');
disp(lambda_zero_beta);

%% Helper function (nested)
function F = compute_beta(lambda, Q_coeffs, m)
    % Construct P(z) = Q(z) - lambda*z^m
    P_coeffs = Q_coeffs;
    P_coeffs(m+1) = P_coeffs(m+1) - lambda;  % subtract lambda*z^m
    roots_z = roots(P_coeffs);
    roots_z = roots_z(abs(roots_z) > 1e-12);  % remove tiny roots
    [~, idx] = min(abs(abs(roots_z) - 1));   % root closest to |z|=1
    z_root = roots_z(idx);
    F = -log(abs(z_root));  % beta
end


%% --- Defining functions ---

function analyze_laurent_roots(a, lambda_range, num_lambda)
    % a: vector of coefficients a_{-m} to a_m (length must be odd)
    % lambda_range: [min_lambda, max_lambda]
    % num_lambda: number of lambda values to evaluate

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

    disp('f(z)')
    disp(Qz)

    % Sweep over lambda values
    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);

    alpha_all = [];
    beta_all = [];
    lambda_all = [];

    tol = 1e-3; % tolerance for checking beta = 0
        Qz_derivative = diff(Qz*z^(-m), z);
        disp(Qz_derivative);

    % Sweep over lambda values
    %lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);

    roots_Qz_derivative = double(solve(Qz_derivative, z));  % Convert to double for plotting
   

figure;
plot(real(roots_Qz_derivative), imag(roots_Qz_derivative), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
hold on;

% Initialize arrays to collect roots across all lambda values
all_positive_imag = [];
all_negative_imag = [];

% Loop through lambda values and collect roots
for k = 1:length(lambda_vals)
    lambda = lambda_vals(k);
    Pz = Qz - lambda * z^m;
    roots_z = roots(sym2poly(Pz));
    
    % Sort roots by magnitude
    [~, sort_idx] = sort(abs(roots_z));
    roots_z_sorted = roots_z(sort_idx);
    
    % Keep the m-th and (m+1)-th roots
    candidate_roots = roots_z_sorted(m:m+1);
    
    % Check if they have approximately the same modulus
    tol = 1e-6;
    mag1 = abs(candidate_roots(1));
    mag2 = abs(candidate_roots(2));
    
    if abs(mag1 - mag2) < tol
        selected_roots = candidate_roots;
        
        % Collect positive and negative imaginary parts
        positive_imag = selected_roots(imag(selected_roots) > 0);
        negative_imag = selected_roots(imag(selected_roots) < 0);
        
        all_positive_imag = [all_positive_imag; positive_imag(:)];
        all_negative_imag = [all_negative_imag; negative_imag(:)];
    end
    
    roots_z = roots_z(~(abs(roots_z) < 1e-12));
    
    % Compute the derivative Q
    alpha = angle(roots_z);
    beta = -log(abs(roots_z));
    
    alpha_all = [alpha_all; alpha(:)];
    beta_all  = [beta_all; beta(:)];
    lambda_all = [lambda_all; lambda * ones(length(alpha), 1)];
end

% Now plot the collected roots with lines
if ~isempty(all_positive_imag)
    % Sort by real part for smooth curve
    [~, idx_pos] = sort(real(all_positive_imag));
    all_positive_imag = all_positive_imag(idx_pos);
    plot(real(all_positive_imag), imag(all_positive_imag), 'k.', 'LineWidth', 1.5)
end

if ~isempty(all_negative_imag)
    [~, idx_neg] = sort(real(all_negative_imag));
    all_negative_imag = all_negative_imag(idx_neg);
    plot(real(all_negative_imag), imag(all_negative_imag), 'k.', 'LineWidth', 1.5)
end

ylim([-7, 7])
xlim([-4.3, 4])
legend('Confluent roots', '$f_m^{-1}(R)$', 'Interpreter', 'latex', 'FontSize', 16);
%legend('$f_m^{-1}(R)$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('$\Re$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\Im$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 

hold off;

    % Plot generic roots: lambda vs alpha and beta
    figure; hold on;
    scatter(alpha_all, lambda_all, 10, 'k', 'filled'); % lambda vs alpha
    scatter(beta_all,  lambda_all, 10, 'r', 'filled'); % lambda vs beta

    % === Find special lambda values where roots coincide ===
    dQz = diff(Qz, z);
    eqn = m * Qz == z * dQz;
    zs_crit = double(vpasolve(eqn, z, 'random', true));
    zs_crit = zs_crit(~isnan(zs_crit) & zs_crit ~= 0); % Clean up

    lambdas_crit = double(subs(Qz / z^m, z, zs_crit));

    % Plot critical points in red and magenta
    for i = 1:length(lambdas_crit)
        lambda = lambdas_crit(i);
        Pz = Qz - lambda * z^m;
        roots_z = roots(sym2poly(Pz));
        roots_z = roots_z(~(abs(roots_z) < 1e-12));

        alpha = angle(roots_z);
        beta  = -log(abs(roots_z));

        %scatter(alpha, lambda * ones(size(alpha)), 50, 'r', 'filled'); % red = alpha
        %scatter(beta,  lambda * ones(size(beta)),  50, 'm', 'filled'); % magenta = beta
    end

    % Plot formatting
    xlabel('\alpha (blue), \beta (green)');
    ylabel('\lambda');
    title('\lambda as a function of \alpha and \beta');
    grid on;
end



%%

function analyze_laurent_roots_hp(a, lambda_range, num_lambda, prec)
    % a: vector of 7 coefficients (a_{-3} to a_3)
    % lambda_range: [min, max]
    % num_lambda: number of lambda samples
    % prec: number of digits for vpa precision

    if nargin < 4
        prec = 50;
    end
    if length(a) ~= 7
        error('Input vector a must have length 7 (a_{-3} to a_3)');
    end

    m = 3;
    digits(prec);  % set vpa precision
    syms z

    % Define symbolic Laurent polynomial Q(z)
    Qz = sym(0);
    for j = -m:m
        Qz = Qz + vpa(a(j + m + 1)) * z^(m - j);
    end

    % Sweep lambda range
    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);
    alpha_all = [];
    beta_all = [];
    lambda_all = [];

    for k = 1:length(lambda_vals)
        lambda = vpa(lambda_vals(k));
        Pz = Qz - lambda * z^m;
        try
            zs = vpasolve(Pz == 0, z);  % high-precision root solver
            zs = double(zs(~(abs(zs) < 1e-12))); % remove small/zero roots

            alpha = angle(zs);
            beta = -log(abs(zs));

            alpha_all = [alpha_all; alpha(:)];
            beta_all  = [beta_all; beta(:)];
            lambda_all = [lambda_all; double(lambda) * ones(length(zs),1)];

        catch
            warning('vpasolve failed for lambda = %.4f', double(lambda));
        end
    end

    % Plotting
    figure; hold on;
    scatter(alpha_all, lambda_all, 10, 'k', 'filled'); % lambda vs alpha
    scatter(beta_all,  lambda_all, 10, 'r', 'filled'); % lambda vs beta
    xlabel('\alpha (blue), \beta (green)');
    ylabel('\lambda');
    title('High-Precision: \lambda vs \alpha and \beta');
    grid on;
end


function analyze_laurent_roots_m2(a, lambda_range, num_lambda)
    % a: 5-element vector of coefficients a_{-2} to a_2
    % lambda_range: [lambda_min, lambda_max]
    % num_lambda: number of lambda samples

    if length(a) ~= 5
        error('Input vector a must have 5 elements (a_{-2} to a_2)');
    end

    m = 2;
    syms z

    % Build polynomial Q(z) = sum a_j z^{m - j}
    Qz = sym(0);
    for j = -m:m
        Qz = Qz + a(j + m + 1) * z^(m - j);
    end

    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);
    alpha_all = [];
    beta_all = [];
    lambda_all = [];

    for k = 1:num_lambda
        lambda = lambda_vals(k);
        Pz = Qz - lambda * z^m;

        % Convert to numeric coefficients and solve explicitly
        coeffsP = sym2poly(Pz);
        zs = roots(coeffsP);
        zs = zs(~(abs(zs) < 1e-12)); % avoid log(0)

        alpha = angle(zs);
        beta = -log(abs(zs));

        alpha_all = [alpha_all; alpha(:)];
        beta_all  = [beta_all; beta(:)];
        lambda_all = [lambda_all; lambda * ones(length(zs), 1)];
    end

    % Plotting
    figure; hold on;
    scatter(alpha_all, lambda_all, 20, 'k', 'filled'); % lambda vs alpha
    scatter(beta_all,  lambda_all, 20, 'r', 'filled'); % lambda vs beta
    xlabel('\alpha (blue), \beta (green)');
    ylabel('\lambda');
    title('\lambda vs \alpha and \beta (explicit roots, m = 2)');
    grid on;
end


