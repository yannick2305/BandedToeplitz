%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [July 2025]
    Description:  [Contour Vanishing imaginary part of symbol]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    n = 20;     % Truncation size for a_k
    p = 1.8;    % Decay rate upwards
    q = 1.8;    % Decay rate downwards
    NP = 200;   % Resolution of grid
    fs = 18;    % Fontsize for annotation

    algebraic = 1; % 1 for algebraic off diagonal decay, 0 for exponential

% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    
    col(1) = 1; % Diagonal element

    % Inverval for coefficients [a, b]
    ai = 1;
    bi = 2;
    
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
    
    a = [ col(end:-1:1)', row];

% --- Generate symbol function ---
    syms z
    Qz = sym(0);
    for j = -n:n
        Qz = Qz + a(j + n + 1) * z^(j);
    end

    % Convert symbolic Qz to function handle
    Qz_fun = matlabFunction(Qz, 'Vars', z);

% --- Evaluate the symbol function ---
    alpha = linspace(-pi-0.15, pi+0.15, NP);
    beta  = linspace(-2, 2, NP);
    [Alpha, Beta] = meshgrid(alpha, beta);
    
    Z = exp(-Beta + 1i * Alpha);
    
    Qval = Qz_fun(Z); % Evaluate Q(z) on the grid
    ImQ = imag(Qval);

% --- Plot contour for imag(Q) = 0 ---
    figure;
    contour(Alpha, Beta, ImQ, [0 0], 'LineWidth', 2, 'LineColor', 'k');
    xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\beta$',  'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);
    set(gcf, 'Position', [100, 100, 500, 300]);
    grid off;
