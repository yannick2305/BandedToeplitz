%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2025]
    Description:  [Illustrate Open spectrum and symbol collapse]
    --------------------------------------------------------------
%}

clear all;
close all;
clc;

% --- Parameters ---
    n = 6;      % Truncation size for a_k
    p = 1;      % Upwards decay rate
    q = 6;      % Downwards decay rate
    fs = 18;
    DimT = 80;  % Dimension of finite Toeplitz matrix to simulate open limit
    video_filename = 'Open_limit_collapse.mp4';

% --- Generate the Toeplitz operator ---
    col = zeros(n,1);
    row = zeros(1,n);
    
    col(1) = 1; % Diagonal element
    
    for k = 2:n
        col(k) = 1 / ((k)^q); % Downward decay (sub-diagonals)
        row(k) = 1 / ((k)^p); % Upward decay   (super-diagonals)
    end

    row(1) = 1;

    col(n+1 : DimT) = 0;
    col_band = col(1 : n);
    
    row(n+1 : DimT) = 0;
    row_band = row(1 : n);
    
    a = [ col_band(end:-1:1)', row_band(2:end) ];
    
    % --- Generate finite Toeplitz matrix ---
    T = toeplitz(col, row);
    eigT = sort(eig(T));
    
    n = n-1;

 % --- Generate a movie for symbol collapse ---
    % --- Generate symbol function ---
    syms z
    Qz = sym(0);
    for j = -n:n
        Qz = Qz + a(j + n + 1) * z^(n - j);
    end
    
    Qcoeffs = sym2poly(Qz);

    % --- Define the unit circle ---
    N = 1000;
    theta = linspace(0, 2*pi, N);
    z_torus = exp(1i * theta);

    % --- Range of r values ---
    r_values = linspace(3.5, 5, 1000);

    % --- Prepare figure ---
    fig = figure('Color','w');
    ylim([-1.3*max(imag(eigT)), 1.3*max(imag(eigT))]);
    grid on
    xlabel('$\Re$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\Im$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);

    hold on;

    % --- Plot "open limit" ---
    eig_plot = plot(real(eigT), imag(eigT), 'b.', 'MarkerSize', 8, 'LineWidth', 1.5);
    % --- Plot the symbol function on the scaled torus ---
    curve_plot = plot(nan, nan, 'k', 'LineWidth', 2);
    title(sprintf('r = %.3f', r_values(1)), 'FontSize', fs, 'Interpreter', 'latex');

    legend('$\lim_{n\to\infty} \sigma(T_n(f))$', '$f(r\mathbf{T})$', 'Interpreter', ...
        'latex', 'FontSize', 16., 'Location','northeast');

    % --- Prepare video writer ---
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 10;  
    open(v);

    % --- Animation loop ---
    for r = r_values
        z_torus = r * exp(1i * theta);
        Qz_vals = zeros(1, N);
    
        for j = 1:N
            z = z_torus(j);
            Qz = 0;
            for k = -n:n
                Qz = Qz + a(k + n + 1) * z^(n - k);
            end
            Qz_vals(j) = Qz / z^n;
        end
    
        % --- Update plot ---
        set(curve_plot, 'XData', real(Qz_vals), 'YData', imag(Qz_vals));
        title(sprintf('r = %.3f', r), 'FontSize', fs, 'Interpreter', 'latex');
        drawnow;
    
        % --- Write frame to video ---
        frame = getframe(fig);
        writeVideo(v, frame);
    end
