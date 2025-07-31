%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [July 2025]
    Description:  [Complex Band Structure complex valued Lambda]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    n = 4;          % Number of bands 
    p = 2.3;        % Decay rate upwards      
    q = 2.8;        % Decay rate downwards       
    algebraic = 0;  % 1 for algebraic off diagonal decay, 0 for exponential 
    fs = 18;

    % --- Add noise to coefficients ---
    ai = 1;         
    bi = 1.4;

% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    col(1) = 1;
    
    for k = 1:n
        r = ai + (bi-ai)*rand;
        if algebraic == 1
            col(k) = 1 / ((k+1)^q);
        else 
            col(k) = r * exp(-q*k);
        end
    end
    
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


% --- Polynomial symbol ---
    syms z
    Qz = sym(0);
    for j = -n:n
        Qz = Qz + a(j + n + 1) * z^(n - j);
    end
    Qcoeffs = sym2poly(Qz);

% --- Complex lambda grid ---

    NP = 400;
    re_vals = linspace(0.75,  1.3,  NP);
    im_vals = linspace(-0.07, 0.07, NP);
    [ReL, ImL] = meshgrid(re_vals, im_vals);
    min_beta = zeros(size(ReL));
    
    % --- Main loop ---
    for i = 1:length(re_vals)
        for j = 1:length(im_vals)
            lambda = ReL(j,i) + 1i*ImL(j,i);
            Pcoeffs = Qcoeffs;
            Pcoeffs(n+1) = Pcoeffs(n+1) - lambda;  % z^m term index is n+1
    
            rts = roots(Pcoeffs);
            rts = rts(abs(rts) > 1e-12);  % Avoid log(0)
    
            beta_vals = -log(abs(rts));
    
            % Find the beta with the smallest absolute value, but keep the sign
            [~, idx] = min(abs(beta_vals));
            min_beta(j,i) = -beta_vals(idx);  % Preserves the sign
        end
    end

% --- Evaluate Symbol function on torus ---
    N = 200;
    theta = linspace(0, 2*pi, N);
    z_torus = exp(1i * theta);
    
    Qz_vals = zeros(1, N);
    for j = 1:N
        z = z_torus(j);
        Qz = 0;
        for k = -n:n
            Qz = Qz + a(k + n + 1) * z^(n - k);
        end
        Qz_vals(j) = Qz / z^n;
    end

% --- Plot CBS heatmap and symbol contour ---
    figure;
    imagesc(re_vals, im_vals, min_beta);
    axis xy;
    xlabel('$Re(\lambda)$ ', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$Im(\lambda)$',  'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);
    cb = colorbar; 
    cb.Label.String = '$\beta$'; 
    cb.Label.Interpreter = 'latex'; 
    cb.Label.FontSize = fs+3; 
    hold on;
    plot(real(Qz_vals), imag(Qz_vals), 'k', 'LineWidth', 2);
    colormap jet;

