function heston_simulation_with_slider()
    seed = 42;
    rng(seed);
    % Default parameters
    S0 = 1; 
    K = 1.1;
    mu = 0.25;
    r = 0.05; v0 = 0.04;
    kappa = 5; theta = 0.4; 
    sigma = 0.3; rho = 0.7; 
    T = 1; dt = 0.001;
    N = T / dt; t = (0:N-1) * dt;
    
    % Create the figure
    fig = figure('Name', 'Heston Model Simulation', 'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);

    % Panel for sliders (Increased height)
    sliderPanel = uipanel('Parent', fig, 'Title', 'Parameters', 'Position', [0.05 0.6 0.3 0.35]);
    autoPanel = uipanel('Parent', fig, 'Title', 'Auto Controls', 'Position',[0.05 0.2 0.3 0.35]);

    % Consistent positioning parameters
    slider_left = 0.1;
    slider_width = 0.7;
    slider_height = 0.15; % Reduced slider height
    text_width = 0.2;
    text_offset = 0.02;
    
    % Slider positions
    theta_slider_bottom = 0.8;
    kappa_slider_bottom = 0.6;
    sigma_slider_bottom = 0.4;
    mu_slider_bottom = 0.2;
    button_height = 0.12; 
    theta_min = 0.01; theta_max = 0.2;
    kappa_min = 0.01; kappa_max = 20;
    sigma_min = 0.01; sigma_max = 0.5;
    mu_min = -1; mu_max = 1;

    % Create slider controls and text labels
    theta_slider = uicontrol('Parent', sliderPanel, 'Style', 'slider', 'Min', theta_min, 'Max', theta_max, 'Value', theta, 'Units', 'normalized', 'Position', [slider_left theta_slider_bottom slider_width slider_height], 'Callback', @update_plot);
    theta_text = uicontrol('Parent', sliderPanel, 'Style', 'text', 'Units', 'normalized', 'Position', [slider_left + slider_width + text_offset theta_slider_bottom text_width slider_height], 'String', sprintf('θ = %0.2f', theta));

    kappa_slider = uicontrol('Parent', sliderPanel, 'Style', 'slider', 'Min', kappa_min, 'Max',kappa_max, 'Value', kappa, 'Units', 'normalized', 'Position', [slider_left kappa_slider_bottom slider_width slider_height], 'Callback', @update_plot);
    kappa_text = uicontrol('Parent', sliderPanel, 'Style', 'text', 'Units', 'normalized', 'Position', [slider_left + slider_width + text_offset kappa_slider_bottom text_width slider_height], 'String', sprintf('κ = %0.2f', kappa));

    sigma_slider = uicontrol('Parent', sliderPanel, 'Style', 'slider', 'Min', sigma_min, 'Max', sigma_max, 'Value', sigma, 'Units', 'normalized', 'Position', [slider_left sigma_slider_bottom slider_width slider_height], 'Callback', @update_plot);
    sigma_text = uicontrol('Parent', sliderPanel, 'Style', 'text', 'Units', 'normalized', 'Position', [slider_left + slider_width + text_offset sigma_slider_bottom text_width slider_height], 'String', sprintf('σ = %0.2f', sigma));

    mu_slider = uicontrol('Parent', sliderPanel, 'Style', 'slider', 'Min', mu_min, 'Max', mu_max, 'Value', mu, 'Units', 'normalized', 'Position', [slider_left mu_slider_bottom slider_width slider_height], 'Callback', @update_plot);
    mu_text = uicontrol('Parent', sliderPanel, 'Style', 'text', 'Units', 'normalized', 'Position', [slider_left + slider_width + text_offset mu_slider_bottom text_width slider_height], 'String', sprintf('μ = %0.2f', mu));
    
    seed_text = uicontrol('Parent', sliderPanel, 'Style', 'text', 'Units', 'normalized', 'Position', [text_offset 0 text_width slider_height], 'String', sprintf('rng(seed) = %d', seed));
    seed_button = uicontrol('Parent', autoPanel, 'Style', 'pushbutton', 'String', 'refresh seed', 'Units', 'normalized', 'Position', [slider_left, 0 , 0.3, button_height], 'Callback', @(src, event) refresh_seed());

    % Adjust button positions below sliders
    % Adjust button positions below sliders
    mu_auto_button = uicontrol('Parent', autoPanel, 'Style', 'pushbutton', 'String', 'Mu Auto', 'Units', 'normalized', 'Position', [slider_left, mu_slider_bottom, 0.3, button_height], 'Callback', @(src, event) toggle_automation('mu'));
    kappa_auto_button = uicontrol('Parent', autoPanel, 'Style', 'pushbutton', 'String', 'Kappa Auto', 'Units', 'normalized', 'Position', [slider_left, kappa_slider_bottom , 0.3, button_height], 'Callback', @(src, event) toggle_automation('kappa'));
    sigma_auto_button = uicontrol('Parent', autoPanel, 'Style', 'pushbutton', 'String', 'Sigma Auto', 'Units', 'normalized', 'Position', [slider_left, sigma_slider_bottom , 0.3, button_height], 'Callback', @(src, event) toggle_automation('sigma'));
    theta_auto_button = uicontrol('Parent', autoPanel, 'Style', 'pushbutton', 'String', 'Theta Auto', 'Units', 'normalized', 'Position', [slider_left , theta_slider_bottom , 0.3, button_height], 'Callback', @(src, event) toggle_automation('theta'));
    % Create plot axes
    ax1 = subplot(2, 1, 1, 'Parent', fig, 'Position', [0.4 0.55 0.55 0.4]);
    title_1 = '$\mathbf{v_t} \hspace{1cm} dv_t = \kappa(\theta-v_t)dt + \xi\sqrt{v_t}dW_{v_t}$';
    title(ax1,title_1,"Interpreter","latex");
    ax2 = subplot(2, 1, 2, 'Parent', fig, 'Position', [0.4 0.05 0.55 0.4]);
    title_2 = "$ \mathbf{S_t} \hspace{1cm} dS_t = \mu S_t dt + \sqrt{v_t} S_t dWs_t $ ";
    title(ax2,title_2,"Interpreter","latex");

    % Create automation flags
    automation = struct('mu', false, 'kappa', false, 'sigma', sigma, 'theta', false);
    max_val = struct('mu', mu_max, 'kappa', kappa_max, 'sigma',sigma_max, 'theta', theta_max);
    min_val =struct('mu', mu_min, 'kappa', kappa_min, 'sigma',sigma_min, 'theta', theta_min);
    increment = 0.05; % 5%
    update_plot(); % Call it once to create the initial plot

    % Nested function to update plot
    function update_plot(~, ~)
        dWt_v_corr = randn(N, 1); dWt_s = randn(N, 1);
        dWt_s_corr = rho * dWt_v_corr + sqrt(1 - rho^2) * dWt_s;
        kappa = get(kappa_slider, 'Value');
        sigma = get(sigma_slider, 'Value');
        mu = get(mu_slider,'Value');
        theta = get(theta_slider,'Value');

        sqrt_v = zeros(N,1); sqrt_v(1) = sqrt(v0);
        for k = 1:N-1
            drift = kappa*((theta-sqrt_v(k)^2)*dt/(2*sqrt_v(k)));
            diffusion = sigma*dWt_v_corr(k)*sqrt(dt);  % Scale by sqrt(dt) for correct time-stepping
            sqrt_v(k+1) = sqrt_v(k) + drift + diffusion;
            sqrt_v(k+1) = max(0.0001, sqrt_v(k+1));
        end
        v_t = sqrt_v.^2;

        % Simulating the stock price (s_t) given by heston model
        % dS_t = mu*S_t*dt + sqrt(v)*S_t*dW_s
        s_t = zeros(N,1); s_t(1) = S0;
        for k = 1:N-1
            drift = mu*s_t(k)*dt;
            diffusion = sqrt(v_t(k))*s_t(k)*dWt_s_corr(k)*sqrt(dt);
            s_t(k+1) = s_t(k) + drift + diffusion;
        end


        % Update plots
        cla(ax1); hold(ax1, 'on');
        plot(ax1, t, v_t, "DisplayName", "v_t");
        yline(ax1, theta, "-b", "LineWidth", 1.5, "DisplayName", "long term mean");
        yline(ax1, mean(v_t), "--r", "LineWidth", 2, "DisplayName", "Empirical mean");
        legend(ax1, 'show'); hold(ax1, 'off');

        cla(ax2); hold(ax2,'on');
        plot(ax2, t, s_t, "DisplayName", "s_t");
        plot(ax2, t, S0 + t*mu, "-b", "LineWidth", 2, "DisplayName", "mu");
        xlim([0,max(t)])
        yline(ax2, mean(s_t), "--r", "LineWidth", 2, "DisplayName", "Empirical Mean");
        legend(ax2, 'show'); hold(ax2, 'off');

        % Update text labels
        theta_text.String = sprintf('θ = %0.2f', theta);
        kappa_text.String = sprintf('κ = %0.2f', kappa);
        sigma_text.String = sprintf('σ = %0.2f', sigma);
        mu_text.String = sprintf('μ = %0.2f', mu);
        seed_text.String = sprintf('rng(seed) = %d', seed);

        drawnow;
    end

    % Toggle automation on/off for each parameter
    function toggle_automation(variable)
        increment_val = (max_val.(variable))*(increment);
        automation.(variable) = ~automation.(variable);
        if automation.(variable)
            % Start automation
            while automation.(variable)
                % Increment the variable value
                current_value = get(eval([variable, '_slider']), 'Value');
                if current_value >= max_val.(variable)
                    current_value = min_val.(variable);
                end
                set(eval([variable, '_slider']), 'Value', current_value + increment_val);
                update_plot();
                pause(0.05); % Pause to slow down the automation
            end
        end
    end
    function refresh_seed()
        seed = randi(10000);
        rng(seed)
        update_plot()
    end
end
