clear variables
seed = 42;
rng(seed)
% Default parameters
S0 = 1; 
K = 1.1;
mu = 0.05;
r = 0.05; v0 = 0.04;
kappa = 5; theta = 0.4; 
sigma = 0.3; rho = 0.7; 
T = 1; dt = 0.01;
N = T / dt; t = (0:N-1) * dt;


Ns = [2,10,50,100,500,1000,10000,50000];
time_s_t = zeros(length(Ns),1);
time_s_t_e = zeros(length(Ns),1);
time_s_t_mils = zeros(length(Ns),1);
time_s_t_mils_l = zeros(length(Ns),1);
time_s_t_fourier = zeros(length(Ns),1);
price_euler_s_t = zeros(length(Ns),1);
price_euler_s_t_e = zeros(length(Ns),1);
price_s_t_mils = zeros(length(Ns),1);
price_s_t_mils_l = zeros(length(Ns),1);
price_s_t_fourier = zeros(length(Ns),1);

%% s_t_(k+1) = s_t(k) * exp((mu - 0.5 * v_t(k)) * dt + diffusion)
for n = 1:length(Ns)
    tic
    payoffs = zeros(Ns(n), 1);
    for c = 1:Ns(n)
        dWt_v_corr = randn(N, 1); dWt_s = randn(N, 1);
        dWt_s_corr = rho * dWt_v_corr + sqrt(1 - rho^2) * dWt_s;
        v_t = zeros(N,1); v_t(1) = v0;
        for k = 1:N-1
            drift = kappa*(theta-v_t(k))*dt;
            diffusion = sigma*sqrt(v_t(k))*dWt_v_corr(k)*sqrt(dt);  % Scale by sqrt(dt) for correct time-stepping
            v_t(k+1) = v_t(k) + drift + diffusion;
            v_t(k+1) = max(0.0001, v_t(k+1));
        end
        sqrt_v = sqrt(v_t);

        s_t_e = zeros(N,1); s_t_e(1) = S0;

        for k = 1:N-1
            diffusion = sqrt_v(k)*dWt_s_corr(k)*sqrt(dt);
            s_t_e(k+1) = s_t_e(k) * exp((mu - 0.5 * v_t(k)) * dt + diffusion);
        end
        payoff = max(s_t_e(end) - K, 0);
        payoffs(c) = payoff ;
        
    end
    price_euler_s_t_e(n) = exp(-r*T)*mean(payoffs);
    time_s_t_e(n) = toc;
end

%% s_t(k+1) = s_t(k) + mu*s_t(k)*dt + sqrt_v(k)*sqrt(dt)*s_t(k)*dWt_s_corr(k);
for n = 1:length(Ns)
    tic
    payoffs = zeros(Ns(n), 1); 
    for c = 1:Ns(n)
        dWt_v_corr = randn(N, 1); dWt_s = randn(N, 1);
        dWt_s_corr = rho * dWt_v_corr + sqrt(1 - rho^2) * dWt_s;
        v_t = zeros(N,1); v_t(1) = v0;
        for k = 1:N-1
            drift = kappa*(theta-v_t(k))*dt;
            diffusion = sigma*sqrt(v_t(k))*dWt_v_corr(k)*sqrt(dt);  % Scale by sqrt(dt) for correct time-stepping
            v_t(k+1) = v_t(k) + drift + diffusion;
            v_t(k+1) = max(0.0001, v_t(k+1));
        end
        sqrt_v = sqrt(v_t);

        % dS_t = mu*S_t*dt + sqrt(v)*S_t*dW_s
        s_t = zeros(N,1); s_t(1) = S0;

        for k = 1:N-1
            drift =  mu*s_t(k)*dt;
            diffusion = sqrt_v(k)*sqrt(dt)*s_t(k)*dWt_s_corr(k);
            s_t(k+1) = s_t(k) + drift + diffusion;
        end
        payoff = max(s_t(end) - K, 0);
        payoffs(c) = payoff ;
    end
    price_euler_s_t(n) = exp(-r*T)*mean(payoffs);
    time_s_t(n) = toc;
end

%% Milstein for St
for n = 1:length(Ns)
    
    tic
    payoffs = zeros(Ns(n), 1); 
    for c = 1:Ns(n)
        dWt_v_corr = randn(N, 1); dWt_s = randn(N, 1);
        dWt_s_corr = rho * dWt_v_corr + sqrt(1 - rho^2) * dWt_s;
        v_t_mils = zeros(N,1); v_t_mils(1) = v0;
        for k = 1:N-1
            v_t_mils(k+1) = (sqrt(v_t_mils(k))+0.5*sigma*sqrt(dt)*dWt_v_corr(k))^2 + kappa*(theta-v_t_mils(k))*dt - 0.25*sigma^2*dt;
            v_t_mils(k+1) = max(0.0001, v_t_mils(k+1));
        end
        sqrt_v_mils = sqrt(v_t_mils);

        % dS_t = mu*S_t*dt + sqrt(v)*S_t*dW_s
        s_t_mils = zeros(N,1); s_t_mils(1) = S0;

        for k = 1:N-1
            s_t_mils(k+1) = s_t_mils(k) + mu*s_t_mils(k)*dt + sqrt_v_mils(k)*sqrt(dt)*s_t_mils(k)*dWt_s_corr(k)+0.25*s_t_mils(k)^2*dt*(dWt_s_corr(k)^2 - 1) ;
        end
        payoff = max(s_t_mils(end) - K, 0);
        payoffs(c) = payoff ;
    end
    price_s_t_mils(n) = exp(-r*T)*mean(payoffs);
    time_s_t_mils(n) = toc;
end

%% Milstein for log(St)
for n = 1:length(Ns)
    tic
    payoffs = zeros(Ns(n), 1); 
    for c = 1:Ns(n)
        dWt_v_corr = randn(N, 1); dWt_s = randn(N, 1);
        dWt_s_corr = rho * dWt_v_corr + sqrt(1 - rho^2) * dWt_s;
        v_t_mils_l = zeros(N,1); v_t_mils_l(1) = v0;
        for k = 1:N-1
            v_t_mils_l(k+1) = (sqrt(v_t_mils_l(k))+0.5*sigma*sqrt(dt)*dWt_v_corr(k))^2 + kappa*(theta-v_t_mils_l(k))*dt - 0.25*sigma^2*dt;
            v_t_mils_l(k+1) = max(0.0001, v_t_mils(k+1));
        end
        sqrt_v_mils_l = sqrt(v_t_mils_l);

        % dS_t = mu*S_t*dt + sqrt(v)*S_t*dW_s
        s_t_mils_l = zeros(N,1); s_t_mils_l(1) = S0;

        for k = 1:N-1
            s_t_mils_l(k+1) = s_t_mils_l(k)*exp((mu - 0.5*v_t_mils_l(k))*dt + sqrt_v_mils_l(k)*sqrt(dt)*dWt_s_corr(k)) ;
        end
        payoff = max(s_t_mils_l(end)-K, 0);
        payoffs(c) = payoff ;
    end
    price_s_t_mils_l(n) = exp(-r*T)*mean(payoffs);
    time_s_t_mils_l(n) = toc;
end

%% Fourier Solution of Heston Model






%% Plots & Tables


figure(1)
hold on
plot(Ns,time_s_t,"DisplayName","Euler S_t")
plot(Ns,time_s_t_e,"DisplayName","Euler ln(S_t)")
plot(Ns,time_s_t_mils,"DisplayName",'Misltein S_t')
plot(Ns,time_s_t_mils_l,"DisplayName",'Misltein ln(S_t)')
%plot(Ns,time_s_t_fourier,"DisplayName","Fourier S_t")
xscale("log")
legend show
hold off
T_time = table(Ns', time_s_t, time_s_t_e,time_s_t_mils,time_s_t_mils_l, ...
    'VariableNames', {'nsamples', 'Euler_S_t', 'Euler_ln(S_t)','Misltein_S_t','Misltein_ln(S_t)'});
disp(T_time)

T_price  = table(Ns',price_euler_s_t,price_euler_s_t_e, price_s_t_mils,price_s_t_mils_l,...
    'VariableNames', {'nsampls', 'Euler_S_t', 'Euler_ln(S_t)','Misltein_S_t','Misltein_ln(S_t)'});
disp(T_price)