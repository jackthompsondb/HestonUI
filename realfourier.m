a = -10;
b = 10;
dt = 0.001;
x = a:dt:b;
N = length(x);
L = 2*b;
xi = (-N/2:N/2-1)*(1/(N*dt));
alpha = 2;

%% Fourier Transform of Gaussian by analytical and numerical 
f_hat = fftshift(fft(ifftshift(exp(-alpha*(x.^2)))))*dt;
f_inv_hat = fftshift(ifft(ifftshift(sqrt(pi/alpha) * exp(-((pi * xi).^2) / alpha))))/(dt);

figure(1)
hold on
plot(x,exp(-alpha*(x.^2)),"--","DisplayName","e^(-ax^2) f(a)","LineWidth",2)
plot(xi, sqrt(pi/alpha) * exp(-((pi * xi).^2) / alpha),"--","DisplayName", "Fourier Transform F(a)","LineWidth",2)
plot(xi,real(f_hat),"DisplayName","fft F(n)")
plot(x,real(f_inv_hat),"DisplayName","ifft f(N)")
xlim([a,b])
legend()
hold off

%% Fourier Transform of a function numerical


function F = fourier(f,x,a,b,N)
    A0 = 1/(b-a) * integral(f,a,b);
    trig_sum = 0;
    for n = 1:N
        An = (2/(b-a)) * quadgk(@(x) f(x) .* cos((2*pi*n*x)/(b-a)),a,b);
        Bn = (2/(b-a)) * quadgk(@(x) f(x) .* sin((2*pi*n*x)/(b-a)),a,b);
        trig_sum = trig_sum + An*cos((2*pi*n*x)/(b - a)) + Bn*sin((2*pi*n*x)/(b - a));
    end
    F = A0 + trig_sum;
end

function Res = residuals(F,f)
       Res = f - F;
end

% Linear Function
a = 0; b = 1;
dt = 0.001;
x = a:dt:b;
f_line = @(x) 7*x+3;

Ns = [8,16,32,64,128];

figure(2)
subplot(2,2,1)
hold on
colors = [linspace(0.8, 0, length(Ns))', linspace(0.9, 0, length(Ns))', ones(length(Ns), 1)];
plot(x,f_line(x),"DisplayName","7x+3","Color","red","LineWidth",2)
for i = 1:length(Ns)
    Fx_line = fourier(f_line,x,a,b,Ns(i));
    plot(x,Fx_line,"DisplayName",sprintf("N=%d",Ns(i)),"Color",colors(i,:))
end
plot(x,residuals(Fx_line,f_line(x)),"DisplayName","Residuals","Color","green")
title("Linear")
hold off

% Parabolic Function
a = -10; b = 10;
dt = 0.001;
x = a:dt:b;
f_parab = @(x) x.^2;
Ns = [8,16,32,64,128];

subplot(2,2,2)
hold on
colors = [linspace(0.8, 0, length(Ns))', linspace(0.9, 0, length(Ns))', ones(length(Ns), 1)];
plot(x,f_parab(x),"DisplayName","x^2","Color","red","LineWidth",2)
for i = 1:length(Ns)
    F_parab = fourier(f_parab,x,a,b,Ns(i));
    plot(x,F_parab,"DisplayName",sprintf("N=%d",Ns(i)),"Color",colors(i,:))
end
plot(x,residuals(F_parab,f_parab(x)),"DisplayName","Residuals","Color","green")
title("Parabolic")
legend show
hold off

% Polynomial
a = -10; b = 10;
dt = 0.001;
x = a:dt:b;
f_quad = @(x) 8*x.^3 -9*x.^2+3*x;
Ns = [8,16,32,64,128];

subplot(2,2,3)
hold on
colors = [linspace(0.8, 0, length(Ns))', linspace(0.9, 0, length(Ns))', ones(length(Ns), 1)];
plot(x,f_quad(x),"DisplayName","x^2","Color","red","LineWidth",2)
for i = 1:length(Ns)
    Fx_poly = fourier(f_quad,x,a,b,Ns(i));
    plot(x,Fx_poly,"DisplayName",sprintf("N=%d",Ns(i)),"Color",colors(i,:))
end
plot(x,residuals(Fx_poly,f_quad(x)),"DisplayName","Residuals","Color","green")
title("Polynomial")
hold off


% Square Function
a = -10; b = 10;
dt = 0.001;
x = a:dt:b;
f_square = @(x) square(x);
Ns = [8,16,32,64,128];

subplot(2,2,4)
hold on
colors = [linspace(0.8, 0, length(Ns))', linspace(0.9, 0, length(Ns))', ones(length(Ns), 1)];
plot(x,f_square(x),"DisplayName","x^2","Color","red","LineWidth",2)
for i = 1:length(Ns)
    Fx_square = fourier(f_square,x,a,b,Ns(i));
    plot(x,Fx_square,"DisplayName",sprintf("N=%d",Ns(i)),"Color",colors(i,:))
end
plot(x,residuals(Fx_square,f_square(x)),"DisplayName","Residuals","Color","green")
title("Square Func")
hold off

