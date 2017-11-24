clear
clc
format long

addpath('../../tools');

%% Parâmetros iniciais da equação da bolha
Pl = 0.195091;
Pv = 0.166391;
el = 4.34802;
ev = 1.79634;

h = 0.1;
x0 = [1 0.0001];
tmin = 0;

%% Calcula tempo de colapso
% Parâmetros da equação posta na forma R'(t)^2 = -b/a + kR(t)^(-2a)
a = 1.2386;
b = 0.0058149;
k = b/a;
t0 = 0;
tmax = sqrt(a)/sqrt(b) * beta((a+1)/(2*a), 0.5)/(2*a) * x0(1);

f = @(t,X) [X(2) (35*(Pv-Pl) - 7*X(2)^2*(5*Pl + 5*el + 4*Pv + 4*ev))/(7*X(1)*(5*Pl+5*el+Pv+ev))]; % [R;U]
x = rk4(f, h, x0, [tmin tmax]); % Resolve sistema f usando Runge-Kutta de quarta ordem

%% Plota as variáveis analisadas (método da lista 2)
hold on
plot(tmin:h:tmax, x(:, 1));
plot(tmin:h:tmax, x(:, 2));
xlabel('Tempo (s)');
ylabel('Variável medida');
grid on;
legend('Raio', 'Velocidade');
hold off

h = 0.1;
R = @(x) sqrt(a*k./(b + a^3*k^3*(x - t0).^2));
t = @(tau) chebyshev(@(x) (R(x)).^(2*a+1), 0, tau, 22);
tau = 0;
ts = [];
rs = [];
while 1
    rs(end+1) = R(tau);
    ts(end+1) = t(tau);
    tau = tau + h;
    if ts(end) > tmax
        break
    end
end


% R = @(tau) sqrt(a*k./(b + a^3*k^2*(tau - t0).^2));
% F21 = @(a,b,c,z) gamma(c)/(gamma(b)*gamma(c-b)) * double(int(x^(b-1)*(1-z*x)^(-a)*(1-x)^(-b+c-1), x, 0, 1));
% t = @(tau) F21(0.5, 1+1/(2*a), 3/2, 1-R(tau)^(-2*a)) * sqrt(R(tau)^(-2*a) - 1)/(sqrt(a)*sqrt(b));
% ts = [];
% rs = [];
% for i = 0:1000
%     rs(end+1) = R(tau);
%     ts(end+1) = t(tau);
% end