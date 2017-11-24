clear
clc
format long

addpath('../../tools');

%% Parâmetros iniciais da equação da bolha
Pl = 0.195091;
Pv = 0.166391;
el = 4.34802;
ev = 1.79634;

h = 0.01;
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

%% Calculo dos valores das equações paramétricas
R = @(x) sqrt(a*k./(b + a^3*k^2*(x - t0).^2));
t = @(tau) chebyshev(@(x) (R(x)).^(2*a+1), 0, tau, 22);

tau = 0:1000;
ts = [];
rs = R(tau);
for tau = tau
    ts(end+1) = t(tau);
end
plot(ts,rs)


% F21 = @(a,b,c,z) gamma(c)/(gamma(b)*gamma(c-b)) * chebyshev(@(x) x.^(b-1).*(1-z*x).^(-a).*(1-x).^(-b+c-1), eps, 1-eps, 22);
% T = @(R) (sqrt(R^(-2*a) - 1) * F21(0.5, 1+1/(2*a), 1.5, 1-R^(-2*a)))/(sqrt(a)*sqrt(b));
% ts = [];
% for i=eps:0.01:1-eps
%     ts(end+1) = T(i);
% end