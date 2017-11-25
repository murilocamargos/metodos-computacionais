clear
clc
format long

addpath('../../tools');

%% Parâmetros iniciais da equação da bolha
Pl = 0.195091;      % Pressão no líquido
Pv = 0.166391;      % Pressão no vapor
el = 4.34802;       % Energia do líquido
ev = 1.79634;       % Energia do vapor

h = 0.1;            % Passo para o método de Runge-Kutta
x0 = [1 0.0001];    % Valores das condições iniciais
tmin = 0;           % Início do intervalo

%% Calcula tempo de colapso
% Parâmetros da equação posta na forma R'(t)^2 = -b/a + kR(t)^(-2a)
a = 1.2386;         % alpha
b = 0.0058149;      % beta
k = b/a;            % constante de integração
t0 = 0;             % tau_0
tc = sqrt(a)/sqrt(b) * beta((a+1)/(2*a), 0.5)/(2*a) * x0(1);

%% Sistema de equações diferenciais a serem resolvidas
f = @(t,X) [X(2) (35*(Pv-Pl) - 7*X(2)^2*(5*Pl + 5*el + 4*Pv + 4*ev))/(7*X(1)*(5*Pl+5*el+Pv+ev))]; % [R;U]

%% Resolve pelo método numérico (aula-2)
x = rk4(f, h, x0, [tmin tc]);

%% Resolve o mesmo problema de forma teórica com equações paramétricas
R = @(x) (a*k./(b + a^3*k^2*(x - t0).^2)) .^ (1/(2*a));
t = @(tau) chebyshev(@(x) (R(x)).^(2*a+1), 0, tau, 22);
taumax = 450;
tabela = zeros(taumax+1, 2);
tabela(:,2) = R(0:taumax);
for tau = 0:taumax
    tabela(tau+1,1) = t(tau);
end

%% Plota resultados numérico e analítico
hold on
plot(tmin:h:tc, x(:, 1), 'b');
plot(tabela(:, 1), tabela(:, 2), 'rx');
legend('Raio (Numérico)', 'Raio (Analítico)');
xlabel('Tempo');
ylabel('Variável medida');
grid on;
hold off;