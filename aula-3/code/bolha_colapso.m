clear
clc
format long

addpath('../../tools');

%% Par�metros iniciais da equa��o da bolha
Pl = 0.195091;      % Press�o no l�quido
Pv = 0.166391;      % Press�o no vapor
el = 4.34802;       % Energia do l�quido
ev = 1.79634;       % Energia do vapor

h = 0.1;            % Passo para o m�todo de Runge-Kutta
x0 = [1 0.0001];    % Valores das condi��es iniciais
tmin = 0;           % In�cio do intervalo

%% Calcula tempo de colapso
% Par�metros da equa��o posta na forma R'(t)^2 = -b/a + kR(t)^(-2a)
a = 1.2386;         % alpha
b = 0.0058149;      % beta
k = b/a;            % constante de integra��o
t0 = 0;             % tau_0
tc = sqrt(a)/sqrt(b) * beta((a+1)/(2*a), 0.5)/(2*a) * x0(1);

%% Sistema de equa��es diferenciais a serem resolvidas
f = @(t,X) [X(2) (35*(Pv-Pl) - 7*X(2)^2*(5*Pl + 5*el + 4*Pv + 4*ev))/(7*X(1)*(5*Pl+5*el+Pv+ev))]; % [R;U]

%% Resolve pelo m�todo num�rico (aula-2)
x = rk4(f, h, x0, [tmin tc]);

%% Resolve o mesmo problema de forma te�rica com equa��es param�tricas
R = @(x) (a*k./(b + a^3*k^2*(x - t0).^2)) .^ (1/(2*a));
t = @(tau) chebyshev(@(x) (R(x)).^(2*a+1), 0, tau, 22);
taumax = 450;
tabela = zeros(taumax+1, 2);
tabela(:,2) = R(0:taumax);
for tau = 0:taumax
    tabela(tau+1,1) = t(tau);
end

%% Plota resultados num�rico e anal�tico
hold on
plot(tmin:h:tc, x(:, 1), 'b');
plot(tabela(:, 1), tabela(:, 2), 'rx');
legend('Raio (Num�rico)', 'Raio (Anal�tico)');
xlabel('Tempo');
ylabel('Vari�vel medida');
grid on;
hold off;