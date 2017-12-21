clear
clc
format long

addpath('../../tools');

%% Parâmetros iniciais da equação da bolha
Pl = 0.17078553413498362;       % Pressão no líquido
Pv = 0.1595106949435665;        % Pressão no vapor
el = 4.275051939160018;         % Energia do líquido
ev = 1.7699955011416335;        % Energia do vapor

h = 0.01;               % Passo para o método de Runge-Kutta
x0 = [0.1 0.1];        % Valores das condições iniciais
tmin = 0;               % Início do intervalo
tmax = 6.67;          % Final do intervalo

%% Integrais em termos das funções hipergeométricas
I1 = @(x) (1/5) * hg2F1(1, 5/2, 7/2, x.^2);
I2 = @(x) (1/1) * hg2F1(1, 1/4, 5/4, x.^2);
I3 = @(x) (2/7) * hg2F1(1, 7/2, 9/2, x.^2);
I4 = @(x) (2/5) * hg2F1(1, 5/4, 9/4, x.^2);

%% Sistema de equações diferenciais a serem resolvidas
% Lembrando que R'(t) = U(t) = X(2)
f = @(t,X) [X(2) -(X(2)^2*(I2(X(1)) * (Pl+el) + 4*I1(X(1))*ev) + Pl + Pv*(4*I1(X(1))*X(2)^2 - 1))/(X(1)*(I4(X(1))*(Pl + el)*X(2)^2 + I2(X(1))*(Pl + el) + (Pv + ev)*(I1(X(1)) + I3(X(1))*X(2)^2)))]; % [R;U]

%% Resolve o sistema acima pelo método numérico
x = rk4(f, h, x0, [tmin tmax]);

%% Plota resultados numérico e analítico
subplot(1,2,1);
plot(tmin:h:tmax, x(:, 1), 'b');
legend('Raio');
xlabel('Tempo');
grid on;
subplot(1,2,2);
plot(tmin:h:tmax, x(:, 2), 'r');
legend('Velocidade');
xlabel('Tempo');
grid on;