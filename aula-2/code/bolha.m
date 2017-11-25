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
tmax = 12.5;

f = @(t,X) [X(2) (35*(Pv-Pl) - 7*X(2)^2*(5*Pl + 5*el + 4*Pv + 4*ev))/(7*X(1)*(5*Pl+5*el+Pv+ev))]; % [R;U]
x = rk4(f, h, x0, [tmin tmax]); % Resolve sistema f usando Runge-Kutta de quarta ordem

%% Plota as variáveis analisadas
hold on
plot(tmin:h:tmax, x(:, 1));
plot(tmin:h:tmax, x(:, 2));
xlabel('Tempo');
ylabel('Variável medida');
grid on;
legend('Raio', 'Velocidade');
hold off