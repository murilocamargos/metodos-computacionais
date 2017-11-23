clear
clc
format long

%% Par�metros iniciais da equa��o da bolha
Pl = 0.195091;
Pv = 0.166391;
el = 4.34802;
ev = 1.79634;

f = @(t,X) [X(2); (35*(Pv-Pl) - 7*X(2)^2*(5*Pl + 5*el + 4*Pv + 4*ev))/(7*X(1)*(5*Pl+5*el+Pv+ev))]; % [R;U]

h = 0.3;                            % Tamanho do passo
d = 2;                              % Quantidade de equa��es no sistema
t = 0:h:12.5;                       % Intervalo da simula��o
x = zeros(d, length(t));            % Pr�-aloca��o da matriz de resultados
x(:,1) = [1; 0.0001];               % Chute inicial

%% Encontra solu��es do sistema de eq. diferenciais utilizando o m�todo de Runge-Kutta de quarta ordem.
for i = 1:length(t)-1
    k1 = h * f(t(i), x(:, i));
    k2 = h * f(t(i) + h/2, x(:, i) + k1/2);
    k3 = h * f(t(i) + h/2, x(:, i) + k2/2);
    k4 = h * f(t(i) + h, x(:, i) + k3);
    x(:, i+1) = x(:, i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

%% Plota as vari�veis analisadas
hold on
plot(t, x(1,:));
plot(t, x(2,:));
xlabel('Tempo (s)');
ylabel('Vari�vel medida');
grid on;
legend('Raio', 'Velocidade');
hold off