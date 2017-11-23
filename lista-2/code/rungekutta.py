import numpy as np
import matplotlib.pyplot as plt

# Parâmetros iniciais da equação da bolha
Pl = 0.195091
Pv = 0.166391
el = 4.34802
ev = 1.79634

f = lambda t,X: np.array([X[1], (35*(Pv-Pl) - 7*X[1]**2*(5*Pl + 5*el + 4*Pv + 4*ev))/(7*X[0]*(5*Pl+5*el+Pv+ev))])

h = 0.2                                 # Tamanho do passo
d = 2                                   # Quantidade de equações no sistema
t = np.arange(0, 22.2, 0.2)             # Intervalo da simulação
x = np.zeros((d, len(t)))               # Pré-alocação da matriz de resultados
x[:, 0] = np.array([0.5, 0.146446])     # Chute inicial

# Encontra soluções do sistema de eq. diferenciais utilizando o método de Runge-Kutta de quarta ordem.
for i in range(len(t)-1):
    k1 = h * f(t[i], x[:, i])
    k2 = h * f(t[i] + h/2, x[:, i] + k1.T/2)
    k3 = h * f(t[i] + h/2, x[:, i] + k2.T/2)
    k4 = h * f(t[i] + h, x[:, i] + k3.T)

    x[:, i+1] = x[:, i] + (k1 + 2*k2 + 2*k3 + k4)/6

# Plota as variáveis analisadas
plt.plot(t, x[0, :], 'b', label='Raio')
plt.plot(t, x[1, :], 'r', label='Velocidade')
plt.grid()
plt.legend()
plt.xlabel('Tempo (s)');
plt.ylabel('Variável medida');
plt.show()