import numpy as np
import matplotlib.pyplot as plt	

step = int(input("step: "))
lines = int(input("lines?: "))

# Загрузка данных из файлов
e = np.loadtxt(f'./out/e_{step}.txt')
t = np.loadtxt(f'./out/t_{step}.txt')
p = np.loadtxt(f'./out/p_{step}.txt')
d = np.loadtxt(f'./out/d_{step}.txt')
x = np.loadtxt(f'./out/k_x_{step}.txt')
y = np.loadtxt(f'./out/k_y_{step}.txt')
v_x = np.loadtxt(f'./out/v_x_{step}.txt')
v_y = np.loadtxt(f'./out/v_y_{step}.txt')

# Создание графиков
fig, axs = plt.subplots(2, 3, figsize=(10, 10))

line_1_1 = [1.1, -0.2]
line_2_1 = [0, 0]
line_1_2 = [1.1, -0.2]
line_2_2 = [0.9, 0.9]
line_1_3 = [0, 0]
line_2_3 = [1.1, -0.2]
line_1_4 = [0.9, 0.9]
line_2_4 = [1.1, -0.2]

# График для значения v
axs[0, 0].scatter(x, y, c = e, cmap='spring')
if lines == 1:
    axs[0, 0].plot(line_1_1, line_2_1, color = 'red')
    axs[0, 0].plot(line_1_2, line_2_2, color = 'red')
    axs[0, 0].plot(line_1_3, line_2_3, color = 'red')
    axs[0, 0].plot(line_1_4, line_2_4, color = 'red')
axs[0, 0].set_title('Энергия')

# График для значения t
axs[0, 1].scatter(x, y, c = t, cmap='spring')
if lines == 1:
	axs[0, 1].plot(line_1_1, line_2_1, color = 'red')
	axs[0, 1].plot(line_1_2, line_2_2, color = 'red')
	axs[0, 1].plot(line_1_3, line_2_3, color = 'red')
	axs[0, 1].plot(line_1_4, line_2_4, color = 'red')
axs[0, 1].set_title('Температура')

# График для значения p
axs[1, 0].scatter(x, y, c = p, cmap='spring')
if lines == 1:
	axs[1, 0].plot(line_1_1, line_2_1, color = 'red')
	axs[1, 0].plot(line_1_2, line_2_2, color = 'red')
	axs[1, 0].plot(line_1_3, line_2_3, color = 'red')
	axs[1, 0].plot(line_1_4, line_2_4, color = 'red')
axs[1, 0].set_title('Плотность')

# График для значения d
axs[1, 1].scatter(x, y, c = d, cmap='spring')
if lines == 1:
	axs[1, 1].plot(line_1_1, line_2_1, color = 'red')
	axs[1, 1].plot(line_1_2, line_2_2, color = 'red')
	axs[1, 1].plot(line_1_3, line_2_3, color = 'red')
	axs[1, 1].plot(line_1_4, line_2_4, color = 'red')
axs[1, 1].set_title('Давление')

axs[0, 2].scatter(x, y, c = v_x, cmap='spring')
if lines == 1:
	axs[0, 2].plot(line_1_1, line_2_1, color = 'red')
	axs[0, 2].plot(line_1_2, line_2_2, color = 'red')
	axs[0, 2].plot(line_1_3, line_2_3, color = 'red')
	axs[0, 2].plot(line_1_4, line_2_4, color = 'red')
axs[0, 2].set_title('Скорость X')

axs[1, 2].scatter(x, y, c = v_y, cmap='spring')
if lines == 1:
	axs[1, 2].plot(line_1_1, line_2_1, color = 'red')
	axs[1, 2].plot(line_1_2, line_2_2, color = 'red')
	axs[1, 2].plot(line_1_3, line_2_3, color = 'red')
	axs[1, 2].plot(line_1_4, line_2_4, color = 'red')
axs[1, 2].set_title('Скорость Y')

# Установка общего заголовка для всех графиков
fig.suptitle('Значения частиц')

# Отображение графиков
plt.show()