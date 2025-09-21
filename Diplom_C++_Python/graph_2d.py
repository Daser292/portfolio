import numpy as np
import matplotlib.pyplot as plt	
from matplotlib import cm
from matplotlib.ticker import LinearLocator

currentStep = 0

def plotAll(step, fig):

	# Загрузка данных из файлов
	e = np.loadtxt(f'./out/e_{step}.txt')
	p = np.loadtxt(f'./out/p_{step}.txt')
	d = np.loadtxt(f'./out/d_{step}.txt')
	x = np.loadtxt(f'./out/k_x_{step}.txt')
	y = np.loadtxt(f'./out/k_y_{step}.txt')
	v_x = np.loadtxt(f'./out/v_x_{step}.txt')
	v_y = np.loadtxt(f'./out/v_y_{step}.txt')
	x_y_num = np.loadtxt(f'./out/x_y_num.txt')

	if x_y_num[2] == 0:
		n_all = int(x_y_num[0])
	else:
		n_all = int(x_y_num[1])

	p = p.reshape((n_all,1))
	e = e.reshape((n_all,1))
	d = d.reshape((n_all,1))
	v_x = v_x.reshape((n_all,1))
	v_y = v_y.reshape((n_all,1))

	fig.axes[0].cla()
	fig.axes[1].cla()
	fig.axes[2].cla()
	fig.axes[3].cla()
	fig.axes[4].cla()
	fig.axes[5].cla()

	axs[0, 0].scatter(x, y, e, c = e, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	axs[0, 0].set_title('Энергия')

	axs[1, 1].scatter(x, y, d, c = d, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	axs[1, 1].set_title('Давление')

	axs[0, 2].scatter(x, y, v_x, c = v_x, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	axs[0, 2].set_title('Скорость X')

	axs[1, 2].scatter(x, y, v_y, c = v_y, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	axs[1, 2].set_title('Скорость Y')

	axs[1, 0].scatter(x, y, p, c = p, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	axs[1, 0].set_title('Плотность')

	fig.suptitle('Значения частиц')

def handle_key_press_event(evt):
	global currentStep
	if (evt.key == "a" or evt.key == "ф") and currentStep > 0:
		currentStep -= 1
		print("currentStep = ", currentStep)
		plotAll(currentStep, evt.canvas.figure)
	if (evt.key == "s" or evt.key == "ы"):
		currentStep += 1
		print("currentStep = ", currentStep)
		plotAll(currentStep, evt.canvas.figure)
	if (evt.key == "f" or evt.key == "а"):
		currentStep += 50
		print("currentStep = ", currentStep)
		plotAll(currentStep, evt.canvas.figure)
	if (evt.key == "d" or evt.key == "в"):
		currentStep = 0
		print("currentStep = ", currentStep)
		plotAll(currentStep, evt.canvas.figure)

plt.rcParams['keymap.save'].remove('s')
fig, axs = plt.subplots(2, 3, figsize=(10, 10), subplot_kw = dict(projection='3d'))
fig.canvas.mpl_connect('key_press_event', handle_key_press_event)
plotAll(currentStep, fig)
plt.show()