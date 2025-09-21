import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np

currentStep = 0
firstplot = True

def plotAll(step, fig):
    global firstplot

    f1 = open('./out/x_y_num.txt', 'r')
    x_num = int(f1.readline())
    y_num = int(f1.readline())
    x_num_f = int(f1.readline())
    y_num_f = int(f1.readline())
    sct = int(f1.readline())
    lines = int(f1.readline())
    x_num_pole = int(f1.readline())
    y_num_pole = int(f1.readline())
    pole = int(f1.readline())

    v_y = np.loadtxt(f'./out/v_y_{step}.txt')
    v_x = np.loadtxt(f'./out/v_x_{step}.txt')
    x = np.loadtxt(f'./out/k_x_{step}.txt')
    y = np.loadtxt(f'./out/k_y_{step}.txt')
    e = np.loadtxt(f'./out/e_{step}.txt')
    d = np.loadtxt(f'./out/d_{step}.txt')
    p = np.loadtxt(f'./out/p_{step}.txt')
    t = np.loadtxt(f'./out/t_{step}.txt')

    if x_num == 1 and y_num == 1:
        print(v_x[-1], v_y[-1])

    zoom = 5

    line_1_1 = [zoom, -zoom]
    line_2_1 = [zoom, zoom]
    line_1_2 = [zoom, -zoom]
    line_2_2 = [-zoom, -zoom]
    line_1_3 = [zoom, zoom]
    line_2_3 = [zoom, -zoom]
    line_1_4 = [-zoom, -zoom]
    line_2_4 = [zoom, -zoom]

    if sct == 0:
        if pole == 1:
            e = e.reshape((y_num_pole, x_num_pole))
            t = t.reshape((y_num_pole, x_num_pole))
            p = p.reshape((y_num_pole, x_num_pole))
            d = d.reshape((y_num_pole, x_num_pole))
            v_x = v_x.reshape((y_num_pole, x_num_pole))
            v_y = v_y.reshape((y_num_pole, x_num_pole))
        else:
            e = e.reshape((y_num + y_num_f * 2, x_num + x_num_f * 2))
            t = t.reshape((y_num + y_num_f * 2, x_num + x_num_f * 2))
            p = p.reshape((y_num + y_num_f * 2, x_num + x_num_f * 2))
            d = d.reshape((y_num + y_num_f * 2, x_num + x_num_f * 2))
            v_x = v_x.reshape((y_num + y_num_f * 2, x_num + x_num_f * 2))
            v_y = v_y.reshape((y_num + y_num_f * 2, x_num + x_num_f * 2))

    firstplot = False
    fig.axes[0].cla()
    fig.axes[1].cla()
    fig.axes[2].cla()
    fig.axes[3].cla()
    fig.axes[4].cla()
    fig.axes[5].cla()
    
    if sct==1 and lines ==1:
        fig.axes[0].plot(line_1_1, line_2_1, color='red')
        fig.axes[0].plot(line_1_3, line_2_3, color='red')
        fig.axes[0].plot(line_1_2, line_2_2, color='red')
        fig.axes[0].plot(line_1_4, line_2_4, color='red')

        fig.axes[1].plot(line_1_1, line_2_1, color='red')
        fig.axes[1].plot(line_1_2, line_2_2, color='red')
        fig.axes[1].plot(line_1_3, line_2_3, color='red')
        fig.axes[1].plot(line_1_4, line_2_4, color='red')

        fig.axes[2].plot(line_1_1, line_2_1, color='red')
        fig.axes[2].plot(line_1_2, line_2_2, color='red')
        fig.axes[2].plot(line_1_3, line_2_3, color='red')
        fig.axes[2].plot(line_1_4, line_2_4, color='red')

        fig.axes[3].plot(line_1_1, line_2_1, color='red')
        fig.axes[3].plot(line_1_2, line_2_2, color='red')
        fig.axes[3].plot(line_1_3, line_2_3, color='red')
        fig.axes[3].plot(line_1_4, line_2_4, color='red')

        fig.axes[4].plot(line_1_1, line_2_1, color='red')
        fig.axes[4].plot(line_1_2, line_2_2, color='red')
        fig.axes[4].plot(line_1_3, line_2_3, color='red')
        fig.axes[4].plot(line_1_4, line_2_4, color='red')

        fig.axes[5].plot(line_1_1, line_2_1, color='red')
        fig.axes[5].plot(line_1_2, line_2_2, color='red')
        fig.axes[5].plot(line_1_3, line_2_3, color='red')
        fig.axes[5].plot(line_1_4, line_2_4, color='red')

    if sct == 0:
        fig.axes[0].imshow(v_y, cmap='spring')
        fig.axes[0].set_title("v_y")

        fig.axes[1].imshow(v_x, cmap='spring')
        fig.axes[1].set_title("v_x")

        fig.axes[2].imshow(e, cmap='spring')
        fig.axes[2].set_title("Energy")

        fig.axes[3].imshow(d, cmap='spring')
        fig.axes[3].set_title("Density")

        fig.axes[4].imshow(p, cmap='spring')
        fig.axes[4].set_title("Pressure")

        fig.axes[5].imshow(t, cmap='spring')
        fig.axes[5].set_title("Temperature")

    else:
        fig.axes[0].scatter(x, y, c=v_y, cmap='spring')
        fig.axes[0].set_title("v_y")

        fig.axes[1].scatter(x, y, c=v_x, cmap='spring')
        fig.axes[1].set_title("v_x")

        fig.axes[2].scatter(x, y, c=e)
        fig.axes[2].set_title("Energy")

        fig.axes[3].scatter(x, y, c=d, cmap='spring')
        fig.axes[3].set_title("Density")

        fig.axes[4].scatter(x, y, c=p, cmap='spring')
        fig.axes[4].set_title("Pressure")

        fig.axes[5].scatter(x, y, c=t, cmap='spring')
        fig.axes[5].set_title("Temperature")

    fig.canvas.draw()

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
fig, axis = plt.subplots(2, 3)
fig.canvas.mpl_connect('key_press_event', handle_key_press_event)
plotAll(currentStep, fig)
plt.show()

