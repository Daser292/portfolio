import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np

currentStep = 0

def plotAll(step, fig):

    anl_x = np.loadtxt(f'./out/k_x_{step}.txt')
    anl_y = np.loadtxt(f'./out/p_{step}.txt')

    fig.axes[0].cla()
    
    fig.axes[0].scatter(anl_x, anl_y, color='blue')   

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
fig, axis = plt.subplots(1, 1)
fig.canvas.mpl_connect('key_press_event', handle_key_press_event)
plotAll(currentStep, fig)
plt.show()