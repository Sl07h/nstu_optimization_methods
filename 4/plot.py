import pylab
import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

DPI = 400

xList = []
yList = []
data = []
text = []

folder = 'pics/'
steps = 'steps/'

# Считывание итогов работы метода
def inputMethodResults(filename):
	global text
	global data
	text = []
	data = []
	with open(steps + filename, 'r') as f:
		for line in f: # read rest of lines
			data.append([ int(x) for x in line.split()])
	for i in range(len(data)):
		text.append('iter:'+str(data[i][0])+'\nfCalc: '+str(data[i][1]))


# Считывание хода метода
def inputSteps(filename):
	global xList
	global yList
	global data
	xList = []
	yList = []
	data = []
	with open(steps + filename, 'r') as f:
		for line in f: # read rest of lines
			data.append([ float(x) for x in line.split()])
	for i in range(len(data)):
		xList.append(data[i][0])
		yList.append(data[i][1])

	


# Целевая функция
# Вариант 5
def f(x, y):
	C = [1, 2, 10, 5, 7, 9]
	a = [0, 0, 3, -7, 6, 6]
	b = [-1, -4, -2, -6, -10, 1]
	result = 0
	for i in range(6):
		result += C[i] / (1 + (x-a[i])**2 + (y-b[i])**2)
	return result


# Рассчёт значений целевой функции на сетке
def buildIsoLines(f):
	x = numpy.linspace (-12, 12, 100)
	y = numpy.linspace (-12, 12, 100)
	xgrid, ygrid = numpy.meshgrid(x, y)
	zgrid = f(xgrid, ygrid)
	return xgrid, ygrid, zgrid


# Рассчёт значений целевой функции на сетке
def buildIsoLines2(f):
	x = numpy.linspace (2.97, 3.03, 100)
	y = numpy.linspace (-2.03, -1.97, 100)
	xgrid, ygrid = numpy.meshgrid(x, y)
	zgrid = f(xgrid, ygrid)
	return xgrid, ygrid, zgrid




def draw2D():
	name = 'var5_2d'
	fig = plt.figure()
	ax = fig.add_subplot(111)
	x, y, z = buildIsoLines(f)
	cs = pylab.contourf(x, y, z, 25, cmap=cm.coolwarm)
	plt.title(name, fontsize=19)
	plt.xlabel('X', fontsize=10)
	plt.ylabel('Y', fontsize=10)
	ax.set_xticks(numpy.linspace(-10, 10, 21))
	ax.set_yticks(numpy.linspace(-10, 10, 21))
	plt.tick_params(axis='both', labelsize=8)
	plt.grid(alpha=0.2)
	plt.savefig(folder + name + '.png', dpi=DPI)
	plt.clf()


def draw2DZoom():
	name = 'var5_2d_zoom'
	fig = plt.figure()
	ax = fig.add_subplot(111)
	x, y, z = buildIsoLines2(f)
	cs = pylab.contourf(x, y, z, 25, cmap=cm.coolwarm)
	plt.title(name, fontsize=19)
	plt.xlabel('X', fontsize=10)
	plt.ylabel('Y', fontsize=10)
	ax.set_xticks(numpy.linspace (2.97, 3.03, 11))
	ax.set_yticks(numpy.linspace (-2.03, -1.97, 11))
	plt.tick_params(axis='both', labelsize=8)
	plt.grid(alpha=0.2)
	plt.savefig(folder + name + '.png', dpi=DPI)
	plt.clf()


def draw3D():
	name = 'var5_3d'
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	X = numpy.linspace (-12, 12, 30)
	Y = numpy.linspace (-12, 12, 30)
	X, Y = numpy.meshgrid(X, Y)
	Z = f(X, Y)
	surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
			linewidth=0, antialiased=True)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.title(name, fontsize=19)
	plt.xlabel('X')
	plt.ylabel('Y')
	plt.savefig(folder + name + '.png', dpi=DPI)


def simpleRandomSearch():
	name = 'simpleRandomSearch'
	inputSteps(name + '.txt')
	fig = plt.figure()
	ax = fig.add_subplot(111)
	x, y, z = buildIsoLines(f)
	cs = pylab.contourf(x, y, z, 25, cmap=cm.coolwarm)
	plt.title(name, fontsize=19)
	plt.xlabel('X', fontsize=10)
	plt.ylabel('Y', fontsize=10)
	ax.set_xticks(numpy.linspace(-10, 10, 21))
	ax.set_yticks(numpy.linspace(-10, 10, 21))
	plt.tick_params(axis='both', labelsize=8)
	plt.grid(alpha=0.2)
	plt.plot(xList, yList, linewidth=0.5, color='grey')
	for i in range(len(xList)):
		plt.scatter(xList[i], yList[i], s=2, color='black')
	plt.savefig(folder + name + '.png', dpi=DPI)
	plt.clf()


if __name__ == '__main__':
	#simpleRandomSearch()
	draw2D()
	draw2DZoom()
	draw3D()