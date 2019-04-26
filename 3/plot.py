import pylab
import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lines

xList = []
yList = []
data = []

folder = 'pics/'
steps = 'steps/'

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
def f(x, y):
	return 4*(y-x)**2 + 3*(x-1)**2


# Рассчёт значений целевой функции на сетке
def buildIsoLines(f):
	x = numpy.linspace (-3, 3, 100)
	y = numpy.linspace (-4, 3, 100)
	xgrid, ygrid = numpy.meshgrid(x, y)
	zgrid = f(xgrid, ygrid)
	return xgrid, ygrid, zgrid


# Закраска областей, где накладывается штраф или барьер уходит в бесконечность
def drawFines():
	x1 = numpy.linspace(-3, 3, 100)
	y1 = -1 - x1
	yBorder = 3
	plt.fill_between(x1,y1,yBorder, color='grey')
	x2 = numpy.linspace(-3, 2, 100)
	y2 = x2 + 1
	plt.plot(x2, y2)


# Отрисовка сходимости метода
def drawMethodConvergence(name, f, x0, y0, xExp, yExp):
	inputSteps(name + '.txt')
	drawFines()
	x, y, z = buildIsoLines(f)
	cs = pylab.contour(x, y, z, 25)
	plt.plot(xList, yList, linewidth=1)
	plt.title(name, fontsize=19)
	plt.xlabel('X', fontsize=10)
	plt.ylabel('Y', fontsize=10)
	plt.tick_params(axis='both', labelsize=8)
	plt.scatter(x0, y0, s=20)
	plt.scatter(xExp, yExp, s=20)
	plt.savefig(folder + name + '.png')
	plt.clf()


# Отрисовка точек в исследовании
def drawPointsAtResearch(name, f):
	inputSteps(name + '.txt')
	drawFines()
	x, y, z = buildIsoLines(f)
	cs = pylab.contour(x, y, z, 25)
	plt.title(name, fontsize=19)
	plt.xlabel('X', fontsize=10)
	plt.ylabel('Y', fontsize=10)
	plt.tick_params(axis='both', labelsize=8)
	for i in range(len(xList)):
		plt.scatter(xList[i], yList[i], s=20)
	plt.savefig(folder + name + '.png')
	plt.clf()



if __name__ == '__main__':
	# парсим начальное приближение
	params = sys.argv[1:5]
	x0Fine = float(params[0])
	y0Fine = float(params[1])
	x0Barrier = float(params[2])
	y0Barrier = float(params[3])
	xExp = 1
	yExp = 1

	drawMethodConvergence('Rosenbrock_Q1', f, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q2', f, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q3', f, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q4', f, x0Barrier, y0Barrier, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q5', f, x0Barrier, y0Barrier, xExp, yExp)
	
	drawPointsAtResearch('tableEFines', f)
	drawPointsAtResearch('tableRFirstFines', f)
	drawPointsAtResearch('tableRMultFines', f)

	drawPointsAtResearch('tableEBarriers', f)
	drawPointsAtResearch('tableRFirstBarriers', f)
	drawPointsAtResearch('tableRMultBarriers', f)
