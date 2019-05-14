import pylab
import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as lines

DPI = 200

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
def drawMethodConvergence(name, f, index, x0, y0, xExp, yExp):
	global text
	inputSteps(name + '.txt')
	drawFines()
	x, y, z = buildIsoLines(f)
	cs = pylab.contour(x, y, z, 25)
	plt.plot(xList, yList, linewidth=1)
	for i in range(len(xList)):
		plt.scatter(xList[i], yList[i], s=2, color='black')
	plt.title(name, fontsize=19)
	plt.xlabel('X', fontsize=10)
	plt.ylabel('Y', fontsize=10)
	plt.tick_params(axis='both', labelsize=8)
	plt.scatter(x0, y0, s=20)
	plt.scatter(xExp, yExp, s=20)
	plt.text(2, 2, text[index], size=12,
         ha="center", va="center",
         bbox=dict(boxstyle="square",
                   ec=(.5, .5, .5),
                   fc=(.8, .8, 0.8),
                   )
         )
	plt.savefig(folder + name + '.png', dpi=DPI)
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
		plt.scatter(xList[i], yList[i], s=5)
	plt.savefig(folder + name + '.png', dpi=DPI)
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

	inputMethodResults('Rosenbrock_Q.txt')
	drawMethodConvergence('Rosenbrock_Q1', f, 0, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q2', f, 1, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q3', f, 2, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q4', f, 3, x0Barrier, y0Barrier, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q5', f, 4, x0Barrier, y0Barrier, xExp, yExp)
	
	
	inputMethodResults('Rosenbrock_Q_1.txt')
	drawMethodConvergence('Rosenbrock_Q1_1', f, 0, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q2_1', f, 1, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q3_1', f, 2, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q4_1', f, 3, x0Barrier, y0Barrier, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q5_1', f, 4, x0Barrier, y0Barrier, xExp, yExp)
	

	inputMethodResults('Rosenbrock_Q_2.txt')
	drawMethodConvergence('Rosenbrock_Q1_2', f, 0, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q2_2', f, 1, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q3_2', f, 2, x0Fine, y0Fine, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q4_2', f, 3, x0Barrier, y0Barrier, xExp, yExp)
	drawMethodConvergence('Rosenbrock_Q5_2', f, 4, x0Barrier, y0Barrier, xExp, yExp)
	


	drawPointsAtResearch('tableEFines', f)
	drawPointsAtResearch('tableEBarriers', f)
	
	drawPointsAtResearch('tableRMultFines', f)
	drawPointsAtResearch('tableRMultBarriers', f)

	drawPointsAtResearch('tableRFirstFines', f)
	drawPointsAtResearch('tableRFirstBarriers', f)