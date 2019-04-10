import pylab
import numpy
import sys
import matplotlib.pyplot as plt


x_number_values = []
y_number_values = []
data = []

folder = "pics/"
steps = "steps/"

# Считывание хода метода
def inputSteps(filename):
    global x_number_values
    global y_number_values
    global data
    x_number_values = []
    y_number_values = []
    data = []
    with open(steps + filename, 'r') as f:
        for line in f: # read rest of lines
            data.append([ float(x) for x in line.split()])
    for i in range(len(data)):
        x_number_values.append(data[i][0])
        y_number_values.append(data[i][1])
        


# Отрисовка хода метода
def draw_line(name):
    plt.plot(x_number_values, y_number_values, linewidth=1)
    plt.title(name, fontsize=19)
    plt.xlabel("X", fontsize=10)
    plt.ylabel("Y", fontsize=10)
    plt.tick_params(axis='both', labelsize=8)

    
    
def f1(x, y):
    A1 = 2
    A2 = 3
    a1 = 1
    a2 = 2
    b1 = 2
    b2 = 3
    c1 = 1
    c2 = 3
    d1 = 1
    d2 = 2
    return -(A1 / (1 + ((x - a1) / b1) ** 2 + ((y - c1) / d1) ** 2) + A2 / (1 + ((x - a2) / b2) ** 2 + ((y - c2) / d2) ** 2));
    
def f2(x, y):
    return 100 * ((y - x)**2) + ((1 - x)**2) 
    
def f3(x, y):
    return 100 * ((y - x**2)**2) + ((1 - x)**2)
    

def makeData(f):
    x = numpy.arange (-1, 2.5, 0.05)
    y = numpy.arange (-1, 3.5, 0.05)
    xgrid, ygrid = numpy.meshgrid(x, y)
    zgrid = f(xgrid, ygrid)
    return xgrid, ygrid, zgrid


def singleTable(name, f, x0, y0, xExp, yExp):
    inputSteps(name + ".txt")
    x, y, z = makeData(f)
    cs = pylab.contour(x, y, z, 25)
    draw_line(name)
    plt.scatter(x0, y0, s=20)
    plt.scatter(xExp, yExp, s=20)
    plt.savefig(folder + name + ".png")
    plt.clf()

    


if __name__ == "__main__":
    # парсим начальное приближение
    params = sys.argv[1:3]
    x0 = float(params[0])
    y0 = float(params[1])
    
    singleTable("Rosenbrock_f1", f1, x0, y0, 1.1733, 1.2109)
    singleTable("Rosenbrock_f2", f2, x0, y0, 1, 1)
    singleTable("Rosenbrock_f3", f3, x0, y0, 1, 1)
    
    singleTable("Broyden_f1", f1, x0, y0, 1.1733, 1.2109)
    singleTable("Broyden_f2", f2, x0, y0, 1, 1)
    singleTable("Broyden_f3", f3, x0, y0, 1, 1)