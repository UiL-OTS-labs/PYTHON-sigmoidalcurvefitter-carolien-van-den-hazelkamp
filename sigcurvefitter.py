# using http://people.duke.edu/~ccc14/pcfb/analysis.html
# En dan is lambda de asymptote, delta de intercept met de x-as (waar accuracy van kans afwijkt) en beta de steilheid van de curve. 

import sys, csv, os, re, math, sys, warnings, itertools
try:
	import matplotlib.pyplot as plt
	from matplotlib.widgets import Button
except:
	print "matplotlib library is not installed. Please install this first"
try:
	import numpy as np
	import numpy.random as npr
except:
	print "Numpy library is not installed. Please install this first"
try:
	from scipy.optimize import leastsq
except:
	print "scipy library is not installed. Please install this first"

class valuesWriter():

	writer = ""

	def __init__(self):
		self.filename = raw_input("Enter the name of the output directory you wish to use. Existing files in the directory will be overwritten:\n")
		try:
			os.stat(self.filename)
		except:
			os.mkdir(self.filename)
		self.f = open(self.filename + "/" + self.filename + "_output.csv", 'w')
		self.writeHeader()

	def writeHeader(self):
		headerrow = "PPId;condition;"
		for i in range(7):
			headerrow += ("X%d;" %(i+1))
			headerrow += ("Y%d-original;" %(i+1))
			headerrow += ("Y%d-fit;" %(i+1))
		headerrow += "lambda;delta;beta;GoF;Keep/Discard\n"
		self.f.write(headerrow)

	def shutdown(self):
		self.f.close()

	def writeRow(self, pp, condition, list1, list2, list3, vals, gof, useful):
		row = "%s;%s;" %(pp, condition)
		for i in range(len(list1)):
			row += "%f;" %list1[i]
			row += "%f;" %list2[i]
			row += "%f;" %list3[i]
		
		# Handle empty Y-values
		missingItems = 7 - len(list1)
		if missingItems > 0:
			for i in range(missingItems):
				row += ";;;"

		row += "%f;%f;%f;%f;%s\n" %(vals['A'], vals['B'], vals['C'],gof,useful)
		self.f.write(row)

	def getFilename(self):
		return self.filename


class data_entry:
	"""Makes it easier to extend."""
	def __init__(self):
		self.ppid = "placeholder"
		self.condition = 0
		self.coordinates = list()

def readfile(name):
	'''
	Read the input file.
	Creates a dictionary with the participants as key
	and a list of 7 items as value
	@arg name 	The filename and relative path for the input file
	@return 	Dictionary of participants and their items
	'''
	if not name.lower().endswith('.csv'): 
		print("The file '%s' does not seem to be a .csv file. Please select a file that ends with '.csv'" %name)
		exit(1)
	elif not os.path.isfile(name):
		print("The file '%s' could not be found. Try again with another file." %name)
		exit(1)
	else:
		print "Now reading file '%s'" %name
		result = dict()
		with open(name, 'rb') as data:
			reader = csv.reader(data, delimiter=';')
			for line,row in enumerate(reader):
				key = row[0]+"_"+row[1]
				if row[0] == "ppid": # Headerline; skip
					continue

				if not key in result:
					result[key] = data_entry()

				result[key].ppid = row[0]
				result[key].condition = row[1]
				try:
					x = float(row[2])
				except ValueError:
					continue
				try:
					y = row[3].replace(',', '.')
					y = float(y)
				except ValueError:
					print "Warning: missing or faulty value (\'%s\') on line: %d"%(row[2], line)
					continue
				result[key].coordinates.append( (x, y) )
		return result

def formula(x, A, B, C):
	return sigmoidCurve(x, A, B, C)

def sigmoidCurve(times, lambda_i, delta, beta):
	"""sigmoid curve"""
	try:
		len(times)
		result = []
		for time in times:
			if time <= delta:
				result.append(0)
			else:
				exponent = -beta*(time-delta)
				if exponent > math.log(sys.float_info.max):
					exponent = math.log(sys.float_info.max)
				result.append(lambda_i*(1-math.exp(exponent)))
	except TypeError:
		if times <= delta:
			result = 0
		else:
			exponent = -beta*(times-delta)
			if exponent > math.log(sys.float_info.max):
				exponent = math.log(sys.float_info.max)
			result = lambda_i*(1-math.exp(exponent))
	
	return result

def residuals(p, y, x):
	"""Deviations of data from fitted curve"""
	A,B,C = p
	with warnings.catch_warnings():
		warnings.filterwarnings('error')
		err = y-formula(x, A, B, C)
		
	return err

def peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C = p
    return formula(x, A, B, C)

def plot(data, sortKeys, index, figure, text_offset, auto_use):
	'''
	Plots the coordinates for a participant.
	@arg 	pp 		Participant ID
	@arg 	coords 	list of coordinates (x = time, y = accurary)
	'''
	pp = sortKeys[index]
	ppid = data[pp].ppid
	coords = data[pp].coordinates
	condition = data[pp].condition

	list1, list2 = [list(t) for t in zip(*coords)]

	# Convert x and y data to numpy array
	x = np.asarray(list1)
	y_meas = np.asarray(list2)

	# Initial guess for parameters
	p0 = [max(y_meas), 2, 0.001]

	# Fit equation using least squares optimization
	plsq = leastsq(residuals, p0, args=(y_meas, x), maxfev=10000)

	x_fitted_plot = np.linspace(1e-8,max(x),200)
	

	theplot = plt.plot(x_fitted_plot, peval(x_fitted_plot, plsq[0]),x,y_meas, 'o')
	
	plt.xlabel('Time (ms)')
	plt.ylabel('accurary')
	plt.title('Least-squares fit to %s data' %pp)
	textX = (max(x)*0.5) + text_offset
	textY = (max(y_meas)*0.4) 
	for i, (param, est) in enumerate(zip('ABC',plsq[0])):
	     plt.text(textX, textY-i*0.2, '%s = %.2f' % (param, est))
	plt.savefig(writer.getFilename() + "/" + pp + '.png')
	plt.draw()

	vals = dict()
	for i, (param, est) in enumerate(zip('ABC',plsq[0])):
	     vals[param] = est
	
	# Calculate new Y-values and Goodness of Fit (GoF)
	list3 = [formula(x, vals['A'],vals['B'],vals['C']) for x in list1]
	gof = sum([(y1 - y2)**2 for (y1, y2) in zip(list2,list3)])
	
	# Ask user how to mark these values
	saveValues(ppid, condition, list1, list2, list3, vals, gof, auto_use)

def selectFile():
	return raw_input("Enter the name of the source csv file\n")

def saveValues(pp, condition, list1, list2, list3, vals, gof, auto_use):
	quest = "For participant %s I found these values: \nlambda: %.2f\tdelta: %.2f\tbeta: %.2f\n" %(pp, vals['A'], vals['B'], vals['C'])
	quest += "The goodness of fit is %.2f\n" %gof
	quest += "The condition was %s\n" %condition
	print quest

	if auto_use:
		useful = "use"
	else:
		useful = raw_input('Do you want to mark these values as useful? The previous condition was automatically marked usefull. Press \'U\'. Otherwise, press \'n\'\nTo exit type "exit"\n')

	u = ""
	if useful in ["u", "U", "USE", "use", "Use","y","Y","yes","YES"]:
		print "Items will be used! On to the next graph.\n\n"
		u = "keep"
	elif useful in ["n", "N", "No", "no", "NO"]:
		u = "discard"
		print "These items will not be used! On to the next graph.\n\n"
	elif useful == "exit":
		#user terminated program
		writer.shutdown()
		exit(0)
	else:
		print "Your input was not recognised."
		return saveValues(pp, condition, list1, list2, list3, vals, gof, auto_use)

	writer.writeRow(pp, condition, list1, list2, list3, vals, gof, u)

if __name__ == "__main__":
	if len(sys.argv) > 1:
		selectedFile = sys.argv[1]
	else:
		selectedFile = selectFile()

	data = readfile(selectedFile)
	sortKeys = sorted(data.keys())

	writer = valuesWriter()
	
	fig1 = plt.figure()
	fig1.show()
	plt.ion()

	for i in range(0,len(sortKeys),2):
		fig1.clear()
		plot(data, sortKeys, i, fig1, 0, True)
		plot(data, sortKeys, i+1, fig1, 350, False)

	print "Done! All the entries in the file have been fitted."

	writer.shutdown()