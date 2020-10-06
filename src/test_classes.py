class RayCalc:

	def __init__(self):
		pass

	def doit(self):
		self.show_info()
		self.calc_params()
		self.form_sources()
		self.multiproc_calc()
		self.show_end_info()

	def show_info(self):
		print('show_info()')

	def calc_params(self):
		print('calc_params()')

	def form_sources(self):
		print('form_sources()')

	def multiproc_calc(self):
		print('multiproc_calc()')

	def show_end_info(self):
		print('show_end_info()')


class RayCalcIdealFull(RayCalc):

	def form_sources(self):
		print('form_sources for Ideal Full')


class RayCalcIdealRibs(RayCalc):

	def form_sources(self):
		print('form_sources for Ideal Ribs off')


def main():
	# RayCalc().doit()
	# RayCalcIdealFull().doit()
	RayCalcIdealRibs().doit()

if __name__ == '__main__':
	main()
