# Profiling manipulations


def main():
	import pstats
	p = pstats.Stats(r'c:\Downloads\Calc\output.prof')
	p.strip_dirs().sort_stats('time').print_stats()
	# p.strip_dirs().sort_stats('cumulative').print_stats()
	

if __name__ == '__main__':
	main()
