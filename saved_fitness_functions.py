def evaluate(self, individual, plot=False):
	if self.random_state is None:
		self.random_state = random.getstate()
	random.setstate(self.random_state)
	step = int(self.config["evolution_live_step"])
	individual.init_life()

	if len(individual.genes) == 0:
		return 0
	# high = 2 / len(individual.genes)
	# mid = 1 / len(individual.genes)
	# low = 1 / 2 * len(individual.genes)
	# goal = [low, low, mid, high, high, mid, low, mid, high]
	concentrations = []
	plot_data = []
	for i in range(10):
		plot_data += individual.live(step)
		conc_list = []
		for gene in individual.genes:
			conc_list.append(gene.protein_concentration)
		concentrations.append(conc_list)
	reward = 0
	switcher = True
	for phase in concentrations:
		try:
			if switcher:
				if phase[0] - phase[1] > 1 / len(individual.genes):
					# print(phase[0], phase[1])
					reward += 1
					switcher = False
			else:
				if phase[1] - phase[0] > 1 / len(individual.genes):
					reward += 1
					switcher = True
		except:
			return 0
	fitness = reward

	# error = 0
	# for index, item in enumerate(goal):
	# 	error += abs(item - concentrations[index])

	# fitness = 1 / (1 + error)

	# if concentrations == [concentrations[0]] * len(concentrations):
	# 	return 0
	#
	# statistics = stats.pearsonr(goal, concentrations)
	# fitness = abs(statistics[0])  # r-value is index 0, and p-value is index 1
	#
	if plot:
		self.plot_individual(plot_data)
	return fitness



	def evaluate(self, individual, plot=False):
		if self.random_state is None:
			self.random_state = random.getstate()
		random.setstate(self.random_state)
		step = int(self.config["evolution_live_step"])
		individual.init_life()

		if len(individual.genes) == 0:
			return 0
		high = 2 / len(individual.genes)
		mid = 1 / len(individual.genes)
		low = 1 / 2 * len(individual.genes)
		goal = [low, low, mid, high, high, mid, low, mid, high]
		concentrations = []
		plot_data = []
		for i in range(len(goal)):
			plot_data += individual.live(step)
			concentrations.append(individual.genes[0].protein_concentration)

		error = 0
		for index, item in enumerate(goal):
			error += abs(item - concentrations[index])

		fitness = 1 / (1 + error)

		if plot:
			self.plot_individual(plot_data)

		return fitness


	def evaluate(self, individual, plot=False):
		if self.random_state is None:
			self.random_state = random.getstate()
		random.setstate(self.random_state)
		step = int(self.config["evolution_live_step"])
		individual.init_life()
		individual.live(10)
		if len(individual.genes) == 0:
			return 0

		concentrations = []
		plot_data = []
		for i in range(200):
			plot_data += individual.live(step)
			concentrations.append(individual.genes[0].protein_concentration)

		error = 0
		for i in range(200):
			error += abs(math.sin(i) - concentrations[i])

		fitness = 1 / (1 + error)

		if plot:
			self.plot_individual(plot_data)

		return fitness


	def evaluate(self, individual, plot=False):
		if self.random_state is None:
			self.random_state = random.getstate()
		random.setstate(self.random_state)
		step = int(self.config["evolution_live_step"])
		individual.init_life()
		conc_array, rate_array = individual.live(step)
		if len(individual.genes) == 0:
			return 0

		error = abs(0.085 - individual.genes[0].protein_concentration)
		fitness = 1 / (1 + error)
		if plot:
			random.setstate(self.random_state)
			individual.init_life()
			conc_array, rate_array = individual.live(500)
			self.plot_individual(conc_array)
		return fitness