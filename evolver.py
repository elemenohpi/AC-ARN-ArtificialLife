import copy
import math
import random
import time

from cell import Cell
import matplotlib.pyplot as plt

from scipy import stats


class Evolver:
	def __init__(self, config):
		self.plot_name = "-"
		self.config = config
		self.gens = int(config["gens"])
		self.tournament_size = int(config["tournament_size"])
		self.random_state = None
		self.plot_counter = 0

	def evolve(self, organisms):
		output_array = []
		for index, individual in enumerate(organisms):
			organisms[index].fake_id = index
		for gen in range(self.gens):
			self.update_fitness(organisms, gen)
			self.sort_population(organisms)
			elite = copy.deepcopy(organisms[0])
			conc_array, rate_array = elite.live(500)
			self.plot_individual(conc_array)
			total_fitness = 0
			for organism in organisms:
				total_fitness += organism.fitness
			avg_fitness = total_fitness / len(organisms)
			print("gen {} best {} avg {}".format(gen, round(elite.fitness, 5), round(avg_fitness, 5)))
			output_array.append(elite.fitness)
			organisms = self.tournament(elite, organisms)
		return output_array

	def mutate(self, offspring):
		rate = float(self.config["mutation_rate"])
		for index, _ in enumerate(offspring.dna):
			if random.random() < rate:
				offspring.dna[index] = random.choice(["A", "G", "C", "T"])

	def crossover(self, p1, p2):
		if random.random() > float(self.config["crossover_rate"]):
			return p1, p2
		destructible_p1 = copy.deepcopy(p1)
		destructible_p2 = copy.deepcopy(p2)
		lower_dna_length = min(len(p1.dna), len(p2.dna))
		crossover_point = random.randint(0, lower_dna_length-1)
		o1, o2 = Cell(self.config), Cell(self.config)
		o1.dna = destructible_p1.dna[0:crossover_point] + destructible_p2.dna[crossover_point:]
		o2.dna = destructible_p2.dna[0:crossover_point] + destructible_p1.dna[crossover_point:]
		return o1, o2

	def tournament(self, elite, organisms):
		next_generation = [copy.deepcopy(elite)]
		organisms_copy = copy.deepcopy(organisms)
		while len(next_generation) < len(organisms):
			tournament_list = []
			for i in range(self.tournament_size):
				tournament_list.append(random.choice(organisms_copy))
			self.sort_population(tournament_list)
			parent_a, parent_b = copy.deepcopy(tournament_list[0]), copy.deepcopy(tournament_list[1])
			offspring_a, offspring_b = self.crossover(parent_a, parent_b)
			self.mutate(offspring_a)
			self.mutate(offspring_b)
			next_generation.append(copy.deepcopy(offspring_a))
			next_generation.append(copy.deepcopy(offspring_b))
		return copy.deepcopy(next_generation)

	def update_fitness(self, organisms, gen):
		for index, organism in enumerate(organisms):
			plot = False
			if index == 0 and gen == self.gens:
				plot = True

			organism.fitness = self.evaluate(organism, plot)

	def sort_population(self, organisms):
		organisms.sort(key=lambda x: x.fitness, reverse=False)

	def evaluate(self, individual, plot=False):
		if self.random_state is None:
			self.random_state = random.getstate()
		random.setstate(self.random_state)
		step = int(self.config["evolution_live_step"])
		individual.init_life()
		reward = 0
		side_flag = True
		if len(individual.genes) < 1:
			return 10
		# for i in range(10):
		# 	conc_array, rate_array = individual.live(step)
		# 	if side_flag and individual.genes[0].protein_concentration > individual.genes[1].protein_concentration:
		# 		reward += 1
		# 		side_flag = not side_flag
		# 	if not side_flag and individual.genes[0].protein_concentration < individual.genes[1].protein_concentration:
		# 		reward += 1
		# 		side_flag = not side_flag

		conc_array, rate_array = individual.live(100)
		deviation = abs(individual.genes[0].protein_concentration - 0.085)

		fitness = deviation
		plot = False
		if plot:
			random.setstate(self.random_state)
			individual.init_life()
			conc_array, rate_array = individual.live(500)
			self.plot_individual(conc_array)

		return fitness

	def plot_protein(self, data, pid):
		output = []
		for item in data:
			output += [item[pid]]
		x = []
		i = 1
		while len(x) < len(data):
			x.append(i)
			i = i + 1
		plt.plot(x, output)
		plt.savefig("Output/plots/" + str(self.plot_counter) + ".png")
		plt.clf()
		self.plot_counter += 1

	def plot_individual(self, data):
		output = []
		for i in range(len(data[0])):
			output.append([])
		for item in data:
			for index, element in enumerate(item):
				output[index].append(element)

		x = []
		i = 1
		while len(x) < len(data):
			x.append(i)
			i = i + 1

		for index, item in enumerate(output):
			plt.plot(x, item, label="Protein " + str(index + 1))
		plt.legend(loc="lower right")
		plt.xlabel("Time (cycle)")
		plt.ylabel("Concentration")
		# print("plotted as ", "Output/plots/" + self.plot_name + ".png")
		plt.savefig("Output/plots/" + self.plot_name + ".png")
		plt.clf()
