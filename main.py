from nature import Nature
import random
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import argparse
from datetime import datetime
import eletility
from tf import TF
from evolver import Evolver
from scipy import stats
from cell import Cell
from statsmodels.tsa.stattools import grangercausalitytests

F = eletility.Files()


def main():
    config_parser = eletility.ConfigParser()
    config = config_parser.read("config.ini")

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-config', help='Runs a set of experiments based on config file parameters', action='store_true')
    arg_parser.add_argument('-plotC', help='Plots the latest conc file', action='store_true')
    arg_parser.add_argument('-plotR', help='Plots the latest rate file', action='store_true')
    arg_parser.add_argument('-xp', help='Runs the given number of experiments')

    args = arg_parser.parse_args()

    if args.config:
        seed_list = eval(config["seed_list"])
        run_xp_seed(seed_list, config)
        return

    if args.plotC:
        df = pd.read_csv(config["concentration_file"])
        fig = px.line(df)
        fig.show()
        return
    
    if args.plotR:
        df = pd.read_csv(config["output_file"])
        fig = px.line(df)
        fig.show()
        return

    if args.xp:
        count = int(args.xp)
        run_experiments(count, config)
        return

    random.seed(int(config["seed"]))

    Pop = Nature(config)
    Pop.create_organisms()

    # Evolution
    all_data = []
    for i in range(10):
        Pop = Nature(config)
        Pop.create_organisms()

        evolver = Evolver(config)
        evolver.plot_name = "exp_" + str(i+1)
        data = evolver.evolve(Pop.organisms)
        all_data.append(data)
    plot_evolution(all_data)

    exit()

    # # Dynamic Generation
    # for index, organism in enumerate(Pop.organisms):
    #     organism.init_life()
    #     conc_array, rate_array = organism.live(int(config["cycles"]))
    #     plot_individual(conc_array)
    #     # plot_save_individual(conc_array, "plot_" + str(index))
    #     # plot_save_rate(rate_array, "rate_" + str(index))
    #
    # exit()
    # Initial State Comparison
    conc_array = []
    comparison_data = []
    protein_id = 0
    for index, organism in enumerate(Pop.organisms):
        start_random_state = random.getstate()
        move_random_state = None
        labels = []
        starting_concentration_conditions = ["(1/len(self.genes))", "0.1", "random.random()"]
        labels = ["conc = 1/N", "0", "random"]
        gene_positions = [[(11, 10), (0, 2)], [(10, 10), (2, 1)], [(10, 11), (1, 3)], [(10, 11), (2, 2)]]

        for i in range(2):
            # random.setstate(start_random_state)

            # organism.delta = 1 - i * 0.1  # change this
            # organism.TF_count = 25 - i * 5  # change this
            # organism.cell_size = 10 + i * 5  # change this
            # mock_mutate_init(organism, i, start_random_state)  # change this
            # organism.starting_concentrations = starting_concentration_conditions[i]
            # random.setstate(start_random_state)
            organism.init_life()
            for g_i, gene in enumerate(organism.genes):
                # print(gene.ipos, gene.epos)
                organism.genes[g_i].ipos, organism.genes[g_i].epos = gene_positions[g_i][0], gene_positions[g_i][1]


            if i == 1:
                print(organism.genes[0].inhibitor)
                for gene in organism.genes:
                    print(gene.protein)
                print(organism.genes[0].epos)
                organism.genes[0].epos = (1, 2)

            #     organism.transcription_factors[0].pos = (3, 5)
            #     organism.transcription_factors[1].pos = (3, 5)
            #     organism.transcription_factors[2].pos = (3, 5)
            #     organism.transcription_factors[3].pos = (3, 5)
            #     for i in range(len(organism.transcription_factors)):
            #         organism.transcription_factors[i].pos = (8, 5)
            # if i <= 2:
              # change this
            # elif i == 2:
            #     continue
            # if move_random_state is None:
            #     move_random_state = random.getstate()
            # else:
            #     random.setstate(move_random_state)
            # random.setstate(start_random_state)
            plot_array = []
            # if i == 1:
            #     for tf_i, tf in enumerate(organism.transcription_factors):
            #         if tf_i < 5:
            #             # print(organism.transcription_factors[tf_i].gene_id)
            #             organism.transcription_factors[tf_i].gene_id = 0
            #             organism.transcription_factors[tf_i].protein = organism.genes[0].protein
            #             # print(organism.transcription_factors[tf_i].gene_id)
            #
            #         elif tf_i < 80:
            #             organism.transcription_factors[tf_i].gene_id = 1
            #             organism.transcription_factors[tf_i].protein = organism.genes[1].protein
            #         elif tf_i < 85:
            #             organism.transcription_factors[tf_i].gene_id = 2
            #             organism.transcription_factors[tf_i].protein = organism.genes[2].protein
            #         else:
            #             organism.transcription_factors[tf_i].gene_id = 5
            #             organism.transcription_factors[tf_i].protein = organism.genes[3].protein
            # g = [0, 0, 0, 0, 0, 0, 0, 0]
            # for tf in organism.transcription_factors:
            #     g[tf.gene_id] += 1
            # print(g)

            # #################################################
            # Changing TF counts for each rather than equal starting values
            # #################################################
            # distribution = [0.17, 0.03, 0.2, 0.6]
            # if i == 1:
            #     for id, gene in enumerate(organism.genes):
            #         organism.genes[id].protein_concentration = 0
            # if i == 2:
            #     for id, gene in enumerate(organism.genes):
            #         organism.genes[id].protein_concentration = distribution[id]/100
            #     organism.transcription_factors = []
            #     index = 0
            #     counter = 0
            #     while True:
            #         # Random Pos
            #         # pos = (random.randint(0, self.cellsize), random.randint(0, self.cellsize))
            #         if counter <= distribution[index]:
            #             tf = TF()
            #             tf.marker = organism.genes[index].protein
            #             tf.pos = (0, 0)
            #             tf.gene_id = organism.genes[index].id
            #             organism.transcription_factors.append(tf)
            #             counter += 1
            #             # print("gene", index, "protein", counter)
            #         elif counter > distribution[index]:
            #             index += 1
            #             counter = 0
            #         if index == 4:
            #             break
            # #################################################
            # END TF Count Change
            # #################################################

            conc_array, rate_array = organism.live(int(config["cycles"]))
            # conc_array, rate_array = organism.live(500)
            # plot_save_individual(conc_array, "plot_conc_" + str(index))
            plot_individual(conc_array)
            # plot_array += conc_array

            random_values = [[(0, 4), (1, 5)], [(5, 5), (7, 2)], [(1, 3), (2, 1)], [(10, 9), (3, 5)]]
            random_values = [[(3, 1), (2, 2)], [(4, 9), (2, 6)], [(2, 7), (1, 6)], [(2, 1), (7, 9)]]
            # print(organism.genes[0].epos)
            # organism.genes[0].epos = (0, 0)
            # for gene_id in range(len(organism.genes)):
            #     organism.genes[gene_id].protein_concentration = random.random()
            #     organism.genes[gene_id].epos = random_values[gene_id][0]
            #     organism.genes[gene_id].ipos = random_values[gene_id][1]
            #     pass
            # for tf_i, tf in enumerate(organism.transcription_factors):
            #     organism.transcription_factors[tf_i].binding_strength = 0
            #     organism.transcription_factors[tf_i].pos = (0, 0)
            #     organism.transcription_factors[tf_i].max_enh = -1
            #     organism.transcription_factors[tf_i].max_inh = -1
            #     organism.transcription_factors[tf_i].bound_gene = -1
            #     if tf_i < 5:
            #         organism.transcription_factors[tf_i].gene_id = 0
            #         organism.transcription_factors[tf_i].protein = organism.genes[0].protein
            #     elif tf_i < 60:
            #         organism.transcription_factors[tf_i].gene_id = 1
            #         organism.transcription_factors[tf_i].protein = organism.genes[1].protein
            #     elif tf_i < 95:
            #         organism.transcription_factors[tf_i].gene_id = 2
            #         organism.transcription_factors[tf_i].protein = organism.genes[2].protein
            #     else:
            #         organism.transcription_factors[tf_i].gene_id = 3
            #         organism.transcription_factors[tf_i].protein = organism.genes[3].protein

            plot_array += conc_array

            conc_array = transform_to_plot_data(plot_array)
            # plot_individual(plot_array)
            # exit()
            comparison_data.append(conc_array[protein_id])
            # labels.append("conc = " + starting_concentration_conditions[i])  # change this
            # labels.append("m = " + str(i))  # change this
        # exit()
        plot_compare(comparison_data, file_name="conc_" + str(index), labels=labels, show=True, loc="upper right")  # change this

        # random.setstate(random_state)
        # organism.beta = 1
        # organism.delta = 1
        # organism.init_life()
        # conc_array[index + 1], rate_array = organism.live(int(config["cycles"]))
        # plot_individual(conc_array[index + 1])
        # correl = get_correlation(conc_array[index], conc_array[index + 1])
        # print("Correlation between two settings: ", correl)

    # conc_array = [[], []]
    # for index, organism in enumerate(Pop.organisms):
    #     random_state = random.getstate()
    #     organism.init_life()
    #     conc_array[index], rate_array = organism.live(int(config["cycles"]))
    #     plot_individual(conc_array[index])
    #
    #     random.setstate(random_state)
    #     organism.beta = 1
    #     organism.delta = 1
    #     organism.init_life()
    #     conc_array[index + 1], rate_array = organism.live(int(config["cycles"]))
    #     plot_individual(conc_array[index+1])
    #     correl = get_correlation(conc_array[index], conc_array[index+1])
    #     print("Correlation between two settings: ", correl)

    # correl = get_correlation(conc_array[0], conc_array[1])
    # print("Correlation between two settings: ", correl)
    # exit()
    # df = pd.read_csv(config["concentration_file"])
    # fig = px.line(df, title='Gene Dynamics')
    # fig.show()


def mock_mutate_init(organism, mutation_count, random_state):
    if mutation_count <= 0:
        organism.init_life()
        return
    organism.identify_genes()
    gene_count = len(organism.genes)
    random_gene_id = random.randint(0, gene_count-1)
    # random_gene_id = 0
    random.seed(107)
    for i in range(mutation_count):
        region = random.choice(["enhancer", "inhibitor", "protein"])
        # region = "protein"
        # print("random mutation region: ", region)
        if region == "enhancer":
            region_size = len(organism.genes[random_gene_id].enhancer)
        elif region == "inhibitor":
            region_size = len(organism.genes[random_gene_id].inhibitor)
        else:  # protein
            region_size = len(organism.genes[random_gene_id].protein)
        random_bp_index = random.randint(0, region_size-1)
        # random_bp_index = 5
        base_pairs = ["A", "G", "C", "T"]
        if region == "inhibitor":
            old_bp = organism.genes[random_gene_id].enhancer[random_bp_index]
            base_pairs.remove(old_bp)
            new_bp = random.choice(base_pairs)
            # new_bp = "T"
            organism.genes[random_gene_id].enhancer[random_bp_index] = new_bp
        elif region == "enhancer":
            old_bp = organism.genes[random_gene_id].inhibitor[random_bp_index]
            base_pairs.remove(old_bp)
            new_bp = random.choice(base_pairs)
            # new_bp = "A"
            organism.genes[random_gene_id].inhibitor[random_bp_index] = new_bp
        elif region == "protein":  # protein
            old_bp = organism.genes[random_gene_id].protein[random_bp_index]
            base_pairs.remove(old_bp)
            new_bp = random.choice(base_pairs)
            # new_bp = "T"
            organism.genes[random_gene_id].protein[random_bp_index] = new_bp
        else:
            raise ValueError("Wrong region specified for mutation")
        print("gene", random_gene_id, "mutation: index", random_bp_index, "of", region, "changed from", old_bp, "to", new_bp)
    random.setstate(random_state)
    organism.phony_init_life()


def transform_to_plot_data(data):
    output = []
    for i in range(len(data[0])):
        output.append([])
    for item in data:
        for index, element in enumerate(item):
            output[index].append(element)
    return output


def get_correlation(s1, s2):
    o1 = []
    for i in range(len(s1[0])):
        o1.append([])
    for item in s1:
        for index, element in enumerate(item):
            o1[index].append(element)
    o2 = []
    for i in range(len(s2[0])):
        o2.append([])
    for item in s2:
        for index, element in enumerate(item):
            o2[index].append(element)

    sum_correl = 0
    sum_pvalue = 0
    for i in range(len(o1)):
        rvalue, pvalue = stats.pearsonr(o1[i], o2[i])
        sum_correl += rvalue
        sum_pvalue += pvalue
    avg_correl = sum_correl/len(o1)
    avg_pvalue = pvalue/len(o1)

    return avg_correl


def plot_compare(data, file_name: str, labels, show=False, loc="center right"):
    line_styles = [
        ('solid', 'solid'),
        ('dotted', 'dotted'),
        ('dashed', 'dashed'),
        ('dashdot', 'dashdot'),
        ('loosely dotted', (0, (1, 10))),
        ('dotted', (0, (1, 1))),
        ('densely dotted', (0, (1, 1))),
        ('loosely dashed', (0, (5, 10))),
        ('dashed', (0, (5, 5))),
        ('densely dashed', (0, (5, 1))),
        ('loosely dashdotted', (0, (3, 10, 1, 10))),
        ('dashdotted', (0, (3, 5, 1, 5))),
        ('densely dashdotted', (0, (3, 1, 1, 1))),
        ('dashdotdotted', (0, (3, 5, 1, 5, 1, 5))),
        ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
        ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]

    colors = ["#000000",
              "#252525",
              "#525252",
              "#969696",
              "#bdbdbd",
              "#d9d9d9",
              "#f0f0f0",
              "#737373"
              ]

    x = []
    i = 1
    while len(x) < len(data[0]):
        x.append(i)
        i = i + 1

    for index, item in enumerate(data):
        plt.plot(x, item, color=colors[index], label=labels[index], linestyle=line_styles[index][1])
    plt.legend(loc=loc)
    plt.xlabel("Time (cycle)")
    plt.ylabel("Concentration")
    if show:
        plt.show()
    else:
        plt.savefig("Output/plots/compare_" + file_name + ".png")
    plt.clf()


def plot_individual(data):
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
    plt.show()
    # plt.savefig("Output/plots/i" + str(self.plot_counter) + ".png")
    # plt.clf()
    # self.plot_counter += 1


def plot_save_individual(data, file_name: str):
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
        plt.plot(x, item, label="Protein " + str(index+1))
    plt.legend(loc="center right")
    plt.xlabel("Time (cycle)")
    plt.ylabel("Concentration")
    plt.savefig("Output/plots/" + file_name + ".png")
    plt.clf()


def plot_save_rate(data, file_name: str):
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
        if index != 4:
            continue
        plt.plot(x, item, color="C4", label="Protein " + str(index + 1))
    plt.legend(loc="center right")
    plt.xlabel("Time (cycle)")
    plt.ylabel("Production Rate")
    plt.savefig("Output/plots/rate_" + file_name + ".png")
    plt.clf()

 
def run_xp_seed(seedlist, config):
    today = datetime.today()
    now = today.strftime("%Y%b%d_%H-%M-%S")
    for i, seed in enumerate(seedlist):
        random.seed(seed)
        config["output_file"] = "experiments/" + \
            now + "_seed{1}_out_{0}.csv".format(i, seed)
        config["concentration_file"] = "experiments/" + \
            now + "_seed{1}_conc_{0}.csv".format(i, seed)
        config["visual_file"] = "experiments/" + \
            now + "_seed{1}_vis_{0}.csv".format(i, seed)

        N = Nature(config)
        N.create_organisms()

        for organism in N.organisms:
            organism.init_life()
            organism.live(int(config["runs"]))


def run_experiments(count, config):
    today = datetime.today()
    now = today.strftime("%Y%b%d_%H-%M-%S")
    for i in range(count):
        seed = random.randint(100, 1000)
        random.seed(seed)

        config["output_file"] = "experiments/" + now + "_seed{1}_out_{0}.csv".format(i, seed)
        config["concentration_file"] = "experiments/" + now + "_seed{1}_conc_{0}.csv".format(i, seed)
        config["visual_file"] = "experiments/" + now + "_seed{1}_vis_{0}.csv".format(i, seed)

        N = Nature(config)
        N.create_organisms()

        for organism in N.organisms:
            organism.init_life()
            organism.live(int(config["runs"]))


def plot_evolution(data):
    # Takes directories, returns a plot of average with q75 and q25. Can plot multiple directories to compare.
    plt.clf()

    plt.xlabel("Generations")
    plt.ylabel("Fitness (1/(1 + error))")

    q25s = []
    q75s = []
    medians = []

    gens = len(data[0])
    for i in range(gens):
        each_gen_array = []
        for j in range(len(data)):
            each_gen_array.append(data[j][i])

        best_fitness_values_at_gen_i_df = pd.DataFrame(each_gen_array)
        quantiles = best_fitness_values_at_gen_i_df.quantile([0.25, 0.75])
        medians.append(best_fitness_values_at_gen_i_df.median().values[0])
        q25s.append(quantiles[0].iloc[0])
        q75s.append(quantiles[0].iloc[1])

    x_axis_data = range(0, gens)
    plt.fill_between(x_axis_data, q25s, q75s, alpha=.1, linewidth=0)
    plt.plot(x_axis_data, medians, linewidth=1.8, label="median fitness", color="#000000")
    plt.legend(loc='lower right')
    plt.show()


if __name__ == "__main__":
    main()
