import copy
import random
import math

import scipy.stats

from gene import Gene
from Helpers import List, Math, Files
from tf import TF
import numpy as np


class Cell:
    # constructor
    def __init__(self, config):
        self.config = config
        self.transcription_factors = None
        self.interaction_matrix = None
        self.genes = None
        self.fitness = 0
        self.info = ""
        try:
            self.starting_concentrations = config["starting_concentrations"]
            self.dna_size = int(config["dna_size"])
            self.promoter_seq = config["promoter_sequence"]
            self.terminator_seq = config["terminator_sequence"]
            self.cell_size = int(config["cell_size"])
            self.TF_count = int(config["TF_count"])
            self.TF_bind_threshold = int(config["TF_bind_threshold"])
            self.max_inh = 0
            self.max_enh = 0
            self.beta = float(config["beta"])
            self.delta = float(config["delta"])
            # self.dissipation_rate = float(config["dissipation_rate"])
            self.visual_file = config["visual_file"]
            self.concentration_file = config["concentration_file"]
            self.rate_file = config["rate_file"]
            self.shock = config["shock"]
            if self.shock == "False":
                self.shock = False
            elif self.shock == "True":
                self.shock = True
            else:
                raise "Shock Error"
        except ValueError as e:
            raise e

        if self.dna_size < 100:
            raise Exception("DNA size is too small")
        if len(self.promoter_seq) < 1 or len(self.terminator_seq) < 1:
            raise Exception("promoter/terminator seq too small")
        for bp in self.promoter_seq:
            if not bp in ["A", "G", "C", "T"]:
                raise Exception("Unknown base-pair in the promoter sequence: " + bp)
        for bp in self.terminator_seq:
            if not bp in ["A", "G", "C", "T"]:
                raise Exception("Unknown base-pair in the terminator sequence: " + bp)

        self.dna = []

    # randomly initialize the genome
    def init_genetic_marker(self):
        base_pairs = ["A", "G", "C", "T"]
        for _ in range(self.dna_size):
            rand = random.randint(0, 3)
            self.dna.append(base_pairs[rand])
    
    def identify_genes(self):
        self.genes = []
        record_flag = False  # used to record the exon pattern into the gene exon variable
        start_index = -1
        gene = None
        for index, bp in enumerate(self.dna):
            if not record_flag and self.is_promoter(index):  # if promoter was found and we're not recording exons
                record_flag = True
                start_index = index + len(self.promoter_seq)  # current dna marker index
                gene = Gene()
            if record_flag and not self.is_terminator(index):  # exon area. record
                gene.exon.append(bp)
            if record_flag and self.is_terminator(index):  # if terminator was found. stop recording
                gene.exon = gene.exon[len(self.promoter_seq):]
                if len(gene.exon) < 6:
                    record_flag = False
                    continue
                gene.promoter = self.promoter_seq
                gene.process_exon(self.dna, start_index)
                gene.id = len(self.genes)
                gene.protein_concentration = -1
                self.genes.append(gene)
                record_flag = False
                start_index = -1
                # if len(gene.exon) < 20:
                #     print(gene.enhancer)
                #     print(gene.inhibitor)
                #     print(gene.protein)
                #     print(gene.exon)
                #     exit()

    def is_promoter(self, start):
        for index, bp in enumerate(self.promoter_seq):
            i = start + index
            if i >= len(self.dna):
                return False
            if self.dna[i] != bp:
                return False
        return True

    def is_terminator(self, start):
        for index, bp in enumerate(self.terminator_seq):
            i = (start + index) % len(self.dna)
            if self.dna[i] != bp:
                return False
        return True

    def make_network(self):
        self.interaction_matrix = List.make_matrix(len(self.genes), len(self.genes))
        for i, g1 in enumerate(self.genes):
            for j, g2 in enumerate(self.genes):
                enh, inh = self.calculate_matching_degree(g1, g2)
                # enh, inh = self.calculate_matching_degree(g1, g2)
                if enh > self.max_enh:
                    self.max_enh = enh
                if inh > self.max_inh:
                    self.max_inh = inh
                self.interaction_matrix[i][j] = (enh, inh)

    def calculate_matching_degree(self, g1: Gene, g2: Gene):
        enh, inh = 0, 0
        min_bp_count = min(len(g1.protein), len(g2.protein))
        for i in range(min_bp_count):
            if g1.protein[i] == "A" and g2.enhancer[i] == "T":
                enh += 1
            elif g1.protein[i] == "G" and g2.enhancer[i] == "C":
                enh += 1
            elif g1.protein[i] == "C" and g2.enhancer[i] == "G":
                enh += 1
            elif g1.protein[i] == "T" and g2.enhancer[i] == "A":
                enh += 1

            if g1.protein[i] == "A" and g2.inhibitor[i] == "T":
                inh += 1
            elif g1.protein[i] == "G" and g2.inhibitor[i] == "C":
                inh += 1
            elif g1.protein[i] == "C" and g2.inhibitor[i] == "G":
                inh += 1
            elif g1.protein[i] == "T" and g2.inhibitor[i] == "A":
                inh += 1
        return enh, inh
    
    def init_life(self):
        self.identify_genes()
        self.make_network()  # creates the interaction matrix of g1 protein vs. g2 enh and inh regions

        # determine entity positions
        self.transcription_factors = []
        for index, gene in enumerate(self.genes):

            start_random_state = random.getstate()
            gene.protein_concentration = eval(self.starting_concentrations)
            random.setstate(start_random_state)

            gene.gene_tf_count = self.TF_count

            i_x = random.randint(int(self.cell_size/2 - 1), int(self.cell_size/2 + 1))
            # i_x = random.randint(int(self.cell_size - 1), int(self.cell_size + 1))
            i_y = random.randint(int(self.cell_size/2 - 1), int(self.cell_size/2 + 1))
            # i_y = random.randint(int(self.cell_size - 1), int(self.cell_size + 1))
            distance = [-2, -1, 0, 1, 2]
            gene.ipos = (i_x,
                         i_y)
            gene.epos = ((random.choice(distance) + gene.ipos[0]) % self.cell_size, (random.choice(distance) + gene.ipos[1]) % self.cell_size)
            for _ in range(self.TF_count):
                # Random Pos 
                # pos = (random.randint(0, self.cellsize), random.randint(0, self.cellsize)) 
                tf = TF()
                tf.marker = gene.protein
                # tf.pos = pos 
                tf.pos = (0, 0)
                tf.gene_id = gene.id
                self.transcription_factors.append(tf)

    def phony_init_life(self):
        self.make_network()  # creates the interaction matrix of g1 protein vs. g2 enh and inh regions
        # determine entity positions
        self.transcription_factors = []
        for index, gene in enumerate(self.genes):
            gene.protein_concentration = eval(self.starting_concentrations)

            gene.gene_tf_count = self.TF_count

            i_x = random.randint(int(self.cell_size - 1), int(self.cell_size + 1))
            i_y = random.randint(int(self.cell_size - 1), int(self.cell_size + 1))
            distance = [-2, -1, 0, 1, 2]
            gene.ipos = (i_x,
                         i_y)
            gene.epos = ((random.choice(distance) + gene.ipos[0]) % self.cell_size,
                         (random.choice(distance) + gene.ipos[1]) % self.cell_size)

            for _ in range(self.TF_count):
                # Random Pos
                # pos = (random.randint(0, self.cellsize), random.randint(0, self.cellsize))
                tf = TF()
                tf.marker = gene.protein
                # tf.pos = pos
                tf.pos = (0, 0)
                tf.gene_id = gene.id
                self.transcription_factors.append(tf)

    def reset_live_genes(self):
        for gene in self.genes:
            gene.transcription_rate = 0
            gene.total_bound_conc = 0
            gene.num_bound = 0

    def live(self, cycles):
        # visual_file = self.visual_file
        # Files.truncate(visual_file)
        #
        # conc_file = self.concentration_file
        # Files.truncate(conc_file)
        #
        # output_file = self.rate_file
        # Files.truncate(output_file)
        rate_array = []
        each_rate_array = []
        conc_array = []
        each_array = []
        for gene in self.genes:
            each_rate_array.append(gene.transcription_rate)
            each_array.append(gene.protein_concentration)
        conc_array.append(each_array)
        rate_array.append(each_rate_array)

        beta = self.beta
        highest_conc_gene_id = -1
        for cycle in range(cycles):
            self.reset_live_genes()

            # Transcription Rate Update Phase
            for tf in self.transcription_factors:
                if tf.binding_strength > 0:
                    # TF is already bound to a gene
                    if tf.max_enh > 0:
                        self.genes[tf.bound_gene].transcription_rate += self.genes[tf.gene_id].protein_concentration * math.exp(beta*(tf.binding_strength - tf.max_enh - 1))
                    elif tf.max_inh > 0: 
                        self.genes[tf.bound_gene].transcription_rate -= self.genes[tf.gene_id].protein_concentration * math.exp(beta*(tf.binding_strength - tf.max_inh - 1))
                    self.genes[tf.bound_gene].num_bound += 1
                    tf.binding_strength -= 1
                    if tf.binding_strength == 0:
                        new_tf_gene_id = highest_conc_gene_id
                        self.genes[new_tf_gene_id].gene_tf_count += 1
                        self.genes[tf.gene_id].gene_tf_count -= 1
                        tf.gene_id = new_tf_gene_id
                        tf.marker = self.genes[new_tf_gene_id].protein
                        tf.max_enh = -1
                        tf.max_inh = -1
                        tf.bound_gene = -1
                        tf.pos = (0, 0)
                    continue

                # Moving Phase
                move = random.choice([-1, 1]) * random.randint(0, int(self.config["TF_step_size"]))
                tf.pos = ((tf.pos[0] + move) % self.cell_size, (tf.pos[1] + move) % self.cell_size)

                # Binding Phase
                for gene in self.genes:
                    if Math.euclidean_distance(gene.epos, tf.pos) <= self.TF_bind_threshold or Math.euclidean_distance(gene.ipos, tf.pos) <= self.TF_bind_threshold:
                        if Math.euclidean_distance(gene.epos, tf.pos) <= self.TF_bind_threshold:
                            enh_inh_index = 0
                        elif Math.euclidean_distance(gene.ipos, tf.pos) <= self.TF_bind_threshold:
                            enh_inh_index = 1
                        else:
                            if random.random() < 0.5:
                                enh_inh_index = 0
                            else:
                                enh_inh_index = 1
                        binding_strength = self.interaction_matrix[tf.gene_id][gene.id][enh_inh_index]
                        if binding_strength == 0:
                            # unable to bind
                            continue
                        tf.binding_strength = binding_strength
                        if enh_inh_index == 0:
                            tf.max_enh = self.max_enh
                            tf.pos = gene.epos
                        elif enh_inh_index == 1:
                            tf.max_inh = self.max_inh
                            tf.pos = gene.ipos
                        tf.bound_gene = gene.id

            total_concentration = 0
            # Production Phase
            for gene in self.genes:
                delta_protein = self.delta * gene.transcription_rate * gene.protein_concentration
                gene.protein_concentration += delta_protein
                if gene.protein_concentration <= 0.001:
                    gene.protein_concentration = 0.001

                total_concentration += gene.protein_concentration
            conc_txt = ""
            rate_txt = ""
            log_txt = str(cycle) + " "
            each_array = []
            each_rate_array = []
            highest_conc = -1
            for gene in self.genes:
                gene.protein_concentration = gene.protein_concentration/total_concentration
                if gene.protein_concentration > highest_conc:
                    highest_conc = gene.protein_concentration
                    highest_conc_gene_id = gene.id
                conc_txt += "{},".format(round(gene.protein_concentration, 5))
                rate_txt += "{},".format(round(gene.transcription_rate, 5))
                log_txt += str(round(gene.protein_concentration, 5)) + " "
                each_array.append(gene.protein_concentration)
                each_rate_array.append(gene.transcription_rate)
            conc_array.append(each_array)
            rate_array.append(each_rate_array)
            # return conc_array

            # print(log_txt)
            # Files.write(conc_file, conc_txt)
            # Files.write(output_file, rate_txt)
            #
            # Files.lbreak(output_file)
            # Files.lbreak(conc_file)
        return conc_array, rate_array


def return_site_from_dna(dna, start, size):
    if len(dna) - start < size:
        site = dna[start:]
        site += dna[0: size - len(site)]
    else:
        site = dna[start:start+size]
    return site
