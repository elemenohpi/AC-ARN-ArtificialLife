import math
import cell


class Gene:
    def __init__(self):
        self.promoter = ""
        self.enhancer = []
        self.inhibitor = []
        self.protein = []
        self.exon = []
        self.epos = ()
        self.ipos = ()
        self.transcription_rate = 0.001
        self.num_bound = 0
        self.total_bound_conc = 0
        self.gene_tf_count = 0
        self.id = -1
        self.protein_concentration = 0

    def print(self):
        for bp in self.exon:
            print(bp, end="")
        print("\n")
        print("Enhancer: ", self.enhancer)
        print("Inhibitor: ", self.inhibitor)
        print("Protein: ", self.protein)

    def process_exon(self, dna, start):
        if len(dna) < 4 or len(self.exon) < 4:
            raise Exception("DNA or EXON too short")

        # regulatory sites length = radical exon length
        reg_length = math.floor(math.sqrt(len(self.exon)))

        locator = self.exon[0:reg_length]

        distance = 0
        for bp in locator:
            if bp == "T":
                distance += -1
            elif bp == "G":
                distance += -2
            elif bp == "C":
                distance += +1
            elif bp == "A":
                distance += +2

        reg_site_starting_point = (start + distance) % len(dna)
        self.enhancer = cell.return_site_from_dna(dna, reg_site_starting_point, reg_length)
        self.inhibitor = cell.return_site_from_dna(dna, reg_site_starting_point + reg_length, reg_length)
        self.protein = self.find_protein(reg_length)

    def find_protein(self, length):
        exon = self.exon[length:]
        chunk_size = math.ceil(len(exon) / length)
        protein = []
        # for each bp in the protein
        for i in range(length):
            if i == length - 1:
                chunk = exon[i * chunk_size:]
            else:
                chunk = exon[i * chunk_size: (i+1) * chunk_size]
            protein.append(self.find_majority(chunk))
        return protein

    def find_majority(self, array):
        # index 0, 1, 2, 3 are A, G, C, T respectively
        distribution = {"A": 0, "G": 0, "C": 0, "T": 0}
        for element in array:
            distribution[element] += 1
        majority = max(distribution, key=distribution.get)
        return majority


