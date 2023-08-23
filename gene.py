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
        if len(dna) < 6 or len(self.exon) < 6:
            raise Exception("DNA or EXON too short")

        # regulatory sites length = radical exon length
        reg_length = math.floor(math.sqrt(len(self.exon)))

        reg_site_starting_point = start % len(dna)
        self.enhancer = cell.return_site_from_dna(dna, reg_site_starting_point, reg_length)
        self.inhibitor = cell.return_site_from_dna(dna, reg_site_starting_point + reg_length, reg_length)
        self.protein = self.find_protein(reg_length)

    def find_protein(self, length):
        # print(f"length of exon: {len(self.exon)}, length of reg site: {length}, ", end="")
        exon = self.exon[length * 2:]
        chunk_size = math.floor(len(exon) / length)
        # print(f"length of exon after trim: {len(exon)}, chunk size: {chunk_size}")
        protein = []
        # for each bp in the protein
        for i in range(length):
            if i == length - 1:
                chunk = exon[i * chunk_size:]
            else:
                chunk = exon[i * chunk_size: (i+1) * chunk_size]
            protein.append(self.find_majority(chunk))
        return protein

    def find_majority(self, sequence):
        base_counts = {}

        for i, base in enumerate(sequence):
            if base in base_counts:
                base_counts[base][0] += 1
            else:
                base_counts[base] = [1, i]

        most_common = None

        for base, (count, index) in base_counts.items():
            if most_common is None or most_common[1] < count or (most_common[1] == count and most_common[2] > index):
                most_common = [base, count, index]

        return most_common[0]


