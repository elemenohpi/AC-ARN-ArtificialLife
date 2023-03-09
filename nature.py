from cell import Cell


class Nature:
    
    def __init__(self, config):
        self.config = config
        try:
            self.organism_count = int(config["organism_count"])
            self.cycles = int(config["cycles"])
        except KeyError as e:
            raise KeyError(e, "Configuration file missing some variables")

        self.organisms = []

    def create_organisms(self):
        # gene_length = 0
        for _ in range(self.organism_count):
            c = Cell(self.config)
            c.init_genetic_marker()
            # c.identify_genes()
            # gene_length += len(c.genes)
            self.organisms.append(c)
        # print(gene_length / self.organism_count)
        # exit()

    def evolve(self):
        pass
