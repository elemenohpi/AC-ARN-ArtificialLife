class ConfigParser:
    # constructor
    def __init__(self):
        pass

    def auto_generate(self):
        # ToDo
        pass

    def read(self, path):
        config = {}
        try:
            file = open(path, 'r')
        except IOError as e:
            raise Exception(e.strerror)
        
        lines = file.readlines()
        for num, line in enumerate(lines):
            line = line.split("#")[0].strip()
            
            if len(line) == 0:
                # empty line, ignore
                continue
            elif line[0] == "#":
                # comment
                continue
            else:
                tokens = line.split("=")
                if len(tokens) != 2:
                    # illegal length
                    raise Exception("Illegal format on line " + repr(num) + " of " + file)
                else:
                    # legal format
                    config[tokens[0].strip()] = tokens[1].strip()
        
        return config
