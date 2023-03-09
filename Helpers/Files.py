
def truncate(file, create=True):
    try:
        f = open(file)
        f.truncate()
    except IOError:
        if create:
            f = open(file, "w+")
            f.close()
            return 1
        return 0
    return 1

def write_line(file, str, create=True):
    try:
        f = open(file, "a")
        f.write(str)
        f.write("\n")
    except IOError:
        if create:
            f = open(file, "w+")
            f.write(str)
            f.write("\n")
            f.close()
            return 1
        return 0
    finally:
        f.close()
    return 1

def write(file, str, create=True):
    try:
        f = open(file, "a")
        f.write(str)
    except IOError:
        if create:
            f = open(file, "w+")
            f.write(str)
            f.close()
            return 1
        return 0
    finally:
        f.close()
    return 1

def lbreak(file, create=False):
    try:
        f = open(file, "a")
        f.write("\n")
    except IOError:
        if create:
            f = open(file, "w+")
            f.write("\n")
            f.close()
            return 1
        return 0
    finally:
        f.close()
    return 1
    pass
