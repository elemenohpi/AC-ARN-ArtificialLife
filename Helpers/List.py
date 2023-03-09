# Returns the element that has the most occurrance in a given list
def majority(list):
    dict = {}
    max = []
    for element in list:
        try:
            dict[element] += 1
            if dict[element] > max[1]:
                max[1] = dict[element]
                max[0] = element
        except:
            dict.update({element: 1})
            if max == []:
                max.append(element)
                max.append(1)
    return max[0]

def make_matrix(x, y):
    list = []
    for i in range(x):
        list.append([])
        for _ in range(y):
            list[i].append([])
    return list