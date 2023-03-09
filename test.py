# class testobject:
# 	def __init__(self):
# 		self.character = 0
#
#
# object_list = []
# changing = []
# for i in range(10):
# 	object_list.append(testobject())
#
# for i in range(5):
# 	changing.append(object_list[i])
#
# changing[3].character = 2
# changing = []
# changing.append(3)
#
# for object in object_list:
# 	print(object.character)
# data = [[1,2,3,45]]
# print(len(data[0]))
# output = [] * len(data[0])
# print(output)

# from scipy import stats
#
# a = [.1, .2, .3, .4, .5, .6]
# b = [.2, .1, .2, .2, .3, .1]
#
# rvalue, pvalue = stats.pearsonr(a, b)
#
# print(rvalue, pvalue)
import random

for i in range(10):
	print(random.randint(0, 1))
