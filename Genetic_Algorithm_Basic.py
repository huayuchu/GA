# coding=utf-8

"""
目标：找到在函数中最大的点
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


DNA_SIZE = 10            # DNA长度
POP_SIZE = 100           # population size 种群大小
CROSS_RATE = 0.8         # mating probability (DNA crossover) 交叉率
MUTATION_RATE = 0.003    # mutation probability 变异概率
N_GENERATIONS = 200      # 代数(迭代次数)
X_BOUND = [0, 5]         # x upper and lower bounds x的上限和下限

# 环境函数
def F(x): return np.sin(10*x)*x + np.cos(2*x)*x     # to find the maximum of this function


# 查找非0适应度值用于选择  find non-zero fitness for selection
def get_fitness(pred): return pred + 1e-3 - np.min(pred)


# convert binary DNA to decimal and normalize it to a range(0, 5)
def translateDNA(pop):
    value = pop.dot(2 ** np.arange(DNA_SIZE)[::-1]) / float(2**DNA_SIZE-1) * X_BOUND[1]
    # print(2 ** np.arange(DNA_SIZE)[::-1])
    # print(value)
    return value


def select(pop, fitness):    # nature selection wrt pop's fitness
    idx = np.random.choice(np.arange(POP_SIZE), size=POP_SIZE, replace=True,
                           p=fitness/fitness.sum())
    # print('idx: \t', idx)
    return pop[idx]


def crossover(parent, pop):     # 交叉
    if np.random.rand() < CROSS_RATE:
        i_ = np.random.randint(0, POP_SIZE, size=1)                             # 在种群里另外选一个交叉的对象 select another individual from pop
        cross_points = np.random.randint(0, 2, size=DNA_SIZE).astype(np.bool)   # 定位交叉的位置 choose crossover points
        parent[cross_points] = pop[i_, cross_points]                            # 开始交叉，其实是替换掉对应位置的值 mating and produce one child
    return parent


def mutate(child):
    for point in range(DNA_SIZE):
        if np.random.rand() < MUTATION_RATE:
            child[point] = 1 if child[point] == 0 else 0
    return child


pop = np.random.randint(2, size=(POP_SIZE, DNA_SIZE))   # initialize the pop DNA
# print(pop)

plt.ion()       # something about plotting
x = np.linspace(*X_BOUND, 200)
# print(x)
plt.plot(x, F(x))

for _ in range(N_GENERATIONS):
    F_values = F(translateDNA(pop))    # 计算各条DNA转换计算之后在F函数中的值 compute function value by extracting DNA
    # print(F_values)

    # something about plotting
    if 'sca' in globals(): sca.remove()
    sca = plt.scatter(translateDNA(pop), F_values, s=200, lw=0, c='red', alpha=0.5); plt.pause(0.05)

    # 进化部分 GA part (evolution)
    fitness = get_fitness(F_values)
    # print(fitness)
    print("Most fitted DNA: ", pop[np.argmax(fitness), :])
    pop = select(pop, fitness) # fitness越大就越有可能被选到，这里整个种群都要进行选择，每一个都会尽可能选择最大fitness的那个DNA
    # print(pop)
    pop_copy = pop.copy()
    for parent in pop:
        child = crossover(parent, pop_copy)
        child = mutate(child)
        parent[:] = child       # parent is replaced by its child

plt.ioff(); plt.show()
