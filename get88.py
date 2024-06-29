import numpy as np
import os
import re
import random
import matplotlib.pyplot as plt

NODE = 88
GEN = 100
N = 200
M = 0.0001

north_east_list = []

def output_png( gene, gen, score ):
    x = []
    y = []

    path = 'N_' + str(N) + '_GEN_' + str(GEN)
    if not os.path.isdir(path):
        os.makedirs(path)

    for i in range(NODE):
        x.append(north_east_list[gene[i]][2])
        y.append(north_east_list[gene[i]][1])
    x.append(north_east_list[gene[0]][2])
    y.append(north_east_list[gene[0]][1])

    #plt.text(35.2, 132.5, 'gen = ' +  str(gen) + ', score = ' + str(score))
    plt.suptitle('gen = ' +  str(gen) + ', score = ' + str(score))
    plt.plot(x, y)
    plt.scatter(x, y, marker='o')
    plt.savefig(path + '/' + str(gen) + '.png', format="png")
    plt.cla()
    return

def distance_p2p(x1, y1, x2, y2):
    return np.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))

def distance_all(gene):
    ret = 0.0
    for i in range(NODE):
        x1 = north_east_list[gene[i]][2]
        y1 = north_east_list[gene[i]][1]
        x2 = north_east_list[gene[(i+1)%NODE]][2]
        y2 = north_east_list[gene[(i+1)%NODE]][1]
        ret += distance_p2p( x1, y1, x2, y2 )
    return ret

def crossover(gene_a, gene_b):
    ret_gene_list = []
    r = random.randint(1, NODE - 1)
    gene_a_1st = np.copy(gene_a[0])
    gene_a_2nd = np.copy(gene_a[0])
    gene_a_1st = gene_a_1st[:r]
    gene_a_2nd = gene_a_2nd[r:]
    gene_b_1st = np.copy(gene_b[0])
    gene_b_2nd = np.copy(gene_b[0])
    # gene_b_1stをgene_a_2ndの要素以外のデータからなるリストにする
    # gene_b_2ndをgene_a_1stの要素以外のデータからなるリストにする
    for value in gene_a_2nd:
        gene_b_1st = gene_b_1st[gene_b_1st != value]
    for value in gene_a_1st:
        gene_b_2nd = gene_b_2nd[gene_b_2nd != value]
    #gene_b_1st = gene_b_1st[:r]
    #gene_b_2nd = gene_b_2nd[r:]

    gene_new1 = np.hstack((gene_a_1st, gene_b_2nd))
    gene_new2 = np.hstack((gene_b_1st, gene_a_2nd))
    #if not check_gene(gene_new1):
    #    raise ValueError('gene_new1' + str(gene_new1))
    #if not check_gene(gene_new2):
    #    raise ValueError('gene_new2' + str(gene_new2))
    ret_gene_list.append(gene_new1)
    ret_gene_list.append(gene_new2)

    return ret_gene_list

def check_gene(gene):
    sorted_gene = [x for x in range(NODE)]
    if sorted_gene == sorted(gene):
        return True
    else:
        return False

def mutation(gene):
    r1 = random.randint(1, NODE - 1)
    r2 = random.randint(1, NODE - 1)
    ret_gene = np.copy(gene)
    ret_gene[r1], ret_gene[r2] = ret_gene[r2], ret_gene[r1]
    return ret_gene

with open('Fudasho88.kml', 'r', encoding='utf-8') as fr:
    for line in fr:
        m = re.search('<name>([0-9]+)番', line.strip())
        if m:
            count = int(m.group(1))
        m = re.search('<span>北緯:.*?\(([0-9.]+)\)<\/span>', line.strip())
        if m:
            north = float(m.group(1))
        m = re.search('<span>東経:.*?\(([0-9.]+)\)<\/span>', line.strip())
        if m:
            east  = float(m.group(1))
            north_east_list.append((count, north, east))

# Nだけgeneを生成
gene_list = []
gene_org = np.array([x for x in range(NODE)])
rng = np.random.default_rng()
for _ in range(N):
    gene = np.copy(gene_org)
    rng.shuffle(gene)
    gene_list.append((gene, distance_all(gene)))

for g in range(1, GEN + 1):
    print('GEN = ' + str(g), end=' ')
    # geneを交叉,突然変異
    for x in range(N):
        for y in range(N):
            new_gene_list = crossover(gene_list[x], gene_list[y])
            if random.random() < M:
                i = random.randint(0, 1)
                new_gene_list[i] = mutation(new_gene_list[i])
            for gene in new_gene_list:
                gene_list.append((gene, distance_all(gene)))

    # geneを評価(ソート)
    gene_list.sort(key=lambda x: x[1])
    del gene_list[N:]
    print(gene_list[0][1], end=' ')
    print(gene_list[0][0])
    if (g - 1) % 10 == 0:
        output_png(gene_list[0][0], g, gene_list[0][1])

print(gene_list[0])
output_png(gene_list[0][0], g, gene_list[0][1])
