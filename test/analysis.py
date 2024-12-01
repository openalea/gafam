from collections import OrderedDict

import numpy as np
import pandas as pd

from openalea.mtg import *


from gafam import analysis
from gafam.data import file, treatment

trtfile = treatment()
df = pd.read_csv(trtfile)
tree_trt = dict(zip(list(df.Tree), list(df.traitement)))


fns=analysis.load_data()

trees = []
trts = []
LA_2018s = []
LA_2019s = []
nb_GU_2018s = []
nb_GU_2019s = []

tree_columns='apple_tree trt LA_2018 LA_2019 nb_GU_2018 nb_GU_2019'
GU_columns='tree trt years GU LeafArea nb_GU nb_leaves nb_flowers nb_fruits'
B_columns='tree trt years'


 GU_trees = []
GU_trts = []
GU_years = []
GU_GU = []
GU_LAs = []
nb_GUs = []
GU_nb_leaves = []
GU_nb_flowers = []
GU_nb_fruits = []
B_trees = []
B_trts = []
B_years = []
B_diameters = []
B_angles = []


for fn in fns:

    g = MTG(fn, has_date=True)

    t = Tree2(g)
    dates = t.dates()
    df = t.dataframe()

    lines = g.property('_line')

    apple_tree = int(g.index(1))
    TRT = tree_trt[apple_tree]

    #print('#'*80)
    #print('Tree ', apple_tree, TRT)



    vs18 = [v for v, d in dates.items() if d =='2018']
    vs19 = [v for v, d in dates.items()
                        if (d=='2019') or
                        (isinstance(lines.get(v), list))
                    ]

    LA_2018 = sum(t.la18(v) for v in vs18) 

    # LA_2019
    LA_2019 = 0
    la19 = g.property('leaf_area_2019')
    for v in vs19:
        la19[v]= la19[v] if v in la19 else sum(la19.get(c,0) for c in g.components(v))
        LA_2019 += la19[v]


    print('LA 2018', LA_2018)
    print ('Tree LA18', g[1]['leaf_area_2018'])
    print('LA 2019', LA_2019)
    print ('Tree LA19', g[1]['leaf_area_2019'])

    labels18 = {}
    nb_GU_18 = 0
    for label in list('ABCDE'):
        labels18[label] = [v for v in vs18 if g.class_name(v) ==label]
        nb_GU_18 += len(labels18[label])
    print(labels18)


    labels19 = {}
    nb_GU_19 = 0
    for label in list('ABCDE'):
        labels19[label] = [v for v in vs19 if g.class_name(v) ==label]
        nb_GU_19 += len(labels19[label])

    print(labels19)

    _la18 = 0
    for label in labels18:
        print('%s18'%label,)
        x = sum(map(t.la18, labels18[label]))
        print(x)
        _la18+=x
    print('LeafSurface 18', _la18)

    _la19 = 0
    for label in labels19:
        print('%s19'%label,)
        x = sum(map(t.la19, labels19[label]))
        print(x)
        _la19+=x
    print('LeafSurface 19', _la19)


    trees.append(apple_tree)
    trts.append(TRT)
    LA_2018s.append(_la18)
    LA_2019s.append(_la19)
    nb_GU_2018s.append(nb_GU_18)
    nb_GU_2019s.append(nb_GU_19)


    ##################
    ## GUs

    _GU_trees = [apple_tree]*10
    _GU_trts = [TRT]*10
    _GU_years = ['2018']*5+['2019']*5
    _GU_GU = list('ABCDE') + list('ABCDE')
    _GU_LAs = [0.]*10
    _nb_GUs = [0]*10
    _GU_nb_leaves = [0]*10
    _GU_nb_flowers = [0]*10
    _GU_nb_fruits = [0]*10


    # Count Number of element/GU ABCDE in 2018 and 2019
    print('#'*80)
    for i, klass in enumerate(list('ABCDE')):
        print('%s18'%klass, len(labels18[klass]) ,)
        _nb_GUs[i] = len(labels18[klass])
        print('%s19'%klass, len(labels19[klass]) ,)
        _nb_GUs[i+5] = len(labels19[klass])
        
    for i, klass in enumerate(list('ABCDE')):
        print('Leaf Area %s18'%klass, sum(map(t.la18, labels18[klass])) ,)
        _GU_LAs[i] = sum(map(t.la18, labels18[klass]))
        print('Leaf Area %s19'%klass, sum(map(t.la19, labels19[klass])) ,)
        _GU_LAs[i+5] = sum(map(t.la19, labels19[klass]))



    # Number of leaves, flowers and fruits (properties have been upscaled)

    nb_leaves = g.property('nb_leaves')
    nb_flowers = g.property('nb_flowers')
    nb_fruits = g.property('nb_fruits')

    i=0
    y=18
    for labels in (labels18, labels19):
        for klass in list('ABCDE'):
        
            vids = labels[klass]
            nbl = sum( nb_leaves.get(vid, 0) for vid in vids )
            nbfl = sum( nb_flowers.get(vid, 0) for vid in vids )
            nbfr = sum( nb_fruits.get(vid, 0) for vid in vids )
            _GU_nb_leaves[i] = nbl
            _GU_nb_flowers[i] = nbfl
            _GU_nb_fruits[i] = nbfr
            if nbl: 
                print('nb leaves %s%d'%(klass, y), nbl ,)
            if nbfl : print('nb flowers %s%d'%(klass, y), nbfl ,)
            if nbfr: print('nb fruits %s%d'%(klass, y), nbfr ,)
            i+=1 

        y=19
            

    GU_trees.extend(_GU_trees)
    GU_trts.extend(_GU_trts)
    GU_years.extend(_GU_years)
    GU_GU.extend(_GU_GU)
    GU_LAs.extend(_GU_LAs)
    nb_GUs.extend(_nb_GUs)
    GU_nb_leaves.extend(_GU_nb_leaves)
    GU_nb_flowers.extend(_GU_nb_flowers)
    GU_nb_fruits.extend(_GU_nb_fruits)

    # Angles and basal diameters
    diameter = g.property('diameter_b')
    diam = g.property('diameter')
    diameter.update(diam)
    angle = g.property('angle')

    bid18 = labels18['B']
    bid19 = labels19['B']
    d18 = [diameter[vid] for vid in bid18 if vid in diameter]
    d19 = [diameter[vid] for vid in bid19 if vid in diameter]
    angles18 = [angle[vid] for vid in bid18 if vid in angle]
    angles19 = [angle[vid] for vid in bid19 if vid in angle]
    d18.sort()
    d19.sort()
    angles18.sort()
    angles19.sort()

    B_trees.append(apple_tree)
    B_trts.append(TRT)
    B_years.append('2018')
    B_diameters.append(d18)
    B_angles.append(angles18)

    if d19 and angles19:
        B_trees.append(apple_tree)
        B_trts.append(TRT)
        B_years.append('2019')
        B_diameters.append(d19)
        B_angles.append(angles19)

    if d18: 
        print (' Diameter(min, max, mean) 18', min(d18), max(d18), np.median(d18) )
    if d19: 
        print (' Diameter(min, max, mean) 19', min(d19), max(d19), np.median(d19) )
    if angles18:
        print (' Angle(min, max, mean) 18', min(angles18), max(angles18), np.median(angles18) )
    if angles19:
        print (' Angle(min, max, mean) 19', min(angles19), max(angles19), np.median(angles19) )


#if max(angles18)
# Mise en forme
#Plant TRT YEAR Branch GUs LA Leaves Flowers Fruits 

#tree_columns='apple_tree trt LA_2018 LA_2019 nb_GU_2018 nb_GU_2019'

tree_df = pd.DataFrame( OrderedDict(
            tree=trees,
            trt=trts,
            LA_2018=LA_2018s,
            LA_2019=LA_2019s,
            nb_GU_2018=nb_GU_2018s,
            nb_GU_2019=nb_GU_2019s,
            ),
            columns=tree_columns.split(' '))

#tree_df.to_csv('tree.csv', sep=';', index=False)

GU_df = pd.DataFrame( OrderedDict(
            tree=GU_trees,
            trt=GU_trts,
            years = GU_years,
            GU = GU_GU,
            LeafArea=GU_LAs,
            nb_GU=nb_GUs,
            nb_leaves=GU_nb_leaves,
            nb_flowers=GU_nb_flowers,
            nb_fruits=GU_nb_fruits,
            ),
            columns=GU_columns.split(' '))
#GU_df.to_csv('GU.csv', sep=';', index=False)

B_df = pd.DataFrame( OrderedDict(
            tree=B_trees,
            trt=B_trts,
            years = B_years,
            diameters=B_diameters,
            angles = B_angles
            ),
            columns=B_columns.split(' '))

#B_df.to_csv('B_dia_angle.csv', sep=';', index=False)

def record(tree, trt, year, data, i):
    l = [str(tree[i]), str(trt[i]), str(year[i])]
    l.extend(map(str, data[i]))
    return ';'.join(l)

def write_angle(fn='B_angle.csv'):
    content = []
    for i in range(len(B_trees)):
        content.append(record(B_trees, B_trts, B_years, B_angles, i))
    data = '\n'.join(content)
    f = open(fn,'w')
    f.write(data)
    f.close()

def write_diameter(fn='B_diams.csv'):
    content = []
    for i in range(len(B_trees)):
        content.append(record(B_trees, B_trts, B_years, B_diameters, i))
    data = '\n'.join(content)
    f = open(fn,'w')
    f.write(data)
    f.close()

