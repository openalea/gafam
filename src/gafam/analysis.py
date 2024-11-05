from math import pi
from pathlib import Path
from collections import OrderedDict

import pandas as pd

from openalea.mtg import *
from openalea.mtg import algo, traversal

from . import data

def load_data():
    "Return a list of MTGs"
    return data.files()

def missing_components(g):
    missings = [v for v in g.vertices(scale=2) if g.nb_components(v) == 0]
    return missings

def check_miss(g):
    m = missing_components(g)
    anchors = [v-1 for v in m]
    #Previous element is either of scale 1 or scale 3
    # If scale is 1, then the main axe is not decomposed: impossible
    # if scale is 2, then 2 successives branches are not decomposed: error?
    # Finally, if scale is 3, this is the anchor where we want to add a ghost node
    scale2 = list(a+1 for a in anchors if g.scale(a)==2)
    return scale2

def check_roots(g, scale=2):
    return [v for v in g.vertices_iter(scale=scale) if v!=scale and g.parent(v) is None]


# Check voir test

def debug():
    fns = load_data()
    for fn in fns:
        print('\n'+'#'*20)
        print(fn)
        g=MTG(fn, has_date=True)
        g.display(max_scale=3)
        #x=raw_input('Yes?')

def tree(g):
    """
    Extract information at tree scale:
        * apple_tree:
            Apple tree number
        * trt:
            treatment or NCI category
        * NCI:
            Neighbourhood crowding index value
        * TSA_2018:
            Trunk section area:
            $\pi*diameter_b_2018^2/4$ (cm2)
        * TSA_2019:
            Trunk section area: $\pi*diameter_b^2/4$ (cm2)
        * L_shoot_V/F_2018/2019:
            Number of vegetative of floral long shoots (>=5cm) in 2018 and 2019
        * S_shoot_V/F_2018/2019:
            Number of vegetative of floral short shoots (<5cm) in 2018 and 2019
        * nb_V_2018/2019:
            Number of vegetative buds in 2018 and 2019
        * nb_F_2018/2019:
            Number of floral buds in 2018 and 2019
        * nb_L_2018/2019:
            Number of latent buds in 2018 and 2019
        * LA_2018/2019:
            Leaf area (cm2) in 2018 and 2019

    """
    pass


###############################################################################

def error_date(g):
    """ Extract all the same vertices at dat 2018 an 2019
    that are not the same  (define with *)

    TODO: Remove
    """

    date = g.property('realisation')

    vids = []
    for v in g.roots_iter(3):
        if g.class_name(v) == g.class_name(v-1):
            if date.get(v) == 2019 and date.get(v-1)==2018:
                vids.append(v)

    roots = g.roots_iter(3)
    roots.next()

    # We search all the patterns that are not A<A/U
    for v in roots:
        p1,p2 = v-1, v-2
        if (g.scale(p1) != g.scale(p2)) or (g.scale(v) != g.scale(p1)+1):
            vids.append(v)

    vids = list(set(vids))
    vids.sort()

    return [g.label(v) for v in vids]

def errors_to_csv(filename='errors.csv'):
    fns = load_data()
    dfs = []
    for fn in fns:
        try:
            g = MTG(fn, has_date=True)
            name = fn.name
            errors = error_date(g)
            df = pd.DataFrame(dict(name=name, vertex=errors))
            dfs.append(df)
        except:
            print('Error: {}'.format(fn))
            continue

    df = pd.concat(dfs, ignore_index=True, axis=0)
    #df.to_csv(filename, sep=';', index=False)
    return df


class A(object):
    def __init__(self, g):
        self.g = g
        self.is_dynamic()

    def __call__(self):
        return self.dataframe()

    def is_dynamic(self):
        """ Check if the MTG is dynamic or not. """
        g = self.g
        ps=g.properties()
        self.dynamic = False
        for p in ps:
            has_list = set(isinstance(x,list) for x in ps[p].values())
            if True in has_list:
                self.dynamic = True
                break

    def att(self, year, name, v, default=None):
        g = self.g
        p = g.property(name).get(v, default)
        if isinstance(p, list):
            p = dict(p)
            return p[year]
        else:
            return p

    def is_veg(self, year, v):
        cat = self.att(year, 'type', v, '')
        return cat.lower() == 'vg'
    def is_flo(self, year, v):
        cat = self.att(year, 'type', v, '')
        return cat.lower() == 'fl'
    def is_latent(self, year, v):
        cat = self.att(year, 'type', v, '')
        return cat.lower() == 'bl'
    def long(self, year, v):
        length = self.att(year, 'length', v, 0)
        return length >= 5
    def short(self, year, v):
        length = self.att(year, 'length', v, 0)
        return 0 <= length < 5
    def la18(self, v):
        return self.att('2018', 'leaf_area_2018', v, 0)
    def la19(self, v):
        return self.att('2019', 'leaf_area_2019', v, 0)

class Tree(A):

    def dates(self):
        g = self.g
        dates = g.property('Date')
        
        _dates = dict()
        for vid in dates:
            date = dates[vid]
            if g.scale(vid) == 2:
                _dates[vid] = date
            elif g.scale(vid) == 3:
                cid = g.complex(vid)
                _dates[cid] = date
            else:
                print('WARNING on dates', vid)
        self._dates = _dates
        return _dates


    def dataframe(self):
        """ Extract Tree level information
        return a dataframe

        Algorithms
        ----------
            - apple_tree:
                Apple tree number
            - trt:
                treatment or NCI category
            - NCI:
                Neighbourhood crowding index value
            - TSA_2018:
                Trunk section area: pi*diameter_b_2018**2/4 (cm2)
            - TSA_2019:
                Trunk section area: pi*diameter_b**2/4 (cm2)
            - L_shoot_V/F_2018/2019:
                Number of vegetative of floral long shoots (>=5cm) in 2018 and 2019
            - S_shoot_V/F_2018/2019:
                Number of vegetative of floral short shoots (<5cm) in 2018 and 2019
            - nb_V_2018/2019:
                Number of vegetative buds in 2018 and 2019
            - nb_F_2018/2019:
                Number of floral buds in 2018 and 2019
            - nb_L_2018/2019:
                Number of latent buds in 2018 and 2019
            - LA_2018/2019:
                Leaf area (cm2) in 2018 and 2019
            - Elongation :
                longueur du tronc (A1+A2+A3+graft point) / mean diameter (base A1 + apex A3)/2
            - Tapper :
                (diameter base (A1) - diameter apex (A3) / longueur trunk

        """
        g = self.g

        dates = self.dates()
        

        apple_tree = g.index(1)
        trt = 'ac'
        NCI = ''

        pnames = g.property_names()

        # Trunk section area: pi*diameter_b_2018**2/4 (cm2)
        A1 = g.node(2)

        if 'diameter_b2018' in g[2]:
            TSA_2018 = pi * (A1.diameter_b2018)**2 / 4.
        else:
            TSA_2018 = pi * (A1.diameter_2018)**2 / 4.

        #Trunk section area: pi*diameter_b**2/4 (cm2)
        TSA_2019 = pi * (A1.diameter_b)**2 / 4.

        #  Number of vegetative of floral long shoots (>=5cm) in 2018 and 2019
        # Colonne type pour le type et length pour la longueur pour 2018
        # mes informations sont a l'echelle B dans la plupart des cas
        # (peut etre des cas a l'echelle S). Pour 2019 a l'echelle S

        # cat = g.property('type')
        # length = g.property('length')
        # leaf_area18 = g.property('leaf_area_2018')
        # leaf_area19 = g.property('leaf_area_2019')

        # 'BL', 'EX', 'FL', 'Fl', 'VG'
        # is_veg = lambda x: cat.get(x, '').lower() =='vg'
        # is_flo = lambda x: cat.get(x, '').lower() =='fl'
        # is_latent = lambda x: cat.get(x, '').lower() =='bl'
        # long = lambda x: length.get(x,0) >= 5
        # short = lambda x: 0<= length.get(x,-1) < 5
        # la18 = lambda x: leaf_area18.get(x,0)
        # la19 = lambda x: leaf_area19.get(x,0)

        #dates = g.property('Date')
        lines = g.property('_line')
        vs18 = [v for v, d in dates.items() if d =='2018']
        vs19 = [v for v, d in dates.items()
                    if (d=='2019') or
                       (isinstance(lines.get(v), list))
               ]
        L_shoot_V_2018 = sum(1 for v in vs18
                                if self.is_veg('2018', v) and
                                   self.long('2018', v)
                            )
        L_shoot_F_2018 = sum(1 for v in vs18
                                if self.is_flo('2018', v) and
                                   self.long('2018', v)
                            )
        L_shoot_V_2019 = sum(1 for v in vs19
                                if self.is_veg('2019', v) and
                                   self.long('2019', v)
                            )
        L_shoot_F_2019 = sum(1 for v in vs19
                                if self.is_flo('2019', v) and
                                   self.long('2019', v)
                            )

        #  Number of vegetative of floral short shoots (<5cm) in 2018 and 2019
        S_shoot_V_2018 = sum(1 for v in vs18
                                if self.is_veg('2018', v) and
                                   self.short('2018', v)
                            )
        S_shoot_F_2018 = sum(1 for v in vs18
                                if self.is_flo('2018', v) and
                                   self.short('2018', v)
                            )
        S_shoot_V_2019 = sum(1 for v in vs19
                                if self.is_veg('2019', v) and
                                   self.short('2019', v)
                            )
        S_shoot_F_2019 = sum(1 for v in vs19
                                if self.is_flo('2019', v) and
                                   self.short('2019', v)
                            )

        # Number of vegetative buds in 2018 and 2019
        nb_V_2018 = sum(1 for v in vs18 if self.is_veg('2018', v))
        nb_V_2019 = sum(1 for v in vs19 if self.is_veg('2019', v))

        # Number of floral buds in 2018 and 2019
        nb_F_2018 = sum(1 for v in vs18 if self.is_flo('2018', v))
        nb_F_2019 = sum(1 for v in vs19 if self.is_flo('2019', v))

        # Number of latent buds in 2018 and 2019
        # Bug: latent bugs have not always a date
        nb_L_2018 = sum(1 for v in vs18 if self.is_latent('2018', v))
        nb_L_2019 = sum(1 for v in vs19 if self.is_latent('2019',v))

        # Leaf area (cm2) in 2018 and 2019
        #
        LA_2018 = sum(self.la18(v) for v in vs18)
        LA_2019 = sum(self.la19(v) for v in vs19)

        # longueur du tronc length graft point +(A1+A2+A3) / diametre base  (A1)
        trunk = g.Trunk(2)
        assert(len(trunk) == 4)

        A1 = g.node(trunk[0])
        A3 = g.node(trunk[2])
        length = g.property('length')
        trunk_len = [length.get(v, 0) for v in trunk]
        # A4 = trunk[-1]
        # if A4 not in length:
        #     first_S = next(g.component_roots_iter(A4))
        #     A4_shoots=list(algo.axis(g, first_S, RestrictedTo='SameComplex'))
        #     trunk_len[-1] = sum(length.get(v, 0) for v in A4_shoots)
        total_length = sum(trunk_len[:-1]) + A1.length_graftpoint_ram

        # Mean diameter: (A1.diameter_b + A3.diameter_a) /2
        mean_diameter = (A1.diameter_b + A3.diameter_a) / 2.
        Elongation = total_length / mean_diameter
        # (diametre base (A1) - diametre sommet) / longueur total

        A3 = g.node(trunk[2])
        Tapper = (A1.diameter_b - A3.diameter_a) / total_length

        # number of node on the trunk
        nb_Trunk_S = len(g.Trunk (2, Scale=3))

        df = pd.DataFrame( OrderedDict(
            apple_tree=[apple_tree],
            trt=[trt],
            NCI=[NCI],
            TSA_2018=[TSA_2018],
            TSA_2019=[TSA_2019],
            L_shoot_V_2018=[L_shoot_V_2018],
            L_shoot_V_2019=[L_shoot_V_2019],
            L_shoot_F_2018=[L_shoot_F_2018],
            L_shoot_F_2019=[L_shoot_F_2019],
            S_shoot_V_2018=[S_shoot_V_2018],
            S_shoot_V_2019=[S_shoot_V_2019],
            S_shoot_F_2018=[S_shoot_F_2018],
            S_shoot_F_2019=[S_shoot_F_2019],
            nb_V_2018=[nb_V_2018],
            nb_V_2019=[nb_V_2019],
            nb_F_2018=[nb_F_2018],
            nb_F_2019=[nb_F_2019],
            nb_L_2018=[nb_L_2018],
            nb_L_2019=[nb_L_2019],
            LA_2018=[LA_2018],
            LA_2019=[LA_2019],
            Elongation=[Elongation],
            Tapper=[Tapper],
            nb_Trunk_S=[nb_Trunk_S],
            ),
            columns='apple_tree trt NCI TSA_2018 TSA_2019 L_shoot_V_2018 L_shoot_V_2019 '
            'L_shoot_F_2018 L_shoot_F_2019 S_shoot_V_2018 S_shoot_V_2019 S_shoot_F_2018 S_shoot_F_2019 '
            'nb_V_2018 nb_V_2019 nb_F_2018 nb_F_2019 nb_L_2018 nb_L_2019 LA_2018 LA_2019 '
            'Elongation Tapper '.split(' '))
        return df

class Branches(Tree):

    def heights(self):
        g = self.g
        v = next(g.component_roots_at_scale_iter(0, scale=3))
        height = dict()
        for v in traversal.pre_order2(g, v):
            height[v] = height.get(g.parent(v), -1)+1

        return height


    def diameter(self, v, year='2018', d18=None, d19=None):
        g = self.g
        dates = self._dates
        default = 0.3
        if year == '2018':
            if v in d18:
                return d18[v]
            elif dates.get(v) == '2019':
                return 0
            else:
                return default
        else:
            return d19.get(v, default)

    def length(self, v, y='2018'):
        g = self.g
        l0 =  self.att(y,'length',v,0)
        if l0 == 0:
            l0 = max(self.att(y,'length', b, 0) for b in g.components(v))
        return l0


    def branch_length(self, v, is_19=False):
        g = self.g
        dates = self._dates
        y= '2019' if is_19 else '2018'
        if is_19:
            bl = sum( self.length(v,'2019') for b in g.Axis(v) )                
        else:
            bl = sum( self.length(v,'2018') for b in g.Axis(v) if dates[b] != '2019')
        return bl                

    def dataframe(self):
        """
        - BSA_2018/2019:
            Branch section area: pi*diameter_b**2/4 (cm2) or diameter if diameter_b has no value
            Information  at scale +B in column diameter_b or diameter if no diamter_b
        - B_dist:
            branch distance on the trunk (0 base of the tree 1 apex of the tree)
        - B_lgth_2018/2019:
            branch length (cm) in 2018 (2017+2018) and 2019 (2017+2018+2019)
        - L_shoot_V/F_2018/2019:
            Number of vegetative of floral long shoots (>=5cm) on the for in 2018 and 2019
        - S_shoot_V/F_2018/2019:
            Number of vegetative of floral short shoots (<5cm) on the branch in 2018 and 2019
        - nb_V_2018/2019:
            Number of vegetative buds on the branch in 2018 and 2019
        - nb_F_2018/2019:
            Number of floral buds on the branch in 2018 and 2019
        - nb_L_2018/2019:
            Number of latent buds on the branch in 2018 and 2019
        - LA_2018/2019:
            Leaf area (cm2) in 2018 and 2019
        - Elongation :
            branch length / base branch diameter 
        - Tapper :
            (base branch diameter (2017 ou 2018) - top diameter (2019)) / branch total length

        """
        #tree = super(Branches, self).dataframe()
        
        g = self.g
        dates = self.dates()

        apple_tree = g.index(1)

        # Extract branches first vertex of each branch beared by the trunk
        trunk = g.Trunk(2)
        brs = [b for v in trunk for b in g.Sons(v, EdgeType='+')]

        vtrunk = g.Trunk(3)
        #vbrs = [b for v in vtrunk for b in g.Sons(v, EdgeType='+')]

        anchors = dict((b, b-1) for b in brs)
        heights = self.heights()

        # Branch Heights
        tip_height = float(heights[vtrunk[-1]])
        b_dists = [(heights[anchors[v]]+1)/tip_height for v in brs]

        # TODO: Missing data for diameter
        default_diameter = 0.3
        

        # Compute diameter
        d18 = g.property('diameter_b2018').copy()
        d18.update(g.property('diameter_2018'))
        d19 = g.property('diameter_b').copy()
        d19.update(g.property('diameter'))

        BSA_2018 = [self.diameter(b, '2018', d18, d19) for b in brs]
        BSA_2019 = [self.diameter(b, '2019', d18, d19) for b in brs]
        
        B_pos = [(heights[anchors[b]]+1) for b in brs]
        B_dist = b_dists
        
        # WIP
        length = g.property('length')
        B_lgth_2018 = [self.branch_length(b, is_19=False) for b in brs]
        B_lgth_2019 = [self.branch_length(b, is_19=True) for b in brs]

        # TODO
        lines = g.property('_line')
        def isdyn(v):
            return isinstance(lines.get(v), list)

        def br18(b):
            axis = g.Axis(b)
            comps = [c for v in axis for c in g.components(v) if dates.get(c)=='2018']
            branch = [a for a in axis if dates.get(a)=='2018']
            branch.extend(comps)
            return branch

        def br19(b):
            
            axis = g.Axis(b)
            comps = [c for v in axis for c in g.components(v) 
                    if dates.get(c)=='2019' or 
                    (dates.get(c) =='2018' and isdyn(v))]
            branch = [a for a in axis if dates.get(a)=='2019']
            branch.extend(comps)
            return branch

        
        L_shoot_V_2018 = [sum(1 for v in br18(b)
                                if self.is_veg('2018', v) and
                                   self.long('2018', v)
                            ) for b in brs]
        L_shoot_F_2018 = [sum(1 for v in br18(b)
                                if self.is_flo('2018', v) and
                                   self.long('2018', v)
                            ) for b in brs]
        L_shoot_V_2019 = [sum(1 for v in br19(b)
                                if self.is_veg('2019', v) and
                                   self.long('2019', v)
                            ) for b in brs]
        L_shoot_F_2019 = [sum(1 for v in br19(b)
                                if self.is_flo('2019', v) and
                                   self.long('2019', v)
                            ) for b in brs]


        #  Number of vegetative of floral short shoots (<5cm) in 2018 and 2019
        S_shoot_V_2018 = [sum(1 for v in br18(b)
                                if self.is_veg('2018', v) and
                                   self.short('2018', v)
                            ) for b in brs]
        S_shoot_F_2018 = [sum(1 for v in br18(b)
                                if self.is_flo('2018', v) and
                                   self.short('2018', v)
                            ) for b in brs]
        S_shoot_V_2019 = [sum(1 for v in br19(b)
                                if self.is_veg('2019', v) and
                                   self.short('2019', v)
                            ) for b in brs]
        S_shoot_F_2019 = [sum(1 for v in br19(b)
                                if self.is_flo('2019', v) and
                                   self.short('2019', v)
                            ) for b in brs]

        # Number of vegetative buds in 2018 and 2019
        nb_V_2018 = [sum(1 for v in br18(b) if self.is_veg('2018', v))
                        for b in brs ]
        nb_V_2019 = [sum(1 for v in br19(b) if self.is_veg('2019', v))
                        for b in brs ]

        # Number of floral buds in 2018 and 2019
        nb_F_2018 = [sum(1 for v in br18(b) if self.is_flo('2018', v))
                        for b in brs ]
        nb_F_2019 = [sum(1 for v in br19(b) if self.is_flo('2019', v))
                        for b in brs ]

        # Number of latent buds in 2018 and 2019
        # Bug: latent bugs have not always a date
        nb_L_2018 = [sum(1 for v in br18(b) if self.is_latent('2018', v))
                        for b in brs ]
        nb_L_2019 = [sum(1 for v in br19(b) if self.is_latent('2019',v))
                        for b in brs ]

        # Leaf area (cm2) in 2018 and 2019
        #
        LA_2018 = [sum(self.la18(v) for v in br18(b)) for b in brs]
        LA_2019 = [sum(self.la19(v) for v in br19(b)) for b in brs]

        Elongation = [(BSA_2019[i] / B_lgth_2019[i]) if B_lgth_2019[i] else 0. for i in range(len(brs))]
        Tapper = []
        n= len(brs)

        columns= """
        apple_tree
        BSA_2018
        BSA_2019
        B_pos
        B_dist
        B_lgth_2018
        B_lgth_2019
        L_shoot_V_2018
        L_shoot_V_2018
        L_shoot_F_2019
        L_shoot_F_2019
        S_shoot_V_2018
        S_shoot_V_2018
        S_shoot_F_2019
        S_shoot_F_2019
        nb_V_2018
        nb_V_2019
        nb_F_2018
        nb_F_2019
        nb_L_2018
        nb_L_2019
        LA_2018
        LA_2019
        Elongation
        Tapper
        """.split()
        df = pd.DataFrame( OrderedDict(
            apple_tree=[apple_tree]*n,
            BSA_2018=BSA_2018,
            BSA_2019=BSA_2019,
            B_pos=B_pos,
            B_dist = B_dist,
            B_lgth_2018=B_lgth_2018,
            B_lgth_2019=B_lgth_2019,
            L_shoot_V_2018=L_shoot_V_2018,
            L_shoot_V_2019=L_shoot_V_2019,
            L_shoot_F_2018=L_shoot_F_2018,
            L_shoot_F_2019=L_shoot_F_2019,
            S_shoot_V_2018=S_shoot_V_2018,
            S_shoot_V_2019=S_shoot_V_2019,
            S_shoot_F_2018=S_shoot_F_2018,
            S_shoot_F_2019=S_shoot_F_2019,
            nb_V_2018=nb_V_2018,
            nb_V_2019=nb_V_2019,
            nb_F_2018=nb_F_2018,
            nb_F_2019=nb_F_2019,
            nb_L_2018=nb_L_2018,
            nb_L_2019=nb_L_2019,
            LA_2018=LA_2018,
            LA_2019=LA_2019,
            Elongation=Elongation,
            Tapper=[0]*n,
            ),
            columns=columns)
        return df


def forest(fns=None, save=False):
    if fns is None:
        fns=load_data()

    dfs = []
    errors=[]
    for fn in fns:
        try:
            print('#'*80)
            print("MTG is ", fn)
            g=MTG(fn, has_date=True)
            df = Tree(g)()
            dfs.append(df)
        except:
            print('#'*80)
            print('Error with file {}'.format(fn))
            print('#'*80)
            print()
            errors.append(fn)
            continue

    all_df = pd.concat(dfs, axis=0, ignore_index=True)
    if save:
        all_df.to_csv('trees.csv', sep=';', index=False)


    return all_df, errors

def branches(fns=None, save=False):
    if fns is None:
        fns=load_data()

    dfs = []
    errors=[]
    for fn in fns:
        try:
            print('#'*80)
            print("MTG is ", fn)
            g=MTG(fn, has_date=True)
            df = Branches(g)()
            dfs.append(df)
        except:
            print('#'*80)
            print('Error with file {}'.format(fn))
            print('#'*80)
            print()
            errors.append(fn)
            continue

    all_df = pd.concat(dfs, axis=0, ignore_index=True)
    if save:
        all_df.to_csv('branches.csv', sep=';', index=False)


    return all_df, errors

# Extract sequences sous forme de liste de liste
#  longueur totale, 

#################################################################
# New analysis

class Tree2(Tree):
    
    def upscale(self, property_name, operator=sum):
        """ Upscale properties from lower scale to upper

        Parameters:
            - property_name : str
                MTG property
            - operator: function
                the upscaling function (default : sum) 
        """
        g = self.g
        max_scale = g.max_scale()

        prop = g.property(property_name)

        for vid in g.vertices(scale=max_scale-1):
            if (vid not in prop):
                try:
                    prop[vid]=operator(prop[v] for v in g.components(vid) if v in prop)
                except:
                    prop[vid]=operator(self.att('2019',property_name,v) for v in g.components(vid) if v in prop)

    def is_flo_19(self, vid):
        nb = sum(1 for v in self.g.components(vid) if self.is_flo('2019',v))
        return bool(nb)

    def is_veg_19(self, vid):
        nb = sum(1 for v in self.g.components(vid) if self.is_veg('2019',v))
        return bool(nb)
    
    def is_latent_19(self, vid):
        nb = sum(1 for v in self.g.components(vid) if self.is_latent('2019',v))
        return bool(nb)

    def long_19(self, vid):
        nb = sum(1 for v in self.g.components(vid) if self.long('2019',v))
        return bool(nb)
        
    def short_19(self, vid):
        nb = sum(1 for v in self.g.components(vid) if self.short('2019',v))
        return bool(nb)

    def dataframe(self):
        """ Extract Tree level information
        return a dataframe

        Algorithms
        ----------
            - apple_tree: 
                Apple tree number
            - LA_2018/2019 : 
                sum of the leaf area 
            - GU_len_2018/2019 : 
                sum of the GU length
            - nb_GU_2018/2019:
                count of the GU
            - nb_FL_2018/2019: (nb_F)
                Number of Floral shoots in 2018 and 2019
            - nb_VG_2018/2019: (nb_V)
                Number of vegetative buds in 2018 and 2019
            - nb_BL_2018/2019:
                Number of latent buds in 2018 and 2019
            - nb_EX_2018/2019:
                Number of latent buds in 2018 and 2019
            - mean_angle_2018/2019:
                Mean insertion angle
            - mean_leaves_2018/2019:
                mean value for all GU
            - mean_fruits_2018/2019: idem
            - mean_flowers_2018/2019: idem


        """
        g = self.g

        dates = self.dates()
        
        self.upscale('length', operator=sum)
        self.upscale('leaf_area_2019', operator=sum)
        self.upscale('nb_leaves', operator=sum)
        self.upscale('nb_flowers', operator=sum)
        self.upscale('nb_fruits', operator=sum)

        apple_tree = g.index(1)
        trt = 'ac'
        NCI = ''

        pnames = g.property_names()

        # Trunk section area: pi*diameter_b_2018**2/4 (cm2)
        A1 = g.node(2)

        if 'diameter_b2018' in g[2]:
            TSA_2018 = pi * (A1.diameter_b2018)**2 / 4.
        else:
            TSA_2018 = pi * (A1.diameter_2018)**2 / 4.

        #Trunk section area: pi*diameter_b**2/4 (cm2)
        TSA_2019 = pi * (A1.diameter_b)**2 / 4.

        lines = g.property('_line')
        vs18 = [v for v, d in dates.items() if d =='2018']
        vs19 = [v for v, d in dates.items()
                    if (d =='2019') or
                       (isinstance(lines.get(v), list))
               ]
        # Leaf Area
        LA_2018 = sum(self.la18(v) for v in vs18) 
        LA_2019 = sum(self.la19(v) for v in vs19) 
        
        ######
        L_shoot_V_2018 = sum(1 for v in vs18
                                if self.is_veg('2018', v) and
                                   self.long('2018', v)
                            )
        L_shoot_F_2018 = sum(1 for v in vs18
                                if self.is_flo('2018', v) and
                                   self.long('2018', v)
                            )
        L_shoot_V_2019 = sum(1 for v in vs19
                                if self.is_veg_19(v) and
                                   self.long_19(v)
                            )
        L_shoot_F_2019 = sum(1 for v in vs19
                                if self.is_flo_19(v) and
                                   self.long_19(v)
                            )

        #  Number of vegetative of floral short shoots (<5cm) in 2018 and 2019
        S_shoot_V_2018 = sum(1 for v in vs18
                                if self.is_veg('2018', v) and
                                   self.short('2018', v)
                            )
        S_shoot_F_2018 = sum(1 for v in vs18
                                if self.is_flo('2018', v) and
                                   self.short('2018', v)
                            )
        S_shoot_V_2019 = sum(1 for v in vs19
                                if self.is_veg_19(v) and
                                   self.short_19(v)
                            )
        S_shoot_F_2019 = sum(1 for v in vs19
                                if self.is_flo_19(v) and
                                   self.short_19(v)
                            )

        # Number of vegetative buds in 2018 and 2019
        nb_V_2018 = sum(1 for v in vs18 if self.is_veg('2018', v))
        nb_V_2019 = sum(1 for v in vs19 if self.is_veg_19(v))

        # Compute nb of inflorescence

        nb_F_2018 = sum(1 for v in vs18 if self.is_flo('2018', v))
        nb_F_2019 = sum(1 for v in vs19 if self.is_flo_19(v))

        # Number of latent buds in 2018 and 2019
        # Bug: latent bugs have not always a date
        nb_L_2018 = sum(1 for v in vs18 if self.is_latent('2018', v))
        nb_L_2019 = sum(1 for v in vs19 if self.is_latent_19(v))

        # longueur du tronc length graft point +(A1+A2+A3) / diametre base  (A1)
        trunk = g.Trunk(2)
        assert(len(trunk) == 4)

        A1 = g.node(trunk[0])
        A3 = g.node(trunk[2])
        length = g.property('length')
        trunk_len = [length.get(v, 0) for v in trunk]
        total_length = sum(trunk_len[:-1]) + A1.length_graftpoint_ram

        # Mean diameter: (A1.diameter_b + A3.diameter_a) /2
        mean_diameter = (A1.diameter_b + A3.diameter_a) / 2.
        Elongation = total_length / mean_diameter
        # (diametre base (A1) - diametre sommet) / longueur total

        A3 = g.node(trunk[2])
        Tapper = (A1.diameter_b - A3.diameter_a) / total_length

        # number of node on the trunk
        nb_Trunk_S = len(g.Trunk (2, Scale=3))

        df = pd.DataFrame( OrderedDict(
            apple_tree=[apple_tree],
            trt=[trt],
            NCI=[NCI],
            TSA_2018=[TSA_2018],
            TSA_2019=[TSA_2019],
            L_shoot_V_2018=[L_shoot_V_2018],
            L_shoot_V_2019=[L_shoot_V_2019],
            L_shoot_F_2018=[L_shoot_F_2018],
            L_shoot_F_2019=[L_shoot_F_2019],
            S_shoot_V_2018=[S_shoot_V_2018],
            S_shoot_V_2019=[S_shoot_V_2019],
            S_shoot_F_2018=[S_shoot_F_2018],
            S_shoot_F_2019=[S_shoot_F_2019],
            nb_V_2018=[nb_V_2018],
            nb_V_2019=[nb_V_2019],
            nb_F_2018=[nb_F_2018],
            nb_F_2019=[nb_F_2019],
            nb_L_2018=[nb_L_2018],
            nb_L_2019=[nb_L_2019],
            LA_2018=[LA_2018],
            LA_2019=[LA_2019],
            Elongation=[Elongation],
            Tapper=[Tapper],
            nb_Trunk_S=[nb_Trunk_S],
            ),
            columns='apple_tree trt NCI TSA_2018 TSA_2019 L_shoot_V_2018 L_shoot_V_2019 '
            'L_shoot_F_2018 L_shoot_F_2019 S_shoot_V_2018 S_shoot_V_2019 S_shoot_F_2018 S_shoot_F_2019 '
            'nb_V_2018 nb_V_2019 nb_F_2018 nb_F_2019 nb_L_2018 nb_L_2019 LA_2018 LA_2019 '
            'Elongation Tapper nb_Trunk_S'.split(' '))
        return df


class Branches2(Tree2):


    def dataframe(self):
        """
        - nb_F_2018/2019:
            Number of floral buds on the branch in 2018 and 2019
        - LA_2018/2019:
            Leaf area (cm2) in 2018 and 2019

        """
        #tree = super(Branches, self).dataframe()
        
        g = self.g
        dates = self.dates()
        lines = g.property('_line')

        self.upscale('length', operator=sum)
        self.upscale('leaf_area_2019', operator=sum)

        apple_tree = g.index(1)

        vs18 = [v for v, d in dates.items() if d !='2019']
        vs19 = [v for v, d in dates.items()
                    if (d =='2019') or
                       (isinstance(lines.get(v), list))
               ]
        setv18 = set(vs18) 
        setv19 = set(vs19)
        
        # Extract branches first vertex of each branch beared by the trunk
        trunk = g.Trunk(2)
        brs = [b for v in trunk for b in g.Sons(v, EdgeType='+')]

        labels18 = {}
        for label in list('ABCDE'):
            labels18[label] = [v for v in vs18 if g.class_name(v) ==label]

        labels19 = {}
        for label in list('ABCDE'):
            labels19[label] = [v for v in vs19 if g.class_name(v) ==label]

        #vtrunk = g.Trunk(3)
        #vbrs = [b for v in vtrunk for b in g.Sons(v, EdgeType='+')]

        #anchors = dict((b, b-1) for b in brs)
        #heights = self.heights()

        # Branch Heights
        #tip_height = float(heights[vtrunk[-1]])
        #b_dists = [(heights[anchors[v]]+1)/tip_height for v in brs]

        # TODO: Missing data for diameter
        #default_diameter = 0.3
        

        # Compute diameter
        #d18 = g.property('diameter_b2018').copy()
        #d18.update(g.property('diameter_2018'))
        #d19 = g.property('diameter_b').copy()
        #d19.update(g.property('diameter'))

        #BSA_2018 = [self.diameter(b, '2018', d18, d19) for b in brs]
        #BSA_2019 = [self.diameter(b, '2019', d18, d19) for b in brs]
        
        #B_pos = [(heights[anchors[b]]+1) for b in brs]
        #B_dist = b_dists
        
        # WIP
        #length = g.property('length')
        #B_lgth_2018 = [self.branch_length(b, is_19=False) for b in brs]
        #B_lgth_2019 = [self.branch_length(b, is_19=True) for b in brs]

        # TODO
        def isdyn(v):
            return isinstance(lines.get(v), list)

        def br18(b):
            brs = g.Descendants(b)
            branch = [a for a in brs if a in setv18]
            return branch

        def br19(b):
            brs = g.Descendants(b)
            branch = [a for a in brs if a in setv19]
            return branch


        # Number of floral buds in 2018 and 2019
        nb_F_2018 = [sum(1 for v in br18(b) if self.is_flo('2018', v))
                        for b in brs ]
        nb_F_2019 = [sum(1 for v in br19(b) if self.is_flo_19(v))
                        for b in brs ]

        #
        LA_2018 = [sum(self.la18(v) for v in br18(b)) for b in brs]
        LA_2019 = [sum(self.la19(v) for v in br19(b)) for b in brs]

        branch_name = [g.label(b) for b in brs]
        n= len(brs)

        columns= """
        apple_tree
        branch
        nb_F_2018
        nb_F_2019
        LA_2018
        LA_2019
        """.split()
        df = pd.DataFrame( OrderedDict(
            apple_tree=[apple_tree]*n,
            branch=branch_name,
            nb_F_2018=nb_F_2018,
            nb_F_2019=nb_F_2019,
            LA_2018=LA_2018,
            LA_2019=LA_2019,
            ),
            columns=columns)
        return df

def forest2(fns=None, save=False):
    if fns is None:
        fns=load_data()

    dfs = []
    errors=[]
    for fn in fns:
        try:
            print('#'*80)
            print("MTG is ", fn)
            g=MTG(fn, has_date=True)
            df = Tree2(g)()
            dfs.append(df)
        except:
            print('#'*80)
            print('Error with file {}'.format(fn))
            print('#'*80)
            print()
            errors.append(fn)
            continue

    all_df = pd.concat(dfs, axis=0, ignore_index=True)
    if save:
        all_df.to_csv('trees2.csv', sep=';', index=False)


    return all_df, errors

def branches2(fns=None, save=False):
    if fns is None:
        fns=load_data()

    dfs = []
    errors=[]
    for fn in fns:
        try:
            print('#'*80)
            print("MTG is ", fn)
            g=MTG(fn, has_date=True)
            df = Branches2(g)()
            dfs.append(df)
        except:
            print('#'*80)
            print('Error with file {}'.format(fn))
            print('#'*80)
            print()
            errors.append(fn)
            continue

    all_df = pd.concat(dfs, axis=0, ignore_index=True)
    if save:
        all_df.to_csv('branches2.csv', sep=';', index=False)


    return all_df, errors
