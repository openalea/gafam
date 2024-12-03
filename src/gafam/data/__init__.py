import os, glob
from pathlib import Path

data_dir = Path(os.path.dirname(__file__)).resolve()

def files(dir=None):
    if dir is None:
        dir = data_dir
    else:
        dir=Path(dir)
    fns = dir.glob('*.txt')
    fns = sorted(fns, key= lambda x:int(x.name.split('.')[0][1:]))
    return fns

def file(fn):
    fns = list(data_dir.glob('%s.txt'%fn))
    return fns[0]

def treatment():
    fns = list(data_dir.glob('treatment.csv'))
    return fns[0]
