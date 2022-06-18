from multiprocessing import Pool
import os
import sys

def f(x):
    print(f'This is a test: {x}', file=sys.stdout)

with Pool(5) as p:
    p.map(f, range(10))
# f(1)
