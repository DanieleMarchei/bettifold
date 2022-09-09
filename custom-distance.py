import bettifold as bf
import matplotlib.pyplot as plt
from ripser import Rips
import numpy as np

def my_pairwise_distances(foldings):
    for f in foldings:
        pass

    return np.zeros((10,10))

foldings = bf.enumerate_foldings("UCCACAGGCAG",
                       n_foldings = "all",
                       bond_constraints = [bf.watson_creek, bf.min_bond_len(4)],
                       folding_constraints = [bf.min_bonds(2), bf.max_bonds(5)],
                       folder = None
                       )

X = my_pairwise_distances(foldings)

rips = Rips(maxdim=2)
rips.fit_transform(X, distance_matrix=True)

rips.plot()
plt.show()