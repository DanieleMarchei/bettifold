import bettifold as bf
import matplotlib.pyplot as plt
from ripser import Rips
import numpy as np
import subprocess

def my_pairwise_distances():
    command = f"java -jar all_pair_aspra_distances.jar my_foldings/"

    out = subprocess.run(command.split(), capture_output=True)
    lines = out.stdout.decode().split("\n")
    distances = np.zeros((50,50))

    del lines[-1]
    
    for line in lines:
        i,j,d = line.split(" ")
        i = int(i)
        j = int(j)
        d = float(d)

        distances[i,j] = d
        distances[j,i] = d
    
    return distances

bf.save_foldings("UCCACAGGCAG",
                       n_foldings = "all",
                       bond_constraints = [bf.watson_creek, bf.min_bond_len(4)],
                       folding_constraints = [bf.min_bonds(2), bf.max_bonds(5)],
                       folder = "my_foldings"
                       )


X = my_pairwise_distances()

rips = Rips(maxdim=2)
rips.fit_transform(X, distance_matrix=True)

rips.plot()
plt.show()