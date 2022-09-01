import subprocess
from itertools import combinations
import os
import numpy as np
from tqdm import tqdm
from sklearn.manifold import MDS
from matplotlib import pyplot as plt
from ripser import Rips

ALLOWED_BONDS = {
    "A" : ["U"],
    "U" : ["A", "G"],
    "C" : ["G"],
    "G" : ["C","U"],
}


def distance_file(file1, file2):

    command = f"java -jar ASPRAlign.jar -a {file1} {file2} -d"

    out = subprocess.run(command.split(), capture_output=True)
    distance = float(out.stdout.decode().replace("Distance = ","").replace("\n",""))
    
    return distance


def min_bond_len(n):
    return lambda bond : abs(bond[0] - bond[1]) >= n

def max_bond_len(n):
    return lambda bond : abs(bond[0] - bond[1]) <= n


def is_valid_folding(rna, folding, constraints = [min_bond_len(4)]):
    consumed_nucleotides = set()
    for bond in folding:
        i,j = bond
        for cons in constraints:
            if not cons(bond):
                return False

        if i in consumed_nucleotides or j in consumed_nucleotides:
            return False

        n1,n2 = rna[i-1], rna[j-1]
        
        if n2 not in ALLOWED_BONDS[n1]:
            return False
        
        consumed_nucleotides.add(i)
        consumed_nucleotides.add(j)
    
    return True


def all_foldings(rna, bond_constraints = [min_bond_len(4)]):
    n = len(rna)
    possible_edges = list(combinations(range(1,n+1), 2))
    for k in range(1, n//2 + 1):
        foldings = combinations(possible_edges, k)
        for f in foldings:
            if is_valid_folding(rna,f, bond_constraints):
                yield f






rna = "ACGUUCAU"
bond_constraints = [min_bond_len(1)]
folder = "test"





if not os.path.exists(folder):
    os.mkdir(folder)



print("COMPUTING FOLDINGS")

for i, folding in tqdm(enumerate(all_foldings(rna, bond_constraints))):
    with open(f"{folder}/{i}", "w") as f:
        f.write(rna + "\n")
        bonds = ";".join([str(b).replace(" ","") for b in folding])
        f.write(bonds)


print("COMPUTING DISTANCES")


files = os.listdir(folder)
n_files = len(files)
distances = np.zeros((n_files, n_files))

pairs = set()
for f1 in files:
    for f2 in files:
        if f1 == f2: continue
        if (f1, f2) in pairs or (f2,f1) in pairs:
            continue

        pairs.add((f1,f2))


max_d = 0
for f1,f2 in tqdm(pairs):
    i,j = int(f1), int(f2)
    path_f1 = f"{folder}/{f1}"
    path_f2 = f"{folder}/{f2}"
    distances[i,j] = distance_file(path_f1,path_f2)
    distances[j,i] = distances[i,j]
    if distances[i,j] > max_d:
        max_d = distances[i,j]

# normalization
distances = distances / max_d

print("COMPUTING MDS")

mds = MDS(metric = True, dissimilarity="precomputed")
X = mds.fit_transform(distances)


print("COMPUTING RIPS")

rips = Rips(maxdim = 2)
#TODO : check if we can feed in the distances directly, without mds
rips.fit_transform(X)
rips.plot()
plt.show()