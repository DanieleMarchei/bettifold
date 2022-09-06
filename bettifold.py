from dataclasses import dataclass
from pathos.multiprocessing import ProcessingPool as Pool
import os
import shutil
import numpy as np
from tqdm import tqdm
from sklearn.manifold import MDS
from ripser import Rips
from functools import partial
import subprocess
from itertools import combinations


@dataclass
class Bond:
    idxs : tuple[int, int]
    nucl : tuple[str, str]

    def __repr__(self):
        return f"({self.idxs[0]},{self.idxs[1]})"
    
    def __hash__(self) -> int:
        return hash(self.idxs)



WATSON_CREEK_WOBBLE = {
    "A" : ["U"],
    "U" : ["A", "G"],
    "C" : ["G"],
    "G" : ["C","U"],
}

WATSON_CREEK = {
    "A" : ["U"],
    "U" : ["A"],
    "C" : ["G"],
    "G" : ["C"],
}

def aspra(file1, file2):

    command = f"java -jar ASPRAlign.jar -a {file1} {file2} -d"

    out = subprocess.run(command.split(), capture_output=True)
    distance = float(out.stdout.decode().replace("Distance = ","").replace("\n",""))
    
    return distance

# ------------ BOND CONSTRAINTS ------------

def min_bond_len(n):
    return lambda bond : abs(bond.idxs[0] - bond.idxs[1]) >= n

def max_bond_len(n):
    return lambda bond : abs(bond.idxs[0] - bond.idxs[1]) <= n

def watson_creek_wobble(bond):
    return bond.nucl[1] in WATSON_CREEK_WOBBLE[bond.nucl[0]]

def watson_creek(bond):
    return bond.nucl[1] in WATSON_CREEK[bond.nucl[0]]


# ------------ FOLDING CONSTRAINTS ------------

def min_bonds(n):
    return lambda folding : len(folding) >= n

def max_bonds(n):
    return lambda folding : len(folding) <= n

# ------------

def _is_valid_folding(bond_constraints, folding_constraints, folding):
    for cons in folding_constraints:
        if not cons(folding):
            return False

    consumed_nucleotides = set()
    for bond in folding:
        for cons in bond_constraints:
            if not cons(bond):
                return False

        i,j = bond.idxs
        if i in consumed_nucleotides or j in consumed_nucleotides:
            return False

        consumed_nucleotides.add(i)
        consumed_nucleotides.add(j)
    
    return True

def _all_foldings_k(seq, bond_constraints, folding_constraints, possible_edges, k):
    foldings = combinations(possible_edges, k)
    vaild_foldings = set()
    for f in foldings:
        _f = tuple([Bond((i,j), (seq[i-1], seq[j-1])) for i,j in f])
        if _is_valid_folding(bond_constraints, folding_constraints, _f):
            vaild_foldings.add(_f)
            
    return vaild_foldings

def _all_foldings(seq, bond_constraints, folding_constraints, n_processes = None):
    n = len(seq)
    possible_edges = list(combinations(range(1,n+1), 2))
    _func = partial(_all_foldings_k, seq, bond_constraints, folding_constraints, possible_edges)
    with Pool(n_processes) as pool:
        results = pool.imap(_func, range(1, n//2 + 1))
        for valid_foldings in results:
            for f in valid_foldings:
                yield f


def _sample_foldings(n_samples):
    raise NotImplementedError("Random sampling of foldings has not been implemented yet.")
    # def _sample(rna, bond_constraints, folding_constraints): 
    #   pass
    # return _sample


def _dist(dist_func, folder, files):
    f1,f2 = files
    i,j = int(f1), int(f2)
    path_f1 = f"{folder}/{f1}"
    path_f2 = f"{folder}/{f2}"

    return i,j,dist_func(path_f1,path_f2)


def bettifold(seq, 
              bond_constraints = [watson_creek_wobble, min_bond_len(4)],
              folding_constraints = [],
              distance_func = aspra,
              folder = "output",
              n_foldings = "all",
              use_mds = False,
              maxdim = 2,
              n_processes = None):

    if n_foldings == "all":
        sample_foldings = _all_foldings
    elif type(n_foldings) == int:
        sample_foldings = _sample_foldings(n_foldings)
    else:
        raise TypeError("The argument 'sampling' can only take values in (\"all\" | int).")


    if os.path.exists(folder):
        shutil.rmtree(folder)
    
    os.mkdir(folder)



    print("COMPUTING FOLDINGS")

    for i, folding in tqdm(enumerate(sample_foldings(seq, bond_constraints, folding_constraints, n_processes))):
        with open(f"{folder}/{i}", "w") as f:
            f.write(seq + "\n")
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
    

    with Pool(n_processes) as p:
        _d = partial(_dist, distance_func, folder)
        results = p.imap(_d, pairs)

        for i,j,d in tqdm(results, total=len(pairs)):
            distances[i,j] = d
            distances[j,i] = d


    rips = Rips(maxdim = maxdim)
    X = None
    if use_mds:
        print("COMPUTING MDS")
        mds = MDS(metric = True, dissimilarity="precomputed")
        X = mds.fit_transform(distances)

        print("COMPUTING RIPS")
        rips.fit_transform(X, distance_matrix = False)

    else:
        print("COMPUTING RIPS")
        rips.fit_transform(distances, distance_matrix = True)


    output = {
        "distances" : distances,
        "mds_emb" : X,
        "rips" : rips
    }

    return output