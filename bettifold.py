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
from random import randint, choice


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

def save_foldings(seq, n_foldings, stop_after_n_misses = 100000, 
            bond_constraints = [watson_creek_wobble, min_bond_len(4)],
            folding_constraints = [],
            folder = None,
            n_processes = None):
    

    if n_foldings == "all":
        sample_foldings = _all_foldings
    elif type(n_foldings) == int:
        sample_foldings = _sample_foldings(n_foldings, stop_after_n_misses)
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
    

def get_foldings(seq, n_foldings, 
            bond_constraints = [watson_creek_wobble, min_bond_len(4)],
            folding_constraints = [],
            n_processes = None):

    if n_foldings == "all":
        sample_foldings = _all_foldings
    elif type(n_foldings) == int:
        sample_foldings = _sample_foldings(n_foldings)
    else:
        raise TypeError("The argument 'sampling' can only take values in (\"all\" | int).")



    for i, folding in tqdm(enumerate(sample_foldings(seq, bond_constraints, folding_constraints, n_processes))):
        yield folding

class _ListDict(object):
    #FROM : https://stackoverflow.com/questions/15993447/python-data-structure-for-efficient-add-remove-and-random-choice

    def __init__(self, items):
        self.item_to_position = {}
        self.items = []
        # for extra efficiency, we are assuming no two items are equal
        for item in items:
            self.items.append(item)
            self.item_to_position[item] = len(self.items)-1

    def remove_item(self, item):
        position = self.item_to_position.pop(item)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position

    def choose_random_item(self):
        return choice(self.items)
    
    def pop_random(self):
        item = self.choose_random_item()
        self.remove_item(item)
        return item

def _random_folding(seq):
    n = len(seq)
    n_bonds = randint(1, n // 2)
    available_nucl = _ListDict(range(1,n+1))
    f = []
    for _ in range(n_bonds):
        n1 = available_nucl.pop_random()
        n2 = available_nucl.pop_random()
        if n1 > n2:
            n1,n2 = n2,n1
        f.append(Bond((n1,n2),(seq[n1-1],seq[n2-1])))
    
    f.sort(key = lambda bond : bond.idxs[0])
    
    return tuple(f)


def _sample_foldings(n_samples, stop_after_n_misses = 100000):
    def _sample(seq, bond_constraints, folding_constraints, n_processes): 
      foldings_found = set()
      _is_valid = partial(_is_valid_folding, bond_constraints, folding_constraints)
      n_misses = 0
      while len(foldings_found) < n_samples and n_misses < stop_after_n_misses:
        f = _random_folding(seq)
        if _is_valid(f) and f not in foldings_found:
            foldings_found.add(f)
            n_misses = 0
            yield f
        else:
            n_misses += 1

    return _sample


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
              stop_after_n_misses = 100000,
              use_mds = False,
              maxdim = 2,
              n_processes = None):

    save_foldings(seq, n_foldings, 
            bond_constraints=bond_constraints, 
            folding_constraints=folding_constraints, 
            folder = folder,
            n_processes=n_processes)


    print("COMPUTING DISTANCES")


    files = os.listdir(folder)
    n_files = len(files)
    distances = np.zeros((n_files, n_files))

    pairs = set()
    for f1 in range(n_files):
        for f2 in range(f1 + 1, n_files):
            pairs.add((str(f1),str(f2)))
    

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