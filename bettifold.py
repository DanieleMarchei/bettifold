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
    import subprocess

    command = f"java -jar ASPRAlign.jar -a {file1} {file2} -d"

    out = subprocess.run(command.split(), capture_output=True)
    distance = float(out.stdout.decode().replace("Distance = ","").replace("\n",""))
    
    return distance

# ------------ BOND CONSTRAINTS ------------

def min_bond_len(n):
    return lambda rna_bond : abs(rna_bond[0][0] - rna_bond[0][1]) >= n

def max_bond_len(n):
    return lambda rna_bond : abs(rna_bond[0][0] - rna_bond[0][1]) <= n

def watson_creek_wobble(rna_bond):
    return rna_bond[1][1] in WATSON_CREEK_WOBBLE[rna_bond[1][0]]

def watson_creek(rna_bond):
    return rna_bond[1][1] in WATSON_CREEK[rna_bond[1][0]]


# ------------ FOLDING CONSTRAINTS ------------

def min_bonds(n):
    return lambda rna, folding : len(folding) >= n

def max_bonds(n):
    return lambda rna, folding : len(folding) <= n

# ------------

def _is_valid_folding(rna, folding, bond_constraints, folding_constraints):

    for cons in folding_constraints:
        if not cons(rna, folding):
            return False

    consumed_nucleotides = set()
    for bond in folding:
        i,j = bond
        rna_bond = (bond, (rna[i-1],rna[j-1]))
        for cons in bond_constraints:
            if not cons(rna_bond):
                return False

        if i in consumed_nucleotides or j in consumed_nucleotides:
            return False

        consumed_nucleotides.add(i)
        consumed_nucleotides.add(j)
    
    return True


def _all_foldings(rna, bond_constraints, folding_constraints):
    from itertools import combinations

    n = len(rna)
    possible_edges = list(combinations(range(1,n+1), 2))
    for k in range(1, n//2 + 1):
        foldings = combinations(possible_edges, k)
        for f in foldings:
            if _is_valid_folding(rna,f, bond_constraints, folding_constraints):
                yield f

def _sample_foldings(n_samples):
    raise NotImplementedError("Random sampling of foldings has not been implemented yet.")
    # def _sample(rna, bond_constraints, folding_constraints): 
    #   pass
    # return _sample



def bettifold(rna, 
              bond_constraints = [watson_creek_wobble, min_bond_len(4)],
              folding_constraints = [],
              distance_func = aspra,
              folder = "output",
              n_samples = "all", #all, number
              maxdim = 2):

    import os
    import numpy as np
    from tqdm import tqdm
    from sklearn.manifold import MDS
    from ripser import Rips

    
    if n_samples == "all":
        sample_foldings = _all_foldings
    elif type(n_samples) == int:
        sample_foldings = _sample_foldings(n_samples)
    else:
        raise TypeError("The argument 'sampling' can only take values in (\"all\" | int).")


    if not os.path.exists(folder):
        os.mkdir(folder)



    print("COMPUTING FOLDINGS")

    for i, folding in tqdm(enumerate(sample_foldings(rna, bond_constraints, folding_constraints))):
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
        distances[i,j] = distance_func(path_f1,path_f2)
        distances[j,i] = distances[i,j]
        if distances[i,j] > max_d:
            max_d = distances[i,j]

    # normalization
    # distances = distances / max_d

    print("COMPUTING MDS")

    mds = MDS(metric = True, dissimilarity="precomputed")
    X = mds.fit_transform(distances)


    print("COMPUTING RIPS")

    rips = Rips(maxdim = maxdim)
    #TODO : check if we can feed in the distances directly, without mds
    rips.fit_transform(X)

    output = {
        "rips" : rips,
        "points" : X,
        "distances" : distances
    }

    return output