import bettifold as bf
import matplotlib.pyplot as plt


def my_distance(file1, file2):
    pass

#            remember : bond.idxs are 1-based
def protein_bond(bond : bf.Bond):
    pass

def my_constraint(folding : tuple[bf.Bond]):
    pass


output = bf.bettifold("MLVPNGEAADYRARKGVLVV",
                       bond_constraints = [protein_bond, bf.min_bond_len(5)],
                       folding_constraints = [my_constraint, bf.max_bonds(10)],
                       distance_func = my_distance,
                       n_foldings = 1000,
                       folder = "protein-test",
                       maxdim = 3
                       )

output["rips"].plot()
plt.show()