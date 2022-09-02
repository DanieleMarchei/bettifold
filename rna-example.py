import bettifold as bf
import matplotlib.pyplot as plt


output = bf.bettifold("UCCACAGGCAG",
                       # remember to specify which type of bonds are allowed
                       bond_constraints = [bf.watson_creek, bf.min_bond_len(4)],
                       folding_constraints = [bf.min_bonds(2), bf.max_bonds(5)],
                       folder = "rna-test",
                       maxdim = 2
                       )

output["rips"].plot()
plt.show()