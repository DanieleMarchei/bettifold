import bettifold as bf
import matplotlib.pyplot as plt

output = bf.bettifold("ACUGCAGUAUC", 
                      bond_constraints = [bf.watson_creek_wobble, bf.min_bond_len(4)],
                      folding_constraints = [bf.min_bonds(3), bf.max_bonds(5)],
                      maxdim = 1,
                      folder = "test"
                    )
output["rips"].plot()
plt.show()