import bettifold as bf
import matplotlib.pyplot as plt
from barcodes import Barcode

output = bf.bettifold("UCCACAGGCAG",
                       # remember to specify which type of bonds are allowed
                       bond_constraints = [bf.watson_creek, bf.min_bond_len(4)],
                       folding_constraints = [bf.min_bonds(2), bf.max_bonds(5)],
                       folder = "rna-test",
                       maxdim = 1
                       )

rips = output["rips"]
rips.plot()
plt.show()
br = Barcode(rips)
br.plot_barcode()