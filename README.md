# BettiFold

This is a code implemetation for the method proposed by the paper ``Topological Analysis of RNA Folding Space" by Daniele Marchei and Emanuela Merelli.

This code allows the user to input a single RNA primary structure and obtain in output a Vietoris-Rips filtration of its folding space.

Altough it was initialy created for RNA-related experiements, this method is completely general and can work with any type of sequences, even non-biological ones.

## Installation
Clone this repository into your working directory and install the following packages:
- numpy
- ripser
- sklearn
- tqdm
- pathos

Then download AspraLign from [here](https://github.com/bdslab/aspralign) and put the executable jar ``ASPRAlign.jar'' into the same folder of bettifold.py.


## Usage


```python
import bettifold as bf
import matplotlib.pyplot as plt

output = bf.bettifold("ACUGCAGUC")
                    
output["rips"].plot()
plt.show()
```

Here are the paramenters that can be accepted by ```bettifold```.
```python
bettifold(seq, 
          bond_constraints    = [watson_creek_wobble, min_bond_len(4)],
          folding_constraints = [],
          distance_func       = aspra,
          folder              = "output",
          n_foldings          = "all",
          use_mds             = False,
          maxdim              = 2,
          n_processes         = None
          )
```
Parametrs:
- ```seq``` *( str )* : A string representing the sequence of letters you want to find the foldings of.
- ```bond_constraints``` *( [func] )* : A list of functions ```f : Bond -> bool``` representing constrains each bond must satisfy.
- ```folding_constraints``` *( [func] )* : A list of functions ```f : tuple[Bond] -> bool``` representing constrains each folding must satisfy.
- ```distance_func``` ( *func* ) : A function between two files describing a folding. The default is the Aspra Distance<sup>1</sup>.
- ```folder``` ( str ) : The name of the folder where you want to save the generated foldings. If a folder with this name is already present it will be deleted.
- ```n_foldings``` ( *"all" | int* ) : The number of folding to be generated. If "all", it will generate all foldings that satisfy the constraints.
- ```use_mds``` ( *bool* ) : If True, MultiDimensional Scaling will be performed on the distance matrix before the Rips filtration. If False, the distance matrix will be fed directly into ripser.
- ```maxdim``` ( *int* ) : The number of homology classes to be calculated. We calculate them using the Ripser package<sup>2</sup>.
- ```n_processes``` ( *None | int* ) : The number of cores to be used. If None, the max number of available cores will be used.


<sup>1</sup>[https://github.com/bdslab/aspralign](https://github.com/bdslab/aspralign)

<sup>2</sup>[https://github.com/scikit-tda/ripser.py](https://github.com/scikit-tda/ripser.py)


Output:

```bettifold``` returns a dict with the following fields:
- ```distances``` : A numpy distance matrix between all enumerated foldings.
- ```mds_emb``` : A numpy array of points generated by using MultiDimensional Scaling. Each point represents a folding. This will be None if ```use_mds``` is False.
- ```rips``` : A ripser object representing the homology classes of the space of points.

## Bond

```Bond``` is a dataclass with the following attributes:
```python
class Bond:
    idxs : tuple[int, int]
    nucl : tuple[str, str]
```
in which ```idxs``` represents the tuple of 1-based indicies where the pair of nucleotides (indicated by ```nucl```) are located in the sequence.

## Constraints

We provide some built-in constraint functions. Each function returns a bool that says if the input satisfies the contraint or not.
### Bond Constraints
- ```min_bond_len``` ( *int* ) : Takes an integer n and returns a function that checks if the input bond has at least length n.
- ```max_bond_len``` ( *int* ) : Takes an integer n and returns a function that checks if the input bond has at most length n.
- ```watson_creek_wobble``` ( *bond* ) : Returns True if the bond is a Watson-Creek-Wobble bond. Returns False otherwise.
- ```watson_creek``` ( *bond* ) : Returns True if the bond is a Watson-Creek bond. Returns False otherwise.

### Folding Constraints
- ```min_bonds``` ( *int* ) : Takes an integer n and returns a function that checks if the input folding has at least n bonds.
- ```max_bonds``` ( *int* ) : Takes an integer n and returns a function that checks if the input folding has at most n bonds.


## Example of a non-RNA sequence
```python
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
```