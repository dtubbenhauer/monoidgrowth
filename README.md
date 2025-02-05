# Code and Erratum for *Growth Problems for Representations of Finite Monoids*

I collected a bit of <a href="https://www.gap-system.org/">GAP</a>, Python and 
<a href="https://magma.maths.usyd.edu.au/magma/">Magma</a> code relevant for the paper *Growth Problems for Representations of Finite Monoids*
<a href="https://arxiv.org/abs/2412.01283">https://arxiv.org/abs/2412.01283</a> on this page.

Here is are two links where to download the data associated to the project (all in .csv files): 
- <a href="https://github.com/dtubbenhauer/monoidgrowth/blob/main/semigroups_results.csv">Click</a>.
- <a href="https://github.com/dtubbenhauer/monoidgrowth/blob/main/semigroups_data.csv">Click</a>.

An Erratum for the paper *Growth Problems for Representations of Finite Monoids* can be found at the bottom of the page.

# Contact

If you find any errors in the paper *Growth Problems for Representations of Finite Monoids* **please email me**:

[dtubbenhauer@gmail.com](mailto:dtubbenhauer@gmail.com?subject=[GitHub]%web-reps)

Same goes for any errors related to this page.

# Monoids of small order

The *Smallsemi* package in GAP contains data for all semigroups of order up to 8. The following GAP code creates a 
.csv file (first link above) with three columns: the semigroup ID, the order of the group of units, and the 
multiplication table. This contains data for all monoids of order 4 up to isomorphism and anti-isomorphism.

```
LoadPackage("Smallsemi");
LoadPackage("Semigroups");

outputFile := OutputTextFile("semigroups_data.csv", false);
AppendTo(outputFile, "Semigroup ID,Multiplication Table,Order of Group of Units\n");

semigroups := AllSmallSemigroups(4, IsMonoidAsSemigroup, true);
for sg in semigroups do
    id := IdSmallSemigroup(sg);
    multTable := RecoverMultiplicationTable(Size(sg), id[2]);
    multTableFormatted := Concatenation("[", JoinStringsWithSeparator(List(multTable, row -> Concatenation
    ("[", JoinStringsWithSeparator(row, ", "), "]")), ", "), "]");
    orderOfUnits := Size(GroupOfUnits(Semigroup(sg)));
    csvLine := Concatenation(
        "\"", String(id), "\"", ",",
        "\"", multTableFormatted, "\"", ",",
        String(orderOfUnits), "\n"
    );
    AppendTo(outputFile, csvLine);
od;

CloseStream(outputFile);
```

# The full transformation monoid

The following function computes the multiplication table for the full transformation monoid in Python.

```
import itertools
import numpy as np

def full_transformation_monoid_matrix(n):
    """
    Generate the multiplication table for the full transformation monoid T_n.

    Parameters:
        n (int): The size of the set X = {1, 2, ..., n}, representing the domain 
                 and codomain of the transformations.

    Returns:
        tuple:
            - numpy.ndarray: The multiplication table for the full transformation monoid T_n, 
                             as a 2D array where each entry contains the index of the resulting function.
            - list[tuple[int]]: The list of all functions from X to X, represented as tuples.
                                Each tuple corresponds to a function in the monoid.
    """
    X = list(range(1, n + 1))
    functions = list(itertools.product(X, repeat=n))
    function_to_index = {f: i + 1 for i, f in enumerate(functions)}
    table_size = len(functions)
    matrix = [[0] * table_size for _ in range(table_size)]

    for i, f in enumerate(functions):
        for j, g in enumerate(functions):
            h = tuple(f[g[x - 1] - 1] for x in X)
            matrix[i][j] = function_to_index[h]

    return np.array(matrix), functions
```

# Direct product of monoids

If we have the multiplication tables for two monoids, \(M_1, M_2\), we can use the following Python code to obtain the multiplication table for \(M_1\times M_2\).

```
def direct_product_table(M1_table, M2_table):
    """
    Compute the direct product of two monoids given their multiplication tables.

    Parameters:
        M1_table (list[list[int]]): Multiplication table for the first monoid.
        M2_table (list[list[int]]): Multiplication table for the second monoid.

    Returns:
        numpy.ndarray: The multiplication table for the direct product of the two monoids.
    """
    m = len(M1_table)
    n = len(M2_table)
    size = m * n
    result = [[0] * size for _ in range(size)]

    def pair_to_index(a, b):
        return (a - 1) * n + (b - 1)

    for a in range(1, m + 1):
        for b in range(1, n + 1):
            left_idx = pair_to_index(a, b)

            for a_prime in range(1, m + 1):
                for b_prime in range(1, n + 1):
                    right_idx = pair_to_index(a_prime, b_prime)
                    x = M1_table[a - 1][a_prime - 1]
                    y = M2_table[b - 1][b_prime - 1]
                    result_idx = pair_to_index(x, y)
                    result[left_idx][right_idx] = result_idx + 1

    return np.array(result)    
```

# Erratum

Empty so far.
