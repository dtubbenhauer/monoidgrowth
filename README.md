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

# Setting up the monoid algebras and computing the growth rate

Using the following code, we obtain the constant \(v_0w_0^T[1]\) (in the notation of <a href="https://arxiv.org/abs/2307.03044">https://arxiv.org/abs/2307.03044</a>) in the periodic expression of 
\(a(n)\) for each projective indecomposable \(kM\)-module \(V\), where \(M\) is a monoid of order 4. We do this by first setting up the monoid algebra as a matrix algebra 
in Magma. Given the multiplication tables of the monoids, we create the matrices which generate this algebra, using the functions *left_regular_representation* 
and *matrix_to_magma*. Then, using the function *generate_magma_code* we set up the corresponding monoid algebra of a multiplication table and computes the constant \(v_0w_0^T[1]\) for each projective indecomposable \(kM\)-module. 
Here is Python code:

```
import numpy as np
import csv
import os
import ast
from scipy.sparse import csr_matrix
import subprocess

def left_regular_representation(table, generators=None):
    """
    Generate left regular representation matrices from the multiplication table.
    Uses the rows of the table to construct matrices.

    Parameters:
        table (numpy.ndarray): The multiplication table of the monoid.
        generators (list[int], optional): Indices of the generators. Defaults to all rows.

    Returns:
        list[numpy.ndarray]: List of matrices representing the left regular representation.
    """
    size = table.shape[0]
    if generators is None:
        generators = list(range(size))  # Use all rows by default
    matrices = []
    for row in generators:
        col_indices = table[row, :] - 1  # Subtract 1 for 0-based indexing
        data = [1] * size  # Non-zero entries
        matrix = csr_matrix((data, (range(size), col_indices)), shape=(size, size))
        matrices.append(matrix.todense())
    return matrices

def matrix_to_magma(matrix, index, prefix="M"):
    """
    Convert a NumPy array to a Magma-compatible string representation.

    Parameters:
        matrix (numpy.ndarray): Matrix to convert.
        index (int): Index of the matrix.
        prefix (str): Prefix for the Magma variable name.

    Returns:
        str: Magma-compatible matrix string.
    """
    mat_list = matrix.tolist()
    rows_str = ["[" + ", ".join(map(str, row)) + "]" for row in mat_list]
    return f"{prefix}_{index}:=A![{',\n '.join(rows_str)}];"


def generate_magma_code(table, prime=0, left_generators=None):
    """
    Generate Magma code for left regular representations of the monoid defined by 'table'.

    Parameters:
        table (numpy.ndarray): Multiplication table of the monoid.
        prime (int): Prime field order (0 for Rationals).
        left_generators (list[int], optional): Indices of left generators.

    Returns:
        str: Magma code for the monoid.
    """
    n = table.shape[0]
    left_reps = left_regular_representation(table, left_generators)
    
    field = f"GF({prime})" if prime > 0 else "Rationals()"
    # We set up the (left) regular representation in Magma using the matrices  
    magma_code = f"A:=MatrixAlgebra({field},{n});\n"
    magma_code += "\n".join(matrix_to_magma(matrix, i + 1, prefix="M") for i, matrix in enumerate(left_reps))
    all_left_gen_labels = [f"M_{i + 1}" for i in range(len(left_reps))]
    magma_code += f"\nM:=sub<A|{','.join(all_left_gen_labels)}>;\n"
    magma_code += f"V:= RModule({field}, {n});\n"
    magma_code += "m:= map< CartesianProduct(M,V) -> V | t :-> t[2]*t[1]>;\n"
    magma_code += "W:=Module(M, m);\n"
    magma_code += "J:=IndecomposableSummands(W);\n"

    # We now compute the tensor products
    magma_code += f"""
for l in [1..#J] do
    W:=J[l];
    V:=RModule({field},1);
    m:=map<CartesianProduct(M,V) -> V | t :-> t[2]>;
    Triv:=Module(M,m);

    X:=[TensorProduct(Triv,Triv)];
    nn:=0;
    max:=40;
    while nn lt #X and nn le max do
        nn+:=1;
        M2:=TensorProduct(X[nn],W);
        XM2:=IndecomposableSummands(M2);
        for j in [1..#XM2] do
            new:=1;
            for i in [1..#X] do
                if IsIsomorphic(XM2[j],X[i]) then new:=0; end if;
            end for;
            if new eq 1 then X:=Append(X,XM2[j]); end if;
        end for;
    end while;

    fus:=[[0 : i in [1..#X]] : j in [1..#X]];
    for i in [1..#X] do
        for j in [1..#X] do
            Z:=IndecomposableSummands(TensorProduct(X[i],W));
            z:=0;
            for k in [1..#Z] do
                if IsIsomorphic(X[j],Z[k]) then z:=z+1; end if;
            end for;
            fus[j][i]:=z;
        end for;
    end for;

    Ma:=Matrix(Rationals(),fus);
    w0:=Eigenspace(Ma,Dimension(W)).1;
    vv0:=Eigenspace(Transpose(Ma),Dimension(W)).1;
    v0:=vv0/ScalarProduct(w0,vv0);
    &+[(Matrix(Rationals(),#X,1,ElementToSequence(v0))*Matrix(Rationals(),1,#X,
ElementToSequence(w0)))[i][1] : i in [1..#X]];
end for;
"""
    return magma_code
```

See also the appendix of <a href="https://arxiv.org/abs/2307.03044">https://arxiv.org/abs/2307.03044</a> to understand the part of the Magma code that computes the tensor product decompositions. Next, we use Python to run the Magma code for each monoid in the second data file above.

```
def run_magma_with_load(magma_code, temp_filename="temp_script.magma", timeout=60):
    """
    Write Magma code to a file and use `load` in Magma to execute it.
    """
    magma_path = "/Applications/Magma/magma"
    with open(temp_filename, "w", encoding="utf-8") as f:
        f.write(magma_code)
    
    load_command = f'load "{temp_filename}";'

    try:
        process = subprocess.Popen(
            [magma_path, '-b', '-e', load_command],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate(timeout=timeout)
        if stderr:
            return f"Error: {stderr.strip()}"
        return stdout.strip()
    except subprocess.TimeoutExpired:
        return "Error: Magma process timed out."
    except FileNotFoundError:
        return "Magma executable not found. Please check the path."
    except Exception as e:
        return f"Error running the Magma script: {e}"
    finally:
        if os.path.exists(temp_filename):
            os.remove(temp_filename)
            
def process_csv(input_csv, output_csv, start_row=0, prime=0, timeout=60):
    with open(input_csv, "r") as infile, open(output_csv, "w", newline="") as outfile:
        reader = list(csv.DictReader(infile))
        fieldnames = ["Monoid ID", "Order of Group of Units", "Magma Output"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader[start_row:]:
            monoid_id = row["Semigroup ID"]
            try:
                table = np.array(ast.literal_eval(row["Multiplication Table"]))
            except Exception:
                continue  
            order_of_units = row["Order of Group of Units"]
            magma_code = generate_magma_code(table, prime=prime, left_generators=None)
            magma_output = run_magma_with_load(magma_code, timeout=timeout)
            writer.writerow({
                "Monoid ID": monoid_id,
                "Order of Group of Units": order_of_units,
                "Magma Output": magma_output
            })
            outfile.flush()

process_csv("semigroups_data.csv", "semigroups_results.csv", start_row=0,prime=11, timeout=100)   
```

The output is written to a file as in the first link above. The file has three columns: the monoid ID, the order of the group of units, and the magma output giving the constants \(v_0w_0^T[1]\) for each projective indecomposable \(kM\)-module.

We note that the *smallsemi* package of GAP gives the multiplication table of all monoids up to isomorphism and anti-isomorphism, so to check for all monoids of order 4 we also have to consider the opposites of 
the monoid algebras considered above. This may be done by working with the right regular representations instead.

# Erratum

Empty so far.
