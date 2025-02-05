# Code and Erratum for *Growth Problems for Representations of Finite Monoids*

I collected a bit of <a href="https://www.gap-system.org/">GAP</a> and 
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

# Erratum

Empty so far.
