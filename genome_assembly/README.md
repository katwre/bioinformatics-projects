Genome assembly
==================

<h2> De Bruijn Graph implementation with Eulerian walk-finder </h2>

Modern short-read assembly algorithms construct a de Bruijn graph by representing all k-mer prefixes and suffixes as nodes and then drawing edges that represent k-mers having a particular prefix and suffix [1].
Eulerian walk allows to reconstruct the DNA sequence from its fragments (k-mers) [2].


<hr>

[1] Phillip E C Compeau, Pavel A Pevzner and Glenn Tesler (2011). How to apply de Bruijn graphs to genome assembly. Nature Biotechnology 29, 987–991

[2] Pavel A. Pevzner, Haixu Tang and Michael S. Waterman (2001). An Eulerian path approach to DNA fragment assembly. Proc Natl Acad Sci U S A., 98(17): 9748–9753
