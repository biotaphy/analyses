======
Phylip
======

Description
===========

The Phylip format [#f1]_ can be used to store character data for a set of taxa.  The
"strict" Phylip format requires that the first line indicate the number of taxa
in the file and the number of characters each taxon will have.  The following
lines should include a label for each taxon (up to 10 characters) and the
character data for that taxon, which is split up into 10 character chunks that may be
padded with spaces (as below).

Example
=======

 ::
 
       3        15
    Rabbit    GCCCAACGTN NNATC
    Squirrel  CGAGGCTANA TATGG
    Lion      AAGCTNGTCC ATTGC

.. rubric :: Footnotes

.. [#f1] http://evolution.genetics.washington.edu/phylip/doc/main.html
