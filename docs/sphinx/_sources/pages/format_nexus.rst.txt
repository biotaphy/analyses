=====
Nexus
=====

Description
===========

Nexus [#f1]_ is a data format containing character state data.  The BiotaPhy
project often uses Nexus to store trees with annotations for the nodes in those
trees.

Example
=======

 ::

    #NEXUS
    
    BEGIN TAXA;
        DIMENSIONS NTAX=3;
        TAXLABELS
            A
            B
            C
      ;
    END;
    
    BEGIN TREES;
        TREE 1 = (A,(B,C));
    END;

.. rubric :: Footnotes

.. [#f1] http://wiki.christophchamp.com/index.php?title=NEXUS_file_format
