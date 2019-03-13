=====
NeXML
=====

Description
===========

NeXML [#f1]_ [#f2]_ is an XML-based data format used for representing phylogenetic
data.

Example
=======

 ::

    <?xml version="1.0" encoding="ISO-8859-1"?>
    <nex:nexml
        version="0.9"
        xsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"
        xmlns="http://www.nexml.org/2009"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:xml="http://www.w3.org/XML/1998/namespace"
        xmlns:nex="http://www.nexml.org/2009"
        xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
    >
        <otus id="d0">
            <otu id="d1" label="A" />
            <otu id="d2" label="B" />
            <otu id="d3" label="C" />
        </otus>
        <trees id="d4" otus="d0">
            <tree id="d5" xsi:type="nex:FloatTree">
                <node id="d6" />
                <node id="d7" otu="d1" />
                <node id="d8" />
                <node id="d9" otu="d2" />
                <node id="d10" otu="d3" />
                <rootedge id="d11" target="d6" />
                <edge id="d12" source="d6" target="d7" />
                <edge id="d13" source="d6" target="d8" />
                <edge id="d14" source="d8" target="d9" />
                <edge id="d15" source="d8" target="d10" />
            </tree>
        </trees>
    </nex:nexml>

.. rubric :: Footnotes

.. [#f1] http://www.nexml.org/
.. [#f2] Rutger A. Vos, James P. Balhoff, Jason A. Caravas, Mark T. Holder, 
           Hilmar Lapp, Wayne P. Maddison, Peter E. Midford, Anurag Priyam, 
           Jeet Sukumaran, Xuhua Xia, Arlin Stoltzfus; NeXML: Rich, Extensible, 
           and Verifiable Representation of Comparative Data and Metadata, 
           Systematic Biology, Volume 61, Issue 4, 1 July 2012, Pages 675–689, 
           https://doi.org/10.1093/sysbio/sys025