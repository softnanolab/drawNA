DrawNA's Documentation
======================
drawNA is a Python Library which can currently be used for generating a system comprised of ssDNA/dsDNA strands. Using a set of vertices, lines (edges) are generated which act as centrelines for single/double strands of DNA. It is modelled using the oxDNA model and can be exported in this format.

The end goal of the package is to autonomously generate many possible DNA origami scaffold configurations and subsequently their complimentary staple strands for a given 2D shape.

About
=====
Polygon
^^^^^^^
* Generate a polygon, given its verticies, where every edge is individually accessible
* Write the polygon to .STL and .PLY file formats

oxDNA
^^^^^
* Create an oxDNA system with tools to control the length, sequence, base-pairs/nucleotides per 2π turn and more.
* Export to oxDNA file format or LAMMPS data format
* Read oxDNA and LAMMPS data and dump formats

Lattice
^^^^^^^
* Create a lattice within a polygon
* Allocate a route for a strand of DNA to take to fill the polygon
* Create an oxDNA system of the scaffold strand which follows the route


Install
=======
To install the code using PyPi::

    git clone https://github.com/softnanolab/drawNA
    cd drawNA
    pip install .

FAQ
===

.. toctree::
   :maxdepth: 2

Examples
========


`drawNA`
========
.. automodule:: drawNA
   :members:

`drawNA.polygons`
-----------------
.. automodule:: drawNA.polygons
   :members:

`drawNA.readers`
-----------------
.. automodule:: drawNA.readers
   :members:

`drawNA.tools`
-----------------
.. automodule:: drawNA.tools
   :members:

`drawNA.oxdna`
==============
.. automodule:: drawNA.oxdna
   :members:

`drawNA.oxdna.nucleotide`
-------------------------
.. automodule:: drawNA.oxdna.nucleotide
   :members:

`drawNA.oxdna.strand`
---------------------
.. automodule:: drawNA.oxdna.strand
   :members:

`drawNA.oxdna.system`
---------------------
.. automodule:: drawNA.oxdna.system
   :members:

`drawNA.oxdna.utils`
--------------------
.. automodule:: drawNA.oxdna.utils
   :members:

`drawNA.lattice`
================
.. automodule:: drawNA.lattice
   :members:

`drawNA.lattice._lattice`
-------------------------
.. automodule:: drawNA.lattice._lattice
   :members: 

`drawNA.lattice.edge`
---------------------
.. automodule:: drawNA.lattice.edge
   :members:

`drawNA.lattice.node`
---------------------
.. automodule:: drawNA.lattice.node
   :members:

`drawNA.lattice.route`
----------------------
.. automodule:: drawNA.lattice.route
   :members:

`drawNA.lattice.utils`
----------------------
.. automodule:: drawNA.lattice.utils
   :members:

* :ref:`genindex`
* :ref:`modindex`
