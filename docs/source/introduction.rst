************
Introduction
************

Here we present an open-source Python library which focuses on
fitting the neutron resonance signal for neutron imaging
measurements. In this package, by defining the sample information
such as elements and thickness in the neutron path, one can extract
elemental/isotopic information of the sample. Various sample
types such as layers of single elements (Ag, Co, etc. in solid form),
chemical compounds (UO\ :sub:`2`, Gd\ :sub:`2`\O\ :sub:`3`, etc.),
or even multiple layers of both types.

The energy dependent cross-section data used in this library are from
`National Nuclear Data Center <http://www.nndc.bnl.gov/>`__, a published
online database. `Evaluated Nuclear Data File
(ENDF/B) <http://www.nndc.bnl.gov/exfor/endf00.jsp>`__ [1] is currently
supported and more evaluated databases will be added in future.

Python packages used are: SciPy [2], NumPy [3], Matplotlib [4], Pandas
[5] lmfit [6] and ImagingReso [7].

Statement of need
#################

Neutron imaging is a powerful tool to characterize material
non-destructively. And based on the unique resonance features,
it is feasible to identify elements and/or isotopes resonance with
incident neutrons. However, a dedicated user-friendly fitting tool
for resonance imaging is missing, and **ResoFit** we presented here
could fill this gap.

Installation instructions
#########################

Install **ResoFit** by typing the following command in Terminal:

``pip install ResoFit``

or by typing the following command under downloaded directory in
Terminal:

``python setup.py``

Example usage
#############

Example of usage is presented in ``tutorial.ipynb`` under ``/notebooks``
directory.

Calculation algorithm
#####################

The neutron transmission calculation algorithm of neutron transmission
*T*\ (*E*), is base on Beer-lambert law [7]-[9]:

.. math:: T\left( E \right) =\frac { I\left( E \right)  }{ { I }_{ 0 }\left( E \right)  } =exp\left[ -\sum\nolimits_i { { N }_{ i }{ d }_{ i } } \sum\nolimits_j { { \sigma  }_{ ij }\left( E \right) { A }_{ ij } }  \right]

:math:`N_i` : number of atoms per unit volume of element :math:`i`,

:math:`d_i` : effective thickness along the neutron path of element :math:`i`,

:math:`\sigma_{ij}\left( E \right)` : energy-dependent neutron total cross-section for the isotope :math:`j` of element :math:`i`,

:math:`A_{ij}` : abundance for the isotope :math:`j` of element :math:`i`.

For solid materials the number of atoms per unit volume can be
calculated from:

.. math:: {N_i} = {N_A}{C_i} = \frac{{{N_A}{\rho _i}}}{{\sum\nolimits_j {{m_{ij}}{A_{ij}}} }}

:math:`N_A` : Avogadro’s number,

:math:`C_i` : molar concentration of element :math:`i`,

:math:`\rho_i` : density of the element :math:`i`,

:math:`m_{ij}` : atomic mass values for the isotope :math:`j` of element :math:`i`.

