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
[5] Periodictable [6], lmfit [7] and ImagingReso [8].

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
*T*\ (*E*), is base on Beer-lambert law [9]-[10]:

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

References
##########

[1] M. B. Chadwick et al., “ENDF/B-VII.1 Nuclear Data for Science and
Technology: Cross Sections, Covariances, Fission Product Yields and
Decay Data,” Nuclear Data Sheets, vol. 112, no. 12, pp. 2887–2996, Dec.
2011.

[2] T. E. Oliphant, “SciPy: Open Source Scientific Tools for Python,”
Computing in Science and Engineering, vol. 9. pp. 10–20, 2007.

[3] S. van der Walt et al., “The NumPy Array: A Structure for Efficient
Numerical Computation,” Computing in Science & Engineering, vol. 13, no.
2, pp. 22–30, Mar. 2011.

[4] J. D. Hunter, “Matplotlib: A 2D Graphics Environment,” Computing in
Science & Engineering, vol. 9, no. 3, pp. 90–95, May 2007.

[5] W. McKinney, “Data Structures for Statistical Computing in Python,”
in Proceedings of the 9th Python in Science Conference, 2010, pp. 51–56.

[6] P. A. Kienzle, “Periodictable V1.5.0,” Journal of Open Source
Software, Jan. 2017.

[7] M. Newville, A. Nelson, A. Ingargiola, T. Stensitzki, R. Otten,
D. Allan, Michał, Glenn, Y. Ram, MerlinSmiles, L. Li, G. Pasquevich,
C. Deil, D.M. Fobes, Stuermer, A. Beelen, O. Frost, A. Stark, T. Spillane,
S. Caldwell, A. Polloreno, stonebig, P.A. Brodtkorb, N. Earl, colgan,
R. Clarken, K. Anagnostopoulos, B. Gamari, A. Almarza, lmfit/lmfit-py 0.9.7,
(2017). doi:10.5281/zenodo.802298.

[8] Y. Zhang and J. C. Bilheux, "ImagingReso".

[9] M. Ooi et al., “Neutron Resonance Imaging of a Au-In-Cd Alloy for
the JSNS,” Physics Procedia, vol. 43, pp. 337–342, 2013.

[10] A. S. Tremsin et al., “Non-Contact Measurement of Partial Gas
Pressure and Distribution of Elemental Composition Using Energy-Resolved
Neutron Imaging,” AIP Advances, vol. 7, no. 1, p. 15315, 2017.

