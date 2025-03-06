# Geometric Structures from Diagrams
Python code for computing 
- the complete hyperbolic structure by giving equations for edge and crossing labels (1 in the menu) and their complex values (2 in the menu) for a hyperbolic link in 3-sphere;
- equations for the canonical component of PSL(2, C)-character variety (3 in the menu)
of a hyperbolic knot in 3-sphere.

The link/knot diagram needs to be taut (e.g. any reduced alternating diagram is taut) and can be given as Dowker-Thistlethwaite (DT) code or planar diagram (PD) code (the latter only for alternating links).

Computing  1 and 2 above are based on this paper by Thistlethwaite and Tsvietkova: www.arxiv.org/abs/1108.0510; 3 is based on an upcoming preprint by K. Petersen and A. Tsvietkova. The equations for 3 can be large, and are hence recorded in a txt file in the folder with the code.

Different people worked on the code at different times: Jaeyun Bae, Mark Bell, Dale Koenig, Alex Lowen, Anastasiia Tsvietkova.

The code uses the spherogram module from SnapPy (www.github.com/3-manifolds/Spherogram), as well as the NumPy, SciPy, and SymPy packages.

Exe file for Windows can be found here: sites.rutgers.edu/anastasiia-tsvietkova/geometric-structures-from-diagrams-code/

If you use the code, please reference it as
Jaeyun Bae, Mark Bell, Dale Koenig, Alex Lowen, Anastasiia Tsvietkova, Geometric structures from diagrams,  www.github.com/Nastia-Tsvietkova/Geometric-Structures
