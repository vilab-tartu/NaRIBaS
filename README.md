# NaRIBaS

NaRIBaS is Unix-based preparation and analysis toolbox for Nanomaterials and Room Temperature Ionic Liquids in Bulk and Slab.

For details see the publication https://doi.org/10.3390/computation6040057 (open access).

NaRIBaS provides a framework that decouples user input parameter and terminal based command lines. NaRIBaS does not replace a simulation software and specific analysis tools like Gromacs, but it allows iterative repetition of tasks while changing specific input parameter.

The toolbox is to be understood as a bash-scripting framework rather than a black-box software.

Task specific input and preparation or analysis function, are easy expandable to meet the scientific characteristic of constant changing properties.

Concepts for efficiently storing and analysing large files in conjunction with graphical visualization and numerical data processing are presented in a Manual.

A documented script collection that allows reproducing all simulations and analysis results presented in the Doctoral thesis of the author.

# Requirements

Any computational code which could be run using a command line.

Bash with local settings like `export LC_NUMERIC="en_US.UTF-8"`

Python (for the analyses)

[Packmol](http://www.ime.unicamp.br/~martinez/packmol/)

# References

E. Roos Nerut, K. Karu, I. Voroshylova, K. Kirchner, T. Kirchner, M. Fedorov, V. Ivaništšev, NaRIBaS — A Scripting Framework for Computational Modeling of Nanomaterials and Room Temperature Ionic Liquids in Bulk and Slab, Computation. 6 (2018) 57. doi:10.3390/computation6040057.
