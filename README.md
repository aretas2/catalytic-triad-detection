# catalytic-triad-detection
The script automates detection of the cannonical catalytic triad (Ser-His-Asp) in serine proteases based on the residue
distances and angles. In theory, this should work with any strucutre with the triad e.g. subtilisin.

The script WILL NOT detect the triad if any of the triad residues are mutated.

The script was written and tested with Python 3 and should work with Python 2.
The script accepts a PDB file as a sole input file and can be easily used with "for f in *.pdb; do ..." shell loop for multiple structures.

The output is provided in the terminal window as atom and residue numbers.
