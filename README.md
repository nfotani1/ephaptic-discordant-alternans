# ephaptic-discordant-alternans

Determine the effect of ephaptic coupling on discordant alternans.

- Version history (in reverse chronological order):
    - Added a spatially detailed, cylindrically-symmetric model of the effect of cleft-cell coupling during action potential propagation (in a new folder, cyl_model_of_cell_cleft_display_currents):
        * As a check, the code displays all the branches, and all the nodes to which they are connected.
        * The code produces movies containing a colorplot of the potential, and a quiver plot of the current traveling through the (linear) resistors in the model.
    - Added the taem (twice-abbreviated ephaptic model) fiber model:
        * Still has a couple of flaws: 
            (1) xi is just given a value (should be calculated).
            (2) line colors are not consistent.
    - Changed local directory location.
    - Just starting this git-based project.  So, start by modeling the subcircuit representing the dynamics of wave propagation in the cleft region (between the cells).