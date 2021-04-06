# ephaptic-discordant-alternans

Determine the effect of ephaptic coupling on discordant alternans.

- Version history (in reverse chronological order):
    - Added a four-node version of the model to the code (in the form of a new folder, four_node_cell_fiber_model), since intracellular resistance comparable to the intercell gap junction resistance should be important and non-negligible.   The code can be run in either of two modes: in the main file, four_node_cell_ephaptic_fiber.m:
        * When velocity_measurement_mode = true and plot_extended_diagnostics = false, the code executes several runs to produce a velocity plot for different values of the gap-junction coupling and cleft width.  
        * When the reverse, the code generates a membrane potential colorplot, showing the discordant alternans pattern in time and space, along with other diagnostics associated with discorant alternans.
    - Added a spatially detailed, cylindrically-symmetric model of the effect of cleft-cell coupling during action potential propagation (in a new folder, cyl_model_of_cell_cleft_display_currents):
        * As a check, the code displays all the branches, and all the nodes to which they are connected.
        * The code produces movies containing a colorplot of the potential, and a quiver plot of the current traveling through the (linear) resistors in the model.
    - Added the taem (twice-abbreviated ephaptic model) fiber model:
        * Still has a couple of flaws: 
            (1) xi is just given a value (should be calculated).
            (2) line colors are not consistent.
    - Changed local directory location.
    - Just starting this git-based project.  So, start by modeling the subcircuit representing the dynamics of wave propagation in the cleft region (between the cells).