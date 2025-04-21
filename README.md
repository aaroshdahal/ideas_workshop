# ideas_workshop

## Workflow:
1. create a mesh generator file (.geo file) and generate .msh using gmsh gui
2. convert .msh to .xml using dolfin-convert
3. .xml is the input to the python file based on finite element platform fenics -- run using batch file
4. output files: .txt and xdmf -- visualised in paraview and occasionally in spreadsheet

## References:
Aditya Kumar, Blaise Bourdin, Gilles A. Francfort, Oscar Lopez-Pamies,
Revisiting nucleation in the phase-field approach to brittle fracture,
Journal of the Mechanics and Physics of Solids,
Volume 142,
2020,
104027,
ISSN 0022-5096,
https://doi.org/10.1016/j.jmps.2020.104027.
(https://www.sciencedirect.com/science/article/pii/S0022509620302623)

Abstract: Twenty years in since their introduction, it is now plain that the regularized formulations dubbed as phase-field of the variational theory of brittle fracture of Francfort and Marigo (1998) provide a powerful macroscopic theory to describe and predict the propagation of cracks in linear elastic brittle materials under arbitrary quasistatic loading conditions. Over the past ten years, the ability of the phase-field approach to also possibly describe and predict crack nucleation has been under intense investigation. The first of two objectives of this paper is to establish that the existing phase-field approach to fracture at large — irrespectively of its particular version — cannot possibly model crack nucleation. This is so because it lacks one essential ingredient: the strength of the material. The second objective is to amend the phase-field theory in a manner such that it can model crack nucleation, be it from large pre-existing cracks, small pre-existing cracks, smooth and non-smooth boundary points, or within the bulk of structures subjected to arbitrary quasistatic loadings, while keeping undisturbed the ability of the standard phase-field formulation to model crack propagation. The central idea is to implicitly account for the presence of the inherent microscopic defects in the material — whose defining macroscopic manifestation is precisely the strength of the material — through the addition of an external driving force in the equation governing the evolution of the phase field. To illustrate the descriptive and predictive capabilities of the proposed theory, the last part of this paper presents sample simulations of experiments spanning the full range of fracture nucleation settings.
Keywords: Fracture; Strength; Energy methods; Configurational forces; Brittle materials

