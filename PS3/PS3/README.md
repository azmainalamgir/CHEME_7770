# ChemE 7770 PS3
Open attached PDF to view matrix calculations and discussion.

## Part a
The stoichiometric array was constructed using a total of 18 metabolites and 15 external fluxes, in addition to the 5 fluxes of interest in the system. Flux *v5* of the reversible reaction was split into a forward and backward reaction, and the external water flux was accordingly split to either enter or exit the system, such that the direction of *v5* determined if there was net flux of water either into or out of the system.

## Part b
Running the following script:
  Balance.jl
calculates the output of the matrix multiplication of transpose of the atom matrix multiplied by the stoichiometric matrix. The first 6 columns of the output matrix are 0, indicating that the system reconstruction was elementally balanced.

## Part c
Running the following scripts:
  FBA.jl
  Flux.jl
determines the maximum rate of urea production. Setting appropriate bounds for metabolite and flux values and optimizing for urea production, the maximized value for *b4* is calculated to be **1.27 mmol/gDW-hr**.
