using CSV

#importing stoichiometric and atom matrices
S = Matrix{Float64}(CSV.read("Stoichiometry.csv", header=0))
A = Matrix{Float64}(CSV.read("Atom.csv", header=0))

#calculating A^T*S and analyzing if first 6 columns are zero
epsilon = transpose(A)*S
