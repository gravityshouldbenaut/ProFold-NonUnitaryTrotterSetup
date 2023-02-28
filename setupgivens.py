from basissetsetup import getHFCoefficients
from basissetselector import find_optimal_basis_sets
from numpy import linalg as LA
import numpy as np
import rdkit
from rdkit.Chem import AllChem as Chem
from vqecircuit import fourQubitVQE
import cirq

#need to use obabel to get the xyz using obabel because RDKit xyz converter leaves out hydrogens 
def getHFCoefficientsBasedOnBestBasisSets(molFile, name):
	outputJSON = find_optimal_basis_sets(molFile, name, 95, "False")
	basisSet = outputJSON["Best Basis Sets Based on Ground State Precision"][0]
	print(basisSet)
	fileName = name + ".xyz"
	ouptutCoefficients = getHFCoefficients(fileName, basisSet)
	return ouptutCoefficients

def getGuessEigenvalue(matrix):
	w, v = LA.eig(matrix)
	return w

def getDeterminant(matrix):
	d = LA.det(matrix)
	return d

def getAngleForVQEFromHamiltonMatrixDeterminant(matrix, noOfAtoms): 
	d = getDeterminant(matrix)
	wfAngle = 1/math.sqrt(math.factorial(noOfAtoms)) * d
	return wfAngle

#gets the exponential setup for a given Hamiltonian matrix via Trotter decomposition
def getReciprocalDiagonal(matrix):
	return 1/(np.sum(np.diag(matrix)))


def setupVQEForHamiltonianMatrix(fileName, proteinName): 
	hfMatrix = getHFCoefficientsBasedOnBestBasisSets(fileName, proteinName)
	reciprocalDiagonal = getReciprocalDiagonal(hfMatrix)
	if(".pdb" in fileName):
		mol = Chem.MolFromPDBFile(fileName)
	else:	
		mol = Chem.MolFromSmiles(fileName)
	numAtoms = 0 
	for atom in mol.GetAtoms():
		numAtoms +=1 
	#addition of the energies found from each eigenvalue of each partition would equal total energy
	#since each subdivision is essentially running on the same circuits, you could multiply the outputs
	potentialPartitions = numAtoms / 4 
	vqe_cir = cirq.Circuit()
	q0, q1, q2, q3 = cirq.LineQubit.range(4)

	vqe_cir.append(fourQubitVQE(theta = reciprocalDiagonal).on(q0, q1, q2, q3)._decompose_())
    
	return vqe_cir



print(setupVQEForHamiltonianMatrix("proline.pdb", "proline"))
