# ProFold-NonUnitaryTrotterSetup
A tool for finding the right basis set for a protein, finding its Hamiltonian matrix, and embedding onto a quantum circuit for VQE via a Givens Rotation of the Trotter decomposition 

#Motivation#
Inspired by the team at Boehringer-Ingelheim who developed a SAPT-VQE method for modeling chemical systems (https://pubs.rsc.org/en/content/articlelanding/2022/sc/d1sc05691c?ref=banner), the members of the Quantum Computing Club at Davis and their affiliates have been working on a protein folding system that utilizes a state-of-the-art quantum circuit for quantum mechanical modeling along with state-of-the-art classical molecular dynamics in a combined QM/MM setup (https://github.com/QC-at-Davis/ProFold).

#Implementation#
Towards that goal, we wanted to make it easier to go from protein structures to Hamiltonian setup on a quantum device by utilzing the Givens Rotation setup used in the Boehringer paper (https://arxiv.org/pdf/2004.04174v4.pdf) for someone who does not know where to start with finding the ground state of their protein structure. To do this, we: 
  1. Estimated the right basis set between cc-pVTZ, 6-31G, and aug-CC-pVDZ given some protein file input (.pdb). 
  2. Utilized pySCF to find the right Hamiltonian setup for Hartree-Fock representations of the basis set that was picked. 
  3. Decompose matrix to Givens Rotations through Trotter decomposition only since the Givens can handle the QR decomposition into unitary.
  4. Setup Givens Rotation based VQE circuit with pair exchanges 
  
 To run this project, please run the setupVQEForHamiltonianMatrix function within setupGivens for the example file, Proline. 
 
 #Goals#
 We hope to implement this method with a quantum-based optimizer sometime soon in order to actually run the VQE for many iterations, and that implementation on hardware would lead to fewer errors than Trotter-Suzuki due to further padding provided by Givens structures.

