from pyscf import gto, scf, cc, fci



#gets the coefficients for a restricted hartree fock via mean field theory, needs xyz inputs
def getHFCoefficients(mol, basisSet):
    mol = gto.M(
        atom = mol, 
        basis = basisSet, 
        symmetry=False,
    )
    m = scf.RHF(mol)
    m.kernel()
    norb = m.mo_energy.size
    return m.mo_coeff



