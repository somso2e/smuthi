This folder contains the NFM-DS Fortran program written by the authors of:

[1] A. Doicu, T. Wriedt, and Y. Eremin: "Light Scattering by Systems of Particles", Berlin, Heidelberg: Springer-Verlag, 2006.

For more information, see: https://scattport.org

We thank the authors for allowing us to use their routines.

For the purpose of integrating the code for the generation of the T-matrix of axisymmetric particles into the SMUTHI
project, the file TAXSYM.f90 was modified and saved as TAXSYM_SMUTHI.f90.

Summary of changes made to enable using TAXSYM and TNONAXSYM in F2Py module:
in TAXSYM_SMUTHI.f90:
TAXSYM subroutine is modified to allocate T-matrix (F2Py does not support allocatable types well) + F2Py statements are included to enable default arguments and input/output intent
convergence_MrankAXSYM and convergence_MrankDSAXSYM subroutines modified to output the T matrix instead of writting it to a file, Mrank is fixed at Nrank
additional subroutine appendtotmat added to facilitate filling the T matrix
in TNONAXSYM.f90:
subroutine TMatrix_Nrank_MrankNONAXSYM is modified to output the T matrix instead of writting it to a file
TNONAXSYM subroutine is modified to allocate T-matrix (F2Py does not support allocatable well) + F2Py statements are included to enable default arguments and input/output intent
Input parameters for integration, interpolation and linear system solver are know fixed to enable using the F2Py without reading input from file. This required modifying Integr.f90, Interp.f90 and MatrixSolv.f90 accordingly. EFMED.f90 is unmodified as we do not use the NFMDS implementation of the effective medium theory.

Summary of changes made to enable using TLAY in F2Py module:
additional subroutine tlayappendtotmat added to facilitate filling the T matrix
subroutines are modified to output the T matrix instead of writing it to file
the subroutine that reads the inputfile is modified to fill the default parameters
several variables are modified from allocatable to fixed size for handling with F2Py

Summary of changes made to facilitate future changes:
several unused NFMDS mode are moved to auxfiles folder
convergence testing is removed altogether
input reading is removed

