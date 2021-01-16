!***********************************************************************************
subroutine printinputNONAXSYM (FileGeom, TypeGeom, FileFEM, Nsurf, Nparam,          &
           Nazimutsym, dNint1, dNint2, ind_refMed, wavelength, anorm, Rcirc, surf,  &
           kb, epsNint, epsNrank, epsMrank, ind_refRel, ind_refRelZ, alphaPR,       &
           betaPR, anisotropic, miror, perfectcond, chiral)
  use parameters
  implicit none
  integer        :: TypeGeom, Nsurf, Nparam, Nazimutsym, dNint1, dNint2, i, LenString                     
  real(O)        :: ind_refMed, wavelength, anorm, surf(Nsurf), kb, epsNint,        &
                    epsNrank, epsMrank, Rcirc, alphaPR, betaPR
  complex(O)     :: ind_refRel, ind_refRelZ
  character(80)  :: FileFEM, FileFEMWrite
  logical        :: FileGeom, miror, perfectcond, chiral, anisotropic
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium,ind_refMed = ', ind_refMed, ';'
  if (.not. anisotropic) then
    write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                             &
   'relative refractive index of the particle, ind_refRel = (', ind_refRel, ');'
  else
    write (iOutput,"(2x,'relative refractive indices of the uniaxial particle:')")
    write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a, 1pe10.3, ',', 1pe10.3, a)")   &
   'ind_refRel = (', ind_refRel, '), ind_refRelZ = (', ind_refRelZ, ');' 
  end if
  write (iOutput,*)
  if (FileGeom) then
    FileFEMWrite = FileFEM(14:LenString(FileFEM))
    write (iOutput,"(2x, a, a30)")                                                  &
   'name of the file containing the particle geometry, FileFEM = ',     FileFEMWrite
  else
    write (iOutput,"(2x,'index of the particle geometry, TypeGeom = ',i2,';')")     &
           TypeGeom
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'ellipsoid;')")       
    else if (TypeGeom == 2) then
      write (iOutput,"(2x,'quadratic prism;')") 
    else if (TypeGeom == 3) then
      write (iOutput,"(2x,'regular N-hedral prism;')")
    end if
    write (iOutput,"(2x,'number of surface parameters, Nsurf = ',i2,';')") Nsurf
    write (iOutput,"(2x,'surface parameters:')")
    do i = 1, Nsurf
      write (iOutput,"(2x,'surf(',i2,') = ',1pe10.3,',')") i, surf(i)
    end do       
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'where surf(1) is the semi-axis along the x-axis,')")
      write (iOutput,"(2x,'surf(2) is the semi-axis along the y-axis and')")
      write (iOutput,"(2x,'surf(3) is the semi-axis along the z-axis;')")
    else if (TypeGeom == 2 .or. TypeGeom == 3) then
      write (iOutput,"(2x,'where surf(1) is the half-length of the prism,')")
      write (iOutput,"(2x,'and   surf(2) is the half-length of the basis side;')")
    end if 
    write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')")       &
           Nparam
  end if
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'characteristic length of the particle, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc = ', Rcirc, ';'
  if (anisotropic) then
    write (iOutput,"(2x, a)")                                                       &
   'Euler angles specifying the orientation of the principal coordinate system:'
    write (iOutput,"(2x, a, 1pe10.3, a, 1pe10.3, a, 1pe10.3, a)")                   &
   'alpha = ', alphaPR * 180._O / Pi, ', beta = ', betaPR * 180._O / Pi,            &
   ', gamma = ', 0._O, ';'        
  end if
  if (miror) write (iOutput,"(2x,'mirror symmetric particle;')")
  write (iOutput,"(2x, a, i2, a)")                                                  &
 'number of azimuthal symmetric sections, Nazimutsym = ', Nazimutsym, ';'
  if (Nazimutsym >=2) then
    write (iOutput,"(2x,'azimuthal symmetry is taken into account;')")
  else
    write (iOutput,"(2x,'azimuthal symmetry is disregarded;')")
  end if
  write (iOutput,*)
  if (perfectcond) then
    write (iOutput,"(2x,'perfectly conducting particle;')")
  else if (chiral) then
    write (iOutput,"(2x,'chiral particle;')")
    write (iOutput,"(2x,'characteristic of chirality, kb = ',1pe10.3,';')") kb  
  else if (anisotropic) then
    write (iOutput,"(2x,'uniaxial anisotropic particle;')")
  else 
    write (iOutput,"(2x,'dielectric particle;')")
  end if 
  write (iOutput,*)    
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'
  if (.not. FileGeom) then               
    write (iOutput,"(2x,'integration tolerance, epsNint = ',1pe10.3,';')") epsNint                         
    write (iOutput,"(2x,'integration steps, dNint1 = ',i4,', dNint2 = ',i4,'.')")   &
           dNint1, dNint2 
  end if
  write (iOutput,"(/)")                       
end subroutine printinputNONAXSYM
!***********************************************************************************
subroutine TMatrix_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_ref, Nsurf,    &
           surf,rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, miror,     &
           Nazimutsym, perfectcond, chiral, kb, FileTmat, PrnProgress,Nmax,b)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,      &
                   Nazimutsym                   
  real(O)       :: k, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD), kb
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, chiral, PrnProgress 
  character(80) :: FileTmat
!      
  integer       :: Nmax, i, j, NintAL
  integer,allocatable    :: Nintparam(:) 
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:)
  complex(O) :: a(2*Nmax,2*Nmax), b(2*Nmax,2*Nmax)
!  
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NintAL = max(Nint1,Nint2)  
  allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),           &
            weightsG(Nparam,NintAL*NintAL))
  allocate (Nintparam(Nparam))
  if (.not. FileGeom) then    
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
         Nintparam, paramG1, paramG2, weightsG, miror, Nazimutsym)                            
  else
    do i = 1, Nparam
      do j = 1, NintAL*NintAL
        paramG1(i,j)  = 0._O
        paramG2(i,j)  = 0._O
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do   
  end if
  if (PrnProgress) call write_progress (.true., 1, 4)
  call matrix_Q (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area,   &
       Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,      &
       weightsG, miror, Nazimutsym, perfectcond, chiral, kb, a, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 2, 4)  
  call matrix_Q (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area,   &
       Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,      &
       weightsG, miror, Nazimutsym, perfectcond, chiral, kb, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)  
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (paramG1, paramG2, weightsG, Nintparam)      
end subroutine TMatrix_Nrank_MrankNONAXSYM 
!***********************************************************************************
subroutine TMatrix_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_ref, ind_refZ,       &
           alfaPR, betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, &
           Nrank, Nbeta, Nint1, Nint2, FileTmat, PrnProgress,Nmax,b)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, Nbeta
  real(O)       :: k, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD),     &
                   alfaPR, betaPR, gamaPR
  complex(O)    :: ind_ref, ind_refZ
  character(80) :: FileTmat
  logical       :: FileGeom, PrnProgress
!        
  integer       :: Nmax, NintAL  
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:)
  complex(O)             :: a(2*Nmax,2*Nmax), b(2*Nmax,2*Nmax) 
!  
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)   
  NintAL = max(Nint1,Nint2)  
!  open (unit = iTmat, file = FileTmat, status = 'replace')
!  call  write_HeadFileTmat (Nmax, Nmax)    
!  allocate (a(2*Nmax,2*Nmax), b(2*Nmax,2*Nmax))                          
  allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),           &
            weightsG(Nparam,NintAL*NintAL))
  allocate (Nintparam(Nparam))
  call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam,   &
       Nintparam, paramG1, paramG2, weightsG,.false., 1)   
  if (PrnProgress) call write_progress (.true., 1, 4)
  call matrix_Q_anis (FileGeom, TypeGeom, 3, 1, k, ind_ref, ind_refZ, alfaPR,       &
       betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax, Nbeta, &
       NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, a, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 2, 4)   
  call matrix_Q_anis (FileGeom, TypeGeom, 1, 1, k, ind_ref, ind_refZ, alfaPR,       &
       betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax, Nbeta, &
       NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)  
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
!  call write_FileTmat (Nmax, Nmax, b)
!  close (unit = iTmat)
!  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., .false.)
!  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., .false.)
  print "(/,2x,'The T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,';')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (paramG1, paramG2, weightsG, Nintparam)      
end subroutine TMatrix_Nrank_MrankAnis
subroutine TNONAXSYM (wavelength, ind_refMed, ind_refRel, ind_RefRelZ,            &
       alphaPR, betaPR, perfectcond, anisotropic, chiral, kb, FileGeom,             &
       TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm, Rcirc, miror, Nazimutsym,     &   
       DoConvTest, ExtThetaDom, Nbeta, Nint1, Nint2, Nrank, Mrank, epsNint,         &
       epsNrank, epsMrank, dNint1, dNint2, PrnProgress, gammaPR, Nface,Nmax,b) 
  !k, snorm, rp, np, area,  TypeConvTest      
  use parameters  
  use derived_parameters  
  implicit none 
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nazimutsym, TypeConvTest, Mrank, &
                   Nrank, Nbeta, Nint1, Nint2, dNint1, dNint2, Nmax                     
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(Nsurf), snorm,            &
                   kb, epsNint, epsNrank, epsMrank, Rcirc, rp(3,NfacePD),           &
                   np(3,NfacePD), area(NfacePD), alphaPR, betaPR, gammaPR
  complex(O)    :: ind_refRel, ind_RefRelZ
  logical       :: FileGeom, miror, perfectcond, chiral, anisotropic, DoConvTest,   &
                   ExtThetaDom, PrnProgress
  character(80) :: FileTmat, FileFEM     
  integer       :: NrankW, i, j, ios, icall                     
  real(O)       :: xpart, x, dp, grd
  logical       :: more, continuare, IntTest, XFindPar
  complex(O), intent(out)    :: b(2*Nmax,2*Nmax) 
  !f2py real(O) :: wavelength  = 1
  !f2py real(O) :: ind_refMed  = 1
  !f2py real(O) :: kb = 0
  !f2py character(80) :: FileFEM  = '../GEOMFILES/cubekc.fem'
  !f2py logical :: miror  = 0 
  !f2py logical :: FileGeom = 1   
  !f2py integer :: TypeGeom = 1  
  !f2py integer :: Nparam = 1 
  !f2py real(O) :: anorm  = 1
  !f2py real(O) :: Rcirc  = 1
  !f2py integer :: Nazimutsym = 0  
  !f2py integer :: Nface = 1
  !f2py complex(O) :: ind_refRelZ = (1.5,0.)
  !f2py complex(O) :: ind_refRel  = (4,0.)  
  !f2py real(O) :: alphaPR = 0
  !f2py real(O) :: betaPR  = 0
  !f2py integer :: Nbeta   = 60
  !f2py integer :: Nrank  = 10
  !f2py integer :: Mrank  = 10
  !f2py integer :: Nint1 = 100
  !f2py integer :: Nint2 = 100          
  !f2py real(O) :: epsNint  = 5.e-2
  !f2py real(O) :: epsNrank = 5.e-2
  !f2py real(O) :: epsMrank = 5.e-2
  !f2py integer :: dNint1 = 4
  !f2py integer :: dNint2 = 4
  !f2py real(O) :: gammaPR = 0 
  !f2py logical :: perfectcond = 0 
  !f2py logical :: anisotropic = 0  
  !f2py logical :: chiral      = 0
  !f2py logical :: DoConvTest   = 0
  !f2py logical :: ExtThetaDom  = 1          
  !f2py logical :: PrnProgress = 1 
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputNONAXSYM                      ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters        
  k   = 2._O * Pi * ind_refMed / wavelength
  grd = Pi / 180._O    
  call check_MatPropNONAXSYM (perfectcond, anisotropic, chiral, kb) 
  if (chiral) call check_chirality (kb)        
  call check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral, anisotropic,          &
       Nazimutsym)
  call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)    
  call check_anorm (anorm)
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
  do i = 1, NfacePD
    do j = 1, 3
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do 
  if (FileGeom) then    
    call read_FileFEM (FileFEM, Nface, rp, np, area) 
    Rcirc = 0._O
    do i = 1, Nface
      dp = sqrt(rp(1,i)**2 + rp(2,i)**2 + rp(3,i)**2)
      if (dp > Rcirc) Rcirc = dp
    end do            
  end if                      
  if (anisotropic) call check_ind_ref3 (ind_refRel, ind_refRelZ)  
  alphaPR = alphaPR * grd
  betaPR  = betaPR  * grd    
  x = k * Rcirc
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O) 
  call check_MrankNrank (Mrank, Nrank)  
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------                                    
  if (.not. anisotropic) then                
    call TMatrix_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_refRel, Nsurf, &
         surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, miror,    &
         Nazimutsym, perfectcond, chiral, kb, FileTmat, PrnProgress,Nmax,b)                                                     
  else
    call TMatrix_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_refRel, ind_refRelZ, &
         alphaPR, betaPR, gammaPR, Nsurf, surf, rp, np, area, Nface, Nparam,      &
         Mrank, Nrank, Nbeta, Nint1, Nint2, FileTmat, PrnProgress,Nmax,b) 
  end if 
end subroutine TNONAXSYM
