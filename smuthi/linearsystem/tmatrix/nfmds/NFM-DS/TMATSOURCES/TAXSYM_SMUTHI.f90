include 'Parameters.f90'
include 'Allocation.f90'
include 'AdditonTh.f90'
include 'BesLeg.f90'
include 'Check.f90'
include 'EFMED.f90'
include 'GeomLib.f90'
include 'GeomTrans.f90'
include 'IncCoeff.f90'
include 'InputOutput.f90'
include 'Integr.f90'
include 'Interp.f90'
include 'MachParam.f90'
include 'MatrixOp.f90'
include 'MatrixQ.f90'
include 'MatrixSolv.f90'
include 'MatrixTrans.f90'
include 'PostProces1.f90'
include 'PostProces2.f90'
include 'PostProces3.f90'
include 'Proces1.f90'
include 'Proces2.f90'
include 'Proces3.f90'
include 'Random.f90'
include 'SCT.f90'
include 'SCTAVRGSPH.f90'
include 'SVWF.f90'
include 'TLAY.f90'
include 'TNONAXSYM.f90'


!***********************************************************************************
subroutine appendtotmat (b,btot,m,Nrank,Nmax,Nmaxmax)
  use parameters
  implicit none   
  integer       :: Nrank, m, Nmax, Nmaxmax, k1, k2, N0, j, &
                   mmax
  complex(O)    :: btot(2*Nmaxmax,2*Nmaxmax)
  complex(O)    :: b(2*Nrank,2*Nrank)
!
  if (m==0) then    
    N0 = 0    
    mmax = Nrank
    do k1 = 1, mmax
      do k2 = 1, Nrank
        btot(k1+N0,k2+N0) = b(k1,k2)
        btot(k1+N0+Nmaxmax,k2+N0+Nmaxmax) = b(k1+Nmax,k2+Nmax)
        btot(k1+N0,k2+N0+Nmaxmax) = b(k1,k2+Nmax)
        btot(k1+N0+Nmaxmax,k2+N0) = b(k1+Nmax,k2)
      end do
    end do
  else
    N0 = Nrank + (m - 1) * (2 * Nrank - m + 2)
    mmax = Nrank - m + 1
    do j = 1,2
      do k1 = 1, mmax
        do k2 = 1, mmax          
          btot(k1+N0,k2+N0) = b(k1,k2)
          btot(k1+N0+Nmaxmax,k2+N0+Nmaxmax) = b(k1+Nmax,k2+Nmax)
          if (j==2) then
            b(k1,k2+Nmax)=-b(k1,k2+Nmax)
            b(k1+Nmax,k2)=-b(k1+Nmax,k2)
          end if
          btot(k1+N0,k2+N0+Nmaxmax) = b(k1,k2+Nmax)
          btot(k1+N0+Nmaxmax,k2+N0) = b(k1+Nmax,k2)
        end do
      end do
      N0 = N0 + Nrank - m + 1
    end do  
  end if
end subroutine appendtotmat
!***********************************************************************************
subroutine tmatrix_MrankAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,    &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,         &
           perfectcond, DS, chiral, kb, epsMrank, FileTmat, PrnProgress,Nmaxmax,btot)                
  use parameters
  implicit none   
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint
  real(O)       :: k, surf(Nsurf), snorm, epsMrank, kb, zRe(Nrank), zIm(Nrank),     &
                   rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
  character(80) :: FileTmat
!      
  integer       :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv, Nprog, j
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext  
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), c(:), c1(:), cc(:)
  complex(O)    :: btot(2*Nmaxmax,2*Nmaxmax)
  complex(O)    :: b(2*Nrank,2*Nrank)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 0
!  Mstart = Nrank
  Mrank  = Nrank      
!  open (unit = iTmat, file = FileTmat, status = 'replace')  
!  call write_HeadFileTmat (Nrank, Nrank) 
!  call write_TypeConvHead (3)
!  call write_2ConvParamAxsym (Nint, Nrank)
  allocate (a(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do  
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  if (.not. FileGeom) then
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)
  else
    do i = 1, Nparam
      do j = 1, Nint
        paramG(i,j)   = 0._O 
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do
  end if        
  if (.not. chiral) then
    Nprog = 4
  else
    Nprog = 7
  end if
  Mrank = - 1
!  Mrank = Mrank - 1
  do m = Mstart, Nrank       
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if  
    if (PrnProgress) call write_progress_m (.true., m, 1, Nprog)    
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank) 
    if (PrnProgress) call write_progress_m (.false., m, 2, Nprog)       
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank) 
    if (PrnProgress) call write_progress_m (.false., m, 3, Nprog)                           
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress_m (.false., m, 4, Nprog)                   
!    call write_FileTmat (Nrank, Nrank, b)
    call appendtotmat (b,btot,m,Nrank,Nmax,Nmaxmax)     
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)       
    if (m /= 0) then
      if (.not. chiral) then
        call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)        
      else 
        call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank,Nmax, Nint, Nparam, Nintparam,        &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
        if (PrnProgress) call write_progress_m (.false., m, 5, Nprog)                     
        call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank,Nmax, Nint, Nparam, Nintparam,        &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)    
        if (PrnProgress) call write_progress_m (.false., m, 6, Nprog)
        call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)       
        if (PrnProgress) call write_progress_m (.false., m, 7, Nprog)
!        call write_FileTmat (Nrank, Nrank, b)
        call appendtotmat (b,btot,-m,Nrank,Nmax,Nmaxmax)
      end if
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c) 
      call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)     
    end if
!    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
!        .false.,.true., h, v)
!    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)     
!    call write_DSCS (Nteta,.false., h, v)
!    if (NthetaConv >= int(0.8*Nteta)) exit 
  end do
!  close (unit = iTmat) 
!  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
!  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
!       alfap, k, snorm, Cext, Qext)
!  call write_Effic (Qscat, Qext)
!  call write_MrankConvRes (NthetaConv, epsMrank)
!  if (NthetaConv >= int(0.8*Nteta)) then
!    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"                                              
!  else
!    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
!  end if
!  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .true., .false., chiral)
!  call ScatCharact (k, FileTmat, Mrank, Nrank, .true., .false., chiral)
!  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"      
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank                     
  deallocate (a, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine tmatrix_MrankAXSYM
!***********************************************************************************
subroutine tmatrix_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,  &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,         &
           perfectcond, DS, chiral, kb, epsMrank, dNintMrank, FileTmat, PrnProgress,Nmaxmax,btot)            
  use parameters
  implicit none   
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint, dNintMrank
  real(O)       :: k, surf(Nsurf), snorm, epsMrank, kb, zRe(Nrank), zIm(Nrank),     &
                   rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
  character(80) :: FileTmat
!      
  integer       :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv, Nprog, j
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext  
  logical       :: ComplexPlane
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:),  c(:), c1(:), cc(:)
  complex(O)    :: btot(2*Nmaxmax,2*Nmaxmax)
  complex(O)    :: b(2*Nrank,2*Nrank)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 0
  Mrank  = Nrank      
  allocate (a(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do  
  allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
  if (.not. FileGeom) then
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)
  else
    do i = 1, Nparam
      do j = 1, Nint
        paramG(i,j)   = 0._O 
        weightsG(i,j) = 0._O
      end do
      Nintparam(i) = 1
    end do
  end if       
  ComplexPlane = .false.        
  do i = 1, Nrank                  
    if (zIm(i) /= 0._O) ComplexPlane = .true.
  end do
  if (.not. chiral) then
    Nprog = 5
  else
    Nprog = 9
  end if                        
  Mrank = - 1
  do m = Mstart, Nrank          
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    if (m == 0) then
      Nmax = Nrank
    else
      Nmax = Nrank - m + 1
    end if  
    if (ComplexPlane) then
      deallocate (paramG, weightsG, Nintparam)
      Nint = Nint + dNintMrank
      allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
      call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,            &
           Nintparam, paramG, weightsG, miror)          
    end if 
    if (PrnProgress) call write_progress_m (.true., m, 1, Nprog)  
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 2, Nprog) 
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,&
         zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,    &
         Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 3, Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress_m (.false., m, 4, Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,&
         2*Nrank)   
    if (PrnProgress) call write_progress_m (.false., m, 5, Nprog)
    call appendtotmat (a,btot,m,Nrank,Nmax,Nmaxmax)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)           
    call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)       
    if (m /= 0) then
      if (.not. chiral) then
        call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
      else 
        call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam,       &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
        if (PrnProgress) call write_progress_m (.false., m, 6, Nprog)
        call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area,   &
             Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,     &
             weightsG, b, Nrank, Nrank)
        if (PrnProgress) call write_progress_m (.false., m, 7, Nprog)
        call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank,   &
             2*Nmax)
        if (PrnProgress) call write_progress_m (.false., m, 8, Nprog)
        call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, &
             area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam,       &
             paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
        call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank,        &
             b, 2*Nrank, 2*Nrank)               
        if (PrnProgress) call write_progress_m (.false., m, 9, Nprog)
        call appendtotmat (a,btot,-m,Nrank,Nmax,Nmaxmax)
      end if
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c)
      call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
    end if
  end do
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank 
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank             
  deallocate (a, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine tmatrix_MrankDSAXSYM
subroutine TAXSYM ( wavelength, ind_refMed, ind_refRel, perfectcond,                &
           chiral, kb, FileGeom, TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm,     &
           Rcirc, miror, DoConvTest, MishConvTest, DS, autGenDS, ComplexPlane,      &
           epsZReIm, Nint, Nrank, epsNint, epsNrank, epsMrank, dNint, dNintMrank,   &
           FileTmat, PrnProgress, Nface,Nmaxmax,b) 
  use parameters 
  use derived_parameters  
  implicit none 
  integer       :: TypeGeom, Nsurf, Nparam, TypeConvTest, dNint, Nrank, Nint,       &
                   Nrank1, dNintMrank, Ndgs, NrankMax, NrankW, Nface, i, j, Mrank,  &
                   Nmaxmax
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(Nsurf), snorm,            &
                   kb, epsNint, epsNrank, epsMrank, zRe(NrankPD), zIm(NrankPD),     &
                   zRe1(NrankPD), zIm1(NrankPD), Rcirc, x, delta, Cscat1, Cext1,    &
                   epsZReIm, rp(2,NfacePD), np(2,NfacePD), area(NfacePD)                          
  complex(O)    :: ind_refRel
  logical       :: FileGeom, miror, perfectcond, DoConvTest, DS, chiral, autGenDS,  &
                   ComplexPlane, MishConvTest, PrnProgress
  character(80) :: FileTmat, FileFEM
  complex(O), intent(out)    :: b(2*Nmaxmax,2*Nmaxmax)
  !f2py real(O) :: wavelength  = 1
  !f2py real(O) :: ind_refMed = 1                                                                 
  !f2py complex(O) :: ind_refRel = (1.5,0) 
  !f2py real(O) :: epsNint  = 5.e-2
  !f2py real(O) :: epsNrank = 5.e-2
  !f2py real(O) :: epsMrank = 5.e-2
  !f2py integer :: dNint    = 4
  !f2py integer :: dNintMrank = 10
  !f2py character(80) :: FileTmat = '../TMATFILES/T.dat'
  !f2py logical :: PrnProgress = 1 
  !f2py logical :: perfectcond = 0  
  !f2py logical :: chiral = 0
  !f2py real(O) :: kb     = 0 
  !f2py logical :: FileGeom = 0
  !f2py character(80) :: FileFEM  = ' '  
  !f2py integer :: TypeGeom = 1 
  !f2py integer :: Nparam = 1
  !f2py real(O) :: anorm  = 1
  !f2py real(O) :: Rcirc  = 1 
  !f2py logical :: miror  = 0
  !f2py integer :: Nface = 1
  !f2py logical :: DoConvTest   = 0
  !f2py logical :: MishConvTest = 0
  !f2py logical :: ComplexPlane = 1
  !f2py real(O) :: epsZReIm = 0.95 
  !f2py integer :: Nint   = 100
  !f2py integer :: Nrank  =  17
  !f2py logical :: DS = 0
  !f2py logical :: autGenDS = 1 
!************************************************************************************
!  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  integer       :: ios                     
  real(O)       :: xpart,dp 
  logical       :: InputDS, XFindPar      
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputAXSYM                         ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters   
  call check_ind_ref (ind_refRel)
  k = 2._O * Pi * ind_refMed / wavelength 
!       
  call check_MatPropAXSYM (perfectcond, chiral, kb)   
  if (chiral) call check_chirality (kb)

!                          
  call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
  call check_geomAXSYMOblate (TypeGeom, Nsurf, surf)  
  call check_anorm (anorm)
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
  do i = 1, NfacePD
    do j = 1, 2
      rp(j,i) = 0._O
      np(j,i) = 0._O
    end do
    area(i) = 0._O
  end do  
  if (FileGeom) then    
    call read_FileFEMAxsym (FileFEM, Nface, rp, np, area)     
    Rcirc = 0._O    
    do i = 1, Nface
      dp = sqrt(rp(1,i)**2 + rp(2,i)**2)
      if (dp > Rcirc) Rcirc = dp
    end do    
  end if          
! 
  if (chiral)   MishConvTest = .false.     
  if (FileGeom) MishConvTest = .false.  
!    
  call check_inputAXSYM (miror, chiral, DS)
  if (FileGeom) autGenDS = .false.  
  InputDS = .false.
  if (DS .and. .not. autGenDS) InputDS = .true. 
  do i = 1, NrankPD
    zRe(i)  = 0._O
    zIm(i)  = 0._O
    zRe1(i) = 0._O
    zIm1(i) = 0._O
  end do                 
! -----------------------------------------------------------------------------------
!         Select the type of convergence test and the values of Nint and Nrank      !
! -----------------------------------------------------------------------------------
  x = k * Rcirc
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  if (.not. DS) then
!    else
    print "(/,2x,'Convergence Test for an Axisymmetric Particle over Mrank')" 
    print "(  2x,'--------------------------------------------------------')"
    print "(/,2x,'Input values:')"
    if (.not. FileGeom) then
      print "(  2x, a, i4, a, i4, a)",                                            & 
     'the input values of Nint and Nrank are ', Nint, ' and ', Nrank,             &
     ', respectively,'  
      print "(  2x, a, i3, a)",                                                   &
     'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW,';'
    else
      print "(  2x,'the input value of Nrank is ', i4,', while')", Nrank
      print "(  2x, a, i3, a)",                                                   &
     'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
    end if                                                                          
    TypeConvTest = 3      
!    end if
    Nrank1 = Nrank - 1 ! redundant            
  else
    if (autGenDS) then                       
!      else
      print "(/,2x,'Convergence Test for an Axisymmetric Particle over Mrank')" 
      print "(  2x,'--------------------------------------------------------')"
      print "(/,2x,'Input values:')"
      print "(  2x, a, i4, a, i4, a)",                                            &
     'the input values of Nint and Nrank are ', Nint, ' and ', Nrank,             &
     ', respectively,' 
      print "(  2x, a, i3, a)",                                                   &
     'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW,';'         
      call check_MaxNrank (Nrank)
      call zDSAXSYM (TypeGeom, Nsurf, surf, Nrank, ComplexPlane, epsZReIm, zRe, zIm) 
      TypeConvTest = 3    
!     end if
      Nrank1 = Nrank - 1 ! redundant  
    else 
!      else
      print "(/,2x,'Convergence Test for an Axisymmetric Particle over Mrank')" 
      print "(  2x,'--------------------------------------------------------')" 
      print "(/,2x,'Input values:')" 
      if (.not. FileGeom) then                  
        print "(  2x, a, i4, a, i4, a)",                                          &
       'the input values of Nint and Nrank are ', Nint, ' and ', Nrank,           &
       ', respectively,' 
        print "(  2x, a, i3, a)",                                                 &
       'while the estimated value of Nrank from Wiscombe''s criterion is ',       &
        NrankW,';' 
      else
        print "(  2x,'the input value of Nrank is ', i4,', while')", Nrank
        print "(  2x, a, i3, a)",                                                 &
       'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
      end if
      TypeConvTest = 3              
!      end if
      Nrank1 = Nrank - 1 ! redundant          
    end if
  end if                                                                            
! -----------------------------------------------------------------------------------
!                               T-matrix calculation                                !
! -----------------------------------------------------------------------------------
!  else      
  if (.not. DS) then    
    call tmatrix_MrankAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm, Nsurf,     &
         surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,                &
           miror, perfectcond, DS, chiral, kb, epsMrank, FileTmat, PrnProgress,   &
           Nmaxmax,b)
  else          
    call tmatrix_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,          &
         Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,  &
         perfectcond, DS, chiral, kb, epsMrank, dNintMrank, FileTmat, PrnProgress,&
         Nmaxmax,b)
  end if  
!  end if 
end subroutine TAXSYM
