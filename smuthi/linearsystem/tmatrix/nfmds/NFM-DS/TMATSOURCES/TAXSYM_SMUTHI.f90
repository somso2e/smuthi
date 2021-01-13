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

!************************************************************************************
subroutine readinputAXSYM ( wavelength, ind_refMed, ind_refRel, perfectcond,        &
           chiral, kb, FileGeom, TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm,     &
           Rcirc, miror, DoConvTest, MishConvTest, DS, autGenDS, ComplexPlane,      &
           epsZReIm, Nint, Nrank, zRe, zIm, zRe1, zIm1, epsNint, epsNrank,          &
           epsMrank, dNint, dNintMrank, FileTmat, PrnProgress, k, snorm, Nface,     &
           rp, np, area )
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, Nsurf, Nparam, dNint, Nrank, Nint, i, dNintMrank, ios, &
                   Nface, j                     
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(NsurfPD), xpart, snorm,   &
                   kb, epsNint, epsNrank, epsMrank, zRe(NrankPD), zIm(NrankPD),     &
                   zRe1(NrankPD), zIm1(NrankPD), Rcirc, epsZReIm, rp(2,NfacePD),    &
                   np(2,NfacePD), area(NfacePD), dp                          
  complex(O)    :: ind_refRel
  logical       :: FileGeom, miror, perfectcond, DoConvTest, DS, chiral, autGenDS,  &
                   ComplexPlane, InputDS, MishConvTest, PrnProgress, XFindPar
  character(80) :: FileTmat, FileFEM, string       
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputAXSYM                         ! 
! -----------------------------------------------------------------------------------    
  call DrvParameters  
  open (unit = iInputAXSYM, file = FileInputAXSYM, status = "old",                  &
        position = "rewind")   
  wavelength = 0.1_O * 2._O * Pi  
  ind_refMed = 1._O                                                                 
  ind_refRel = (1.5_O,0._O) 
  string     = 'OptProp'    
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if      
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if
  call check_ind_ref (ind_refRel)
  k = 2._O * Pi * ind_refMed / wavelength 
!
  perfectcond = .false.  
  chiral = .false.
  kb     = 0._O  
  string = 'MatProp'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) perfectcond
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable perfectcond;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) chiral
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable chiral;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) kb
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable kb;')"
      stop
    end if      
  else
    print "(/,2x,'Group name MatProp not found;')"
    stop  
  end if      
  call check_MatPropAXSYM (perfectcond, chiral, kb)   
  if (chiral) call check_chirality (kb)
!
  FileGeom = .false.
  FileFEM  = ' '  
  TypeGeom = 1  
  Nsurf    = 2
  do i = 1, NsurfPD
    surf(i)= 1._O
  end do
  Nparam = 1
  anorm  = 1._O
  Rcirc  = 1._O 
  miror  = .false.
  string = 'GeomProp'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) FileGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileGeom;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) FileFEM
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileFEM;')"
      stop
    end if       
    read (iInputAXSYM, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if       
    read (iInputAXSYM, *, iostat = ios) Nsurf
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nsurf;')"
      stop
    end if
    if (Nsurf > NsurfPD) then
      print "(/,2x,'Input error: Nsurf exceeds NsurfPD;')"                                    
      stop
    end if            
    do i = 1, Nsurf
      read (iInputAXSYM, *, iostat = ios) surf(i)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable surf;')"
        stop
      end if
    end do    
    read (iInputAXSYM, *, iostat = ios) Nparam
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nparam;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) miror
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable miror;')"
      stop
    end if                                                
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if                          
  call check_geomAXSYM (TypeGeom, Nsurf, Nparam)
  call check_geomAXSYMOblate (TypeGeom, Nsurf, surf)  
  call check_anorm (anorm)
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
  Nface = 1
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
  DoConvTest   = .true.
  MishConvTest = .true.
  string       = 'ConvTest'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) MishConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable MishConvTest;')"
      stop
    end if         
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if   
  if (chiral)   MishConvTest = .false.     
  if (FileGeom) MishConvTest = .false.  
!
  DS = .false.
  autGenDS = .true.
  string   = 'Sources'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) DS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DS;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) autGenDS
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable autGenDS;')"
      stop
    end if         
  else
    print "(/,2x,'Group name Sources not found;')"
    stop  
  end if          
  call check_inputAXSYM (miror, chiral, DS)
  if (FileGeom) autGenDS = .false.  
  InputDS = .false.
  if (DS .and. .not. autGenDS) InputDS = .true. 
!
  ComplexPlane = .false.
  epsZReIm = 0.95_O 
  if (DS .and. autGenDS) then
    string = 'SourcePosAut'
    if (XFindPar (iInputAXSYM, string)) then
      read (iInputAXSYM, *, iostat = ios) ComplexPlane
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ComplexPlane;')"
        stop
      end if
      read (iInputAXSYM, *, iostat = ios) epsZReIm
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable epsZReIm;')"
        stop
      end if         
    else
      print "(/,2x,'Group name SourcePosAut not found;')"
      stop  
    end if              
  end if
!
  Nint   = 100
  Nrank  =  17 
  if (.not. DoConvTest .or. InputDS) then       
    string = 'NintNrank'
    if (XFindPar (iInputAXSYM, string)) then
      read (iInputAXSYM, *, iostat = ios) Nint
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nint;')"
        stop
      end if
      read (iInputAXSYM, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if         
    else
      print "(/,2x,'Group name NintNrank not found;')"
      stop  
    end if                              
  end if    
!
  do i = 1, NrankPD
    zRe(i)  = 0._O
    zIm(i)  = 0._O
    zRe1(i) = 0._O
    zIm1(i) = 0._O
  end do
  if (InputDS) then   
    call check_MaxNrank (Nrank) 
    string = 'SourcePosInp'
    if (XFindPar (iInputAXSYM, string)) then
      do i = 1, Nrank
        read (iInputAXSYM, *, iostat = ios) zRe(i), zIm(i)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variables zRe and zIm;')"
          stop
        end if
      end do      
      read (iInputAXSYM, *)
      do i = 1, Nrank - 1
        read (iInputAXSYM, *, iostat = ios) zRe1(i), zIm1(i)
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variables zRe1 and zIm1;')"
          stop
        end if
      end do
    else
      print "(/,2x,'Group name SourcePosInp not found;')"
      stop  
    end if                       
  end if         
!
  epsNint  = 5.e-2_O
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  dNint    = 4
  dNintMrank = 10
  string   = 'Errors'
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if 
    read (iInputAXSYM, *, iostat = ios) dNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint;')"
      stop
    end if
    read (iInputAXSYM, *, iostat = ios) dNintMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNintMrank;')"
      stop
    end if         
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if
!
  FileTmat = '../TMATFILES/T.dat'
  string   = 'Tmat' 
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) FileTmat
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileTmat;')"
      stop
    end if             
  else
    print "(/,2x,'Group name Tmat not found;')"
    stop  
  end if   
!
  PrnProgress = .true.
  string   = 'PrintProgress' 
  if (XFindPar (iInputAXSYM, string)) then
    read (iInputAXSYM, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if       
  close (unit = iInputAXSYM) 
end subroutine readinputAXSYM   
!************************************************************************************
subroutine printinputAXSYM (ic, FileGeom, TypeGeom, FileFEM, Nsurf, Nparam, Nrank,  &
           Nrank1, dNint, dNintMrank,ind_refMed, wavelength, anorm, Rcirc, surf,    &
           kb, epsNint, epsNrank, epsMrank, zRe, zIm, zRe1, zIm1, ind_refRel,       &
           miror, perfectcond, DS, chiral, autGenDS)
  use parameters
  implicit none
  integer       :: ic, TypeGeom, Nsurf, Nparam, Nrank, Nrank1, dNint, dNintMrank,   &
                   i, LenString                     
  real(O)       :: ind_refMed, wavelength, anorm, Rcirc, surf(Nsurf), kb, epsNint,  &
                   epsNrank, epsMrank, zRe(Nrank), zIm(Nrank), zRe1(Nrank1),        &
                   zIm1(Nrank1)
  complex(O)    :: ind_refRel
  character(80) :: FileFEM, FileFEMWrite
  logical       :: FileGeom, miror, perfectcond, DS, chiral, autGenDS  
!
  write (iOutput,"(/,2x,'Input Parameters of the Scattering Problem:',/)")
  write (iOutput,"(2x,'wavelength of the free space, wavelength = ',1pe13.4,';')")  &
         wavelength
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'refractive index of the ambient medium, ind_refMed = ', ind_refMed, ';'
  write (iOutput,"(2x, a, 1pe10.3, ',', 1pe10.3, a)")                               &
 'relative refractive index of the particle, ind_refRel = (', ind_refRel, ');' 
  write (iOutput,*)
  if (FileGeom) then
    FileFEMWrite = FileFEM(14:LenString(FileFEM))
    write (iOutput,"(2x, a, a30)")                                                  &
   'name of the file containing the particle geometry, FileFEM = ', FileFEMWrite
  else
    write (iOutput,"(2x,'index of the particle geometry, TypeGeom = ',i2,';')")     &
           TypeGeom
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'spheroid;')")        
    else if (TypeGeom == 2) then
      write (iOutput,"(2x,'cylinder;')")        
    else if (TypeGeom == 3) then
      write (iOutput,"(2x,'rounded oblate cylinder;')")
    end if
    write (iOutput,"(2x,'number of surface parameters, Nsurf = ',i2,';')") Nsurf
    write (iOutput,"(2x,'surface parameters:')")
    do i = 1, Nsurf
      write (iOutput,"(2x,'surf(',i2,') = ',1pe10.3,',')") i, surf(i)
    end do
    if (TypeGeom == 1) then
      write (iOutput,"(2x,'where surf(1) is the semi-axis along the symmetry axis,')")
      write (iOutput,"(2x,'and   surf(2) is the second semi-axis;')")
    else if (TypeGeom == 2 .or. TypeGeom == 3) then
      write (iOutput,"(2x,'where surf(1) is the half-length of the cylinder')") 
      write (iOutput,"(2x,'and   surf(2) is the cylinder radius;')")
    end if               
    write (iOutput,"(2x,'number of integration surfaces, Nparam = ',i2,';')") Nparam
  end if  
  write (iOutput,"(2x,a, 1pe10.3, a)")                                              &
 'characteristic length of the particle, anorm = ', anorm, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'radius of the circumscribing sphere, Rcirc = ', Rcirc, ';'
  if (miror) write (iOutput,"(2x,'mirror symmetric particle;')")  
  write (iOutput,*)
  if (perfectcond) then
    write (iOutput,"(2x,'perfectly conducting particle;')")
  else if (chiral) then
    write (iOutput,"(2x,'chiral particle;')")
    write (iOutput,"(2x,'characteristic of chirality, kb = ',1pe10.3,';')") kb  
  else 
    write (iOutput,"(2x,'dielectric particle;')")
  end if 
  write (iOutput,*)
  if (.not.DS) then
    write (iOutput,"(2x,'computation with localized sources;')")
  else
    write (iOutput,"(2x,'computation with distributed sources;')")
    if (autGenDS) then
      write (iOutput,"(2x,'the sources are generated automatically;')") 
    else
      write (iOutput,"(2x,'the sources are user-defined;')")
    end if                      
    write (iOutput,*)       
    write (iOutput,"(2x,'z - coordinates of the sources:')")
    write (iOutput,"(2x,'first configuration:')")
    do i = 1, Nrank
      if (i /= Nrank) then
        write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,',')")&
               i, zRe(i), i, zIm(i)
      else
        write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,';')")&
               i, zRe(i), i, zIm(i)
      end if
    end do
    if (ic == 2) then
      write (iOutput,*)
      write (iOutput,"(2x,'second configuration:')")
      do i = 1,Nrank1
        if (i /= Nrank1) then
          write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,',')")&
                 i, zRe1(i), i, zIm1(i)
        else
          write (iOutput,"(2x,'zRe(',i3,') = ',1pe10.3,', zIm(',i3,') = ',1pe10.3,';')")&
                 i, zRe1(i), i, zIm1(i)
        end if
      end do
    end if
  end if
  write (iOutput,*)
  write (iOutput,"(2x,'integration tolerance, epsNint = ',1pe10.3,';')") epsNint
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum expansion order tolerance, epsNrank = ', epsNrank, ';'
  write (iOutput,"(2x, a, 1pe10.3, a)")                                             &
 'maximum azimuthal order tolerance, epsMrank = ', epsMrank, ';'          
  write (iOutput,"(2x,'integration step, dNint = ',i3,'.')") dNint
  if (DS) then
    write (iOutput,"(2x, a, i3, a)")                                                &
   'integration step for Mrank calculation, dNintMrank = ',     dNintMrank, '.'
  end if
  write (iOutput,"(/)")                       
end subroutine printinputAXSYM
!***********************************************************************************
subroutine convergence_NintAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,     &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, dNint, miror,  &
           perfectcond, DS, chiral, kb, epsNint, PrnProgress)
  use parameters                                                                           
  implicit none 
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint, dNint
  real(O)    :: k, surf(Nsurf), snorm, epsNint, kb, zRe(Nrank), zIm(Nrank),         &
                rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, iNint, NthetaConv,       &
                iprog, Nprog
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O  
  Mstart = 1
  Mrank  = 1
  Nmax   = Nrank        
  m = 1            
  call write_TypeConvHead (1)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (.not. chiral) then
    Nprog = 7
    iprog = 3
  else
    Nprog = 13
    iprog =  6
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)
  do iNint = 1, 2         
    allocate (paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)       
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 2+iprog*(iNint-1), Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 3+iprog*(iNint-1), Nprog)
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+iprog*(iNint-1), Nprog)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
    if (.not. chiral) then
      call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)
    else 
      call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 5+iprog*(iNint-1), Nprog)
      call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 6+iprog*(iNint-1), Nprog)
      call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
      if (PrnProgress) call write_progress (.false., 7+iprog*(iNint-1), Nprog)

    end if
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,    &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)
    call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,         &
         alfap, k, snorm, Cext, Qext)
    call write_3ConvParam (Nint, Nrank, m)
    call write_DSCS(Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)  
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do  
  call write_NintConvRes (NthetaConv, Nteta, epsNint)                 
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv)
end subroutine convergence_NintAXSYM 
!***********************************************************************************
subroutine convergence_NrankAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,    &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,         &
           perfectcond, DS, chiral, kb, epsNrank, PrnProgress)                  
  use parameters
  implicit none  
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint
  real(O)    :: k, surf(Nsurf), snorm, epsNrank, kb, zRe(Nrank), zIm(Nrank),        &
                rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv, iprog,       &
                Nprog, j   
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:), ap(:,:), bp(:,:),   &
                            am(:,:), bm(:,:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1
  Nmax   = Nrank    
  m = 1 
  call write_TypeConvHead (2)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))     
  allocate (ap(2*Nrank,2*Nrank), bp(2*Nrank,2*Nrank))
  allocate (am(2*Nrank,2*Nrank), bm(2*Nrank,2*Nrank))  
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
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
    Nprog = 5
  else
    Nprog = 8
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)             
! --- Nrank configuration ---                    
  call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 2, Nprog)
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, ap, 2*Nrank, 2*Nrank)
  call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 3, Nprog)
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, bp, 2*Nrank, 2*Nrank)
  call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, Nprog)
  iprog = 4
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
  if (.not. chiral) then
    call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)
  else 
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 4, Nprog)
    call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, am, 2*Nrank, 2*Nrank)        
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, b, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 5, Nprog)
    call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, bm, 2*Nrank, 2*Nrank)    
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 6, Nprog)
    iprog = 6
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v) 
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  call write_3ConvParam (Nint, Nrank, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- (Nrank - 1) configuration ---
  call copy_matrix (2*Nmax, 2*Nmax, ap, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank)
  call matrix_Nrank_m_left (Nmax, a, Nrank, Nrank)
  call copy_matrix (2*Nmax, 2*Nmax, bp, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank)
  call matrix_Nrank_m (Nmax, b, Nrank, Nrank)
  call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
  if (PrnProgress) call write_progress (.false., iprog+1, Nprog)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,       &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
  if (.not. chiral) then
    call matrix_m_negativ (Nmax, Nmax, b, Nrank, Nrank)
  else 
    call copy_matrix (2*Nmax, 2*Nmax, am, 2*Nrank, 2*Nrank, a, 2*Nrank, 2*Nrank)
    call matrix_Nrank_m_left (Nmax, a, Nrank, Nrank)
    call copy_matrix (2*Nmax, 2*Nmax, bm, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank)        
    call matrix_Nrank_m (Nmax, b, Nrank, Nrank)
    call LU_SYSTEM (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., iprog+2, Nprog)
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_3ConvParam (Nint, Nrank - 1, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank) 
  deallocate (ap, bp, am, bm)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_NrankAXSYM
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
subroutine convergence_MrankAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,    &
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
end subroutine convergence_MrankAXSYM
!***********************************************************************************
subroutine convergence_NintDSAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,   &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, dNint, miror,  &
           perfectcond, DS, chiral, kb, epsNint, PrnProgress)
  use parameters                                                                           
  implicit none 
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nint, dNint
  real(O)    :: k, surf(Nsurf), snorm, epsNint, kb, zRe(Nrank), zIm(Nrank),         &
                rp(2,NfacePD), np(2,NfacePD), area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, iNint, NthetaConv,       &
                iprog, Nprog 
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1
  Nmax   = Nrank        
  m = 1           
  call write_TypeConvHead (1)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (.not. chiral) then
    Nprog = 9
    iprog = 4
  else
    Nprog = 17
    iprog =  8
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)
  do iNint = 1, 2         
    allocate(paramG(Nparam,Nint), weightsG(Nparam,Nint), Nintparam(Nparam))
    call interpolation_listAXSYM (TypeGeom, Nsurf, surf, Nint, Nparam,              &
         Nintparam, paramG, weightsG, miror)       
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 2+iprog*(iNint-1), Nprog)
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,& 
         zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,    &
         Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 3+iprog*(iNint-1), Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+iprog*(iNint-1), Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,    &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,&
         2*Nrank)    
    if (PrnProgress) call write_progress (.false., 5+iprog*(iNint-1), Nprog)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank,     &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
    if (.not. chiral) then
      call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
    else 
      call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 6+iprog*(iNint-1), Nprog)
      call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area,     &
           Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,       &
           weightsG, b, Nrank, Nrank)
      if (PrnProgress) call write_progress (.false., 7+iprog*(iNint-1), Nprog)
      call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank,     &
           2*Nmax)
      if (PrnProgress) call write_progress (.false., 8+iprog*(iNint-1), Nprog)
      call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,   &
           area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, &
           weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)    
      call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank,          &
           b, 2*Nrank, 2*Nrank)
      if (PrnProgress) call write_progress (.false., 9+iprog*(iNint-1), Nprog)
    end if
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,    &
         Nmax, c)
    call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
    call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
    call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
        .false.,.true., h, v)
    call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
    call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,         &
         alfap, k, snorm, Cext, Qext)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)
    call write_3ConvParam (Nint, Nrank, m)
    call write_DSCS (Nteta,.false., h, v)
    call write_Effic (Qscat, Qext)
    Nint = Nint + dNint
    deallocate (paramG, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)                 
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv)
end subroutine convergence_NintDSAXSYM
!***********************************************************************************
subroutine convergence_NrankDSAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,  &
           surf, rp, np, area, Nface, zRe, zIm, zRe1, zIm1, Nparam, Nrank, Nrank1,  &
           Nint, miror, perfectcond, DS, chiral, kb, epsNrank, PrnProgress)              
  use parameters
  implicit none  
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Nrank, Nrank1, Nint
  real(O)    :: k, surf(Nsurf), snorm, epsNrank, kb, zRe(Nrank), zIm(Nrank),        &
                zRe1(Nrank1), zIm1(Nrank1), rp(2,NfacePD), np(2,NfacePD),           &
                area(NfacePD)
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, DS, chiral, PrnProgress
!      
  integer    :: Nteta, Mstart, Mrank, Nmax, Nmaxmax, i, m, NthetaConv,              &
                iprog, Nprog, j
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG(:,:), weightsG(:,:), h(:), v(:), oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), b(:,:), c(:), c1(:), cc(:)
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = 0._O
  beta   = 0._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Mstart = 1
  Mrank  = 1
  Nmax   = Nrank  
  m = 1
  call write_TypeConvHead (2)
  allocate (a(2*Nrank,2*Nrank), b(2*Nrank,2*Nrank), c(2*Nrank), c1(2*Nrank))
  Nmaxmax = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  allocate (cc(2*Nmaxmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))  
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
    Nprog = 9
  else
    Nprog = 17
  end if
  if (PrnProgress) call write_progress (.true., 1, Nprog)
! --- Nrank configuration ---                    
  call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., 2, Nprog)
  call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,  &
       zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,      &
       Nrank, Nrank)                              
  if (PrnProgress) call write_progress (.false., 3, Nprog)
  call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, Nprog)
  call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe, zIm, m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,  &
       miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,  &
       2*Nrank)    
  if (PrnProgress) call write_progress (.false., 5, Nprog)
  iprog = 5
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank, Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank, Nmax, Nmaxmax)
  if (.not. chiral) then
    call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
  else 
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 6, Nprog)
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,&
         zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,   &
         Nrank, Nrank)
    if (PrnProgress) call write_progress (.false., 7, Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 8, Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe, zIm, -m, Nrank, Nmax, Nint, Nparam, Nintparam, paramG,   &
         weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
    call product_matrices (2*Nmax, 2*Nrank, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank,&
         2*Nrank) 
    if (PrnProgress) call write_progress (.false., 9, Nprog)
        iprog = 9
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)  
  call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,   &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,           &
       alfap, k, snorm, Cext, Qext)
  do i = 1, Nteta
    oldh(i) = h(i)
    oldv(i) = v(i)
  end do  
  call write_3ConvParam (Nint, Nrank, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
! --- Nrank1 configuration ---
  Nmax   = Nrank1
  Nmaxmax= Nrank1 + Mrank*(2*Nrank1 - Mrank + 1)  
  call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe1, zIm1, m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG,         &
       weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., iprog+1, Nprog)
  call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,  &
       zRe1, zIm1, m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG, weightsG, b,   &
       Nrank, Nrank)
  if (PrnProgress) call write_progress (.false., iprog+2, Nprog)
  call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank1, 2*Nmax)
  if (PrnProgress) call write_progress (.false., iprog+3, Nprog)
  call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
       Nface, zRe1, zIm1, m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG,         &
       weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)
  call product_matrices (2*Nmax, 2*Nrank1, 2*Nmax, a, 2*Nrank, 2*Nrank, b, 2*Nrank, &
       2*Nrank)    
  if (PrnProgress) call write_progress (.false., iprog+4, Nprog)
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m, Nrank1,      &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_positive (c1, cc, m, Mstart, Nrank1, Nmax, Nmaxmax)  
  if (.not. chiral) then 
    call matrix_m_negativ (Nmax, Nmax, a, Nrank, Nrank)
  else 
    call matrix_Q_m (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe1, zIm1, -m, Nrank1, Nmax, Nint, Nparam, Nintparam,        &
         paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)    
    if (PrnProgress) call write_progress (.false., iprog+5, Nprog)
    call incident_matrix_m (FileGeom, TypeGeom, k, Nsurf, surf, rp, np, area, Nface,& 
         zRe1, zIm1, -m, Nrank1, Nmax, Nint, Nparam, Nintparam, paramG, weightsG,   &
         b, Nrank, Nrank)    
    if (PrnProgress) call write_progress (.false., iprog+6, Nprog)
    call LU_SYSTEM_DIRECT (a, 2*Nrank, 2*Nrank, b, 2*Nrank, 2*Nrank, 2*Nrank1,      &
         2*Nmax)    
    if (PrnProgress) call write_progress (.false., iprog+7, Nprog)
    call matrix_Q_m (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np,     &
         area, Nface, zRe1, zIm1, -m, Nrank1, Nmax, Nint, Nparam, Nintparam,        &
         paramG, weightsG, miror, perfectcond, DS, chiral, kb, a, Nrank, Nrank)    
    call product_matrices (2*Nmax, 2*Nrank1, 2*Nmax, a, 2*Nrank, 2*Nrank,           &
         b, 2*Nrank, 2*Nrank)       
    if (PrnProgress) call write_progress (.false., iprog+8, Nprog)
  end if
  call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank1,     &
       Nmax, c)
  call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
  call extend_vector_negative (c1, cc, m, Nrank1, Nmax, Nmaxmax)  
  call DSCS (cc, Mrank, Nrank1, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm,  &
      .false.,.true., h, v)
  call CQscat (cc, Mrank, Nrank1, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, Nrank1, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,          &
       alfap, k, snorm, Cext, Qext)
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConv)  
  call write_3ConvParam (Nint, Nrank1, m)
  call write_DSCS (Nteta,.false., h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConv, Nteta, epsNrank)
  deallocate (a, b, c, c1, cc, h, v, oldh, oldv, paramG, weightsG, Nintparam)
end subroutine convergence_NrankDSAXSYM 
!***********************************************************************************
subroutine convergence_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,  &
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
!  open (unit = iTmat, file = FileTmat, status = 'replace') 
!  call write_HeadFileTmat (Nrank, Nrank) 
!  call write_TypeConvHead (3)
!  call write_2ConvParamAxsym (Nint, Nrank)
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
!    call write_FileTmat (Nrank, Nrank, a)
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
!        call write_FileTmat (Nrank, Nrank, a)  
        call appendtotmat (a,btot,-m,Nrank,Nmax,Nmaxmax)
      end if
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m, Nrank,  &
           Nmax, c)
      call product_matrix_vector (2*Nmax, 2*Nmax, a, 2*Nrank, 2*Nrank, c, c1)
      call extend_vector_negative (c1, cc, m, Nrank, Nmax, Nmaxmax)
    end if
 !   call DSCS (cc, Mrank, Nrank, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k, snorm, &
 !       .false.,.true., h, v)
 !   call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)     
 !   call write_DSCS (Nteta,.false., h, v)
 !   if (NthetaConv >= int(0.8*Nteta)) exit 
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
end subroutine convergence_MrankDSAXSYM
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
    if (DoConvTest) then                      
      if (MishConvTest) then                         
        print "(/,2x, a)",                                                          &
       'Estimates of Nint and Nrank Using Mishchenko''s Convergence Procedure' 
        print "(  2x, a)",                                                          &
       '---------------------------------------------------------------------'                                         
        call estimateNrankMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,         &
             zRe, zIm, Nparam, miror, perfectcond, DS, ComplexPlane,                &
             epsZReIm, x, delta, Ndgs, Nint, Nrank, NrankMax, Cscat1, Cext1)
        call estimateNintMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,          &
             zRe, zIm, Nparam, miror, perfectcond, DS, x, delta, Ndgs, Nint,        &
             dNint, Nrank, NrankMax, Cscat1, Cext1)
        print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
        print "(  2x,'---------------------------------------------')"	       
        print "(/,2x,'- enter the estimated values of Nint and Nrank;')"            
        call read_integer2 (Nint, Nrank)                                
      else 
        print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
        print "(  2x,'---------------------------------------------')"             
        print "(/,2x,'Nrank estimate:')"                                                            
        print "(  2x, a, i3, a)",                                                   &  
       'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
        if (.not. FileGeom) then
          print "(/,2x, a)",                                                        &
         '- enter the estimated values of Nint and Nrank, where Nint = Ndgs * Nrank'
          print "(  2x,'  and Ndgs = 3,4,...;')"
          call read_integer2 (Nint, Nrank)                                                                      
        else
          print "(/,2x,'- enter the estimated value of Nrank;')"              
          call read_integer (Nrank)     
        end if  
      end if
      if (.not. FileGeom) then        
        print "(/,2x, a)",                                                          &
       '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;' 
        call read_integerbound (TypeConvTest, 1, 3)
      else
        print "(/,2x,'- enter the type of convergence test: 2 - Nrank, 3 - Mrank;')"
        call read_integerbound (TypeConvTest, 2, 3)
      end if
    else
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
    end if
    Nrank1 = Nrank - 1 ! redundant            
  else
    if (autGenDS) then
      if (DoConvTest) then                        
        if (MishConvTest) then           
          print "(/,2x, a)",                                                        &
         'Estimates of Nint and Nrank Using Mishchenko''s Convergence Procedure'
          print "(  2x, a)",                                                        &
         '---------------------------------------------------------------------'	  	 	 	                    
          call estimateNrankMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,       &
               zRe, zIm, Nparam, miror, perfectcond, DS, ComplexPlane,              &
               epsZReIm, x, delta, Ndgs, Nint, Nrank, NrankMax, Cscat1, Cext1)
          call estimateNintMishchenko (TypeGeom, k, ind_refRel, Nsurf, surf,        &
               zRe, zIm, Nparam, miror, perfectcond, DS, x, delta, Ndgs, Nint,      &
               dNint, Nrank, NrankMax, Cscat1, Cext1)   
          print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
          print "(  2x,'---------------------------------------------')"	       
          print "(/,2x,'- enter the estimated values of Nint and Nrank;')"           
          call read_integer2 (Nint, Nrank)                               
        else
          print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
          print "(  2x,'---------------------------------------------')"	
          print "(/,2x,'Nrank estimate:')"                                                                
          print "(  2x, a, i3, a)",                                                 &
         'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
          print "(/,2x, a)",                                                        &
         '- enter the estimated values of Nint and Nrank, where Nint = Ndgs * Nrank'
          print "(  2x,'  and Ndgs = 5,6,...;')"
          call read_integer2 (Nint, Nrank)                                                                                 
        end if      
        call check_MaxNrank (Nrank)             
        call zDSAXSYM (TypeGeom, Nsurf, surf, Nrank, ComplexPlane, epsZReIm, zRe, zIm)     
        print "(/,2x, a)",                                                          &
       '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;'
        call read_integerbound (TypeConvTest, 1, 3)
        if (TypeConvTest == 2) then
          Nrank1 = Nrank - 1
          call zDSAXSYM (TypeGeom, Nsurf, surf, Nrank1, ComplexPlane, epsZReIm,     &
               zRe1, zIm1)         
        end if
      else
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
      end if
      Nrank1 = Nrank - 1 ! redundant  
    else 
      if (DoConvTest) then      
        print "(/,2x,'Convergence Test for an Axisymmetric Particle')"
        print "(  2x,'---------------------------------------------')"
        print "(/,2x,'Input values:')"
        print "(  2x,'the input value of Nrank is ', i4,', while')", Nrank
        print "(  2x, a, i3, a)",                                                   &
       'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';'
        if (.not. FileGeom) then                                            
          print "(/,2x, a)",                                                        &
         '- enter the estimated value of Nint, where Nint = Ndgs * Nrank'
          print "(  2x,'  and Ndgs = 5,6,...;')"
          call read_integer (Nint) 
          print "(/,2x, a)",                                                        &
         '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;'
          call read_integerbound (TypeConvTest, 1, 3)
        else
          print "(/,2x,'- enter the type of convergence test: 2 - Nrank, 3 - Mrank;')"
          call read_integerbound (TypeConvTest, 2, 3)
        end if
        if (TypeConvTest == 2) Nrank1 = Nrank - 1
      else
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
      end if
      Nrank1 = Nrank - 1 ! redundant          
    end if
  end if                                                                            
! -----------------------------------------------------------------------------------
!                               Convergence test                                    !
! -----------------------------------------------------------------------------------
!  open (unit = iOutput, file = FileOutput, status = "replace") 
!  call printinputAXSYM (TypeConvTest, FileGeom, TypeGeom, FileFEM, Nsurf, Nparam,   &
!       Nrank, Nrank1, dNint, dNintMrank, ind_refMed, wavelength, anorm, Rcirc,      &
!       surf, kb, epsNint, epsNrank, epsMrank, zRe, zIm, zRe1, zIm1, ind_refRel,     &
!       miror, perfectcond, DS, chiral, autGenDS)                     
  if (DoConvTest) then              
    if (TypeConvTest == 1) then
      if (.not. DS) then                  
        call convergence_NintAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,       &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             dNint, miror, perfectcond, DS, chiral, kb, epsNint, PrnProgress)
      else 
        call convergence_NintDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,     &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             dNint, miror, perfectcond, DS, chiral, kb, epsNint, PrnProgress)
      end if
    else if (TypeConvTest == 2) then
      if (.not. DS) then        
        call convergence_NrankAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,      &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             miror, perfectcond, DS, chiral, kb, epsNrank, PrnProgress)
      else      
        call convergence_NrankDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,    &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, zRe1, zIm1, Nparam, Nrank, &
             Nrank1, Nint, miror, perfectcond, DS, chiral, kb, epsNrank,            &
             PrnProgress)
      end if
    else 
      if (.not. DS) then    
        call convergence_MrankAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,    &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             miror, perfectcond, DS, chiral, kb, epsMrank, FileTmat, PrnProgress,   &
             Nmaxmax,b)
      else 
        call convergence_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,    &
             Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,       &
             miror, perfectcond, DS, chiral, kb, epsMrank, dNintMrank, FileTmat,    &
             PrnProgress,Nmaxmax,b)
      end if
    end if  
  else      
    if (.not. DS) then    
      call convergence_MrankAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm, Nsurf, &
           surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint,                &
             miror, perfectcond, DS, chiral, kb, epsMrank, FileTmat, PrnProgress,   &
             Nmaxmax,b)
    else          
      call convergence_MrankDSAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,      &
           Nsurf, surf, rp, np, area, Nface, zRe, zIm, Nparam, Nrank, Nint, miror,  &
           perfectcond, DS, chiral, kb, epsMrank, dNintMrank, FileTmat, PrnProgress,&
           Nmaxmax,b)
    end if  
  end if 
end subroutine TAXSYM
