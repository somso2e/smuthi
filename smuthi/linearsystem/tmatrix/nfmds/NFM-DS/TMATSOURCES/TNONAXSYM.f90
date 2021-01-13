!***********************************************************************************
subroutine readinputNONAXSYM (wavelength, ind_refMed, ind_refRel, ind_RefRelZ,      &
           alphaPR, betaPR, perfectcond, anisotropic, chiral, kb, FileGeom,         &
           TypeGeom, FileFEM, Nsurf, surf, Nparam, anorm, Rcirc, miror, Nazimutsym, &   
           DoConvTest, ExtThetaDom, Nbeta, Nint1, Nint2, Nrank, Mrank, epsNint,     &
           epsNrank, epsMrank, dNint1, dNint2, FileTmat, PrnProgress, k, gammaPR,   &
           snorm, rp, np, area, Nface, TypeConvTest )
  use parameters
  use derived_parameters
  implicit none 
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Nazimutsym, TypeConvTest, Mrank, &
                   Nrank, NrankW, Nbeta, Nint1, Nint2, dNint1, dNint2, i, j, ios,   &
                   icall                     
  real(O)       :: k, ind_refMed, wavelength, anorm, surf(NsurfPD), xpart, snorm,   &
                   kb, epsNint, epsNrank, epsMrank, Rcirc, x, rp(3,NfacePD),        &
                   np(3,NfacePD), area(NfacePD), dp, alphaPR, betaPR, gammaPR, grd
  complex(O)    :: ind_refRel, ind_RefRelZ
  logical       :: FileGeom, miror, perfectcond, chiral, anisotropic, DoConvTest,   &
                   ExtThetaDom, PrnProgress, more, continuare, IntTest, XFindPar
  character(80) :: FileTmat, FileFEM, string    
! -----------------------------------------------------------------------------------
!                        Read the input file FileInputNONAXSYM                      ! 
! ----------------------------------------------------------------------------------- 
  call DrvParameters   
  open (unit = iInputNONAXSYM, file = FileInputNONAXSYM, status = "old",            &
        position = "rewind")   
  wavelength  = 0.1_O * 2._O * Pi
  ind_refMed  = 1._O
  ind_refRel  = (1.5_O,0._O)  
  string  = 'OptProp'
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) wavelength
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable wavelength;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) ind_refMed
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refMed;')"
      stop
    end if    
    read (iInputNONAXSYM, *, iostat = ios) ind_refRel
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ind_refRel;')"
      stop
    end if        
  else
    print "(/,2x,'Group name OptProp not found;')"
    stop  
  end if                          
  k   = 2._O * Pi * ind_refMed / wavelength
  grd = Pi / 180._O 
!
  perfectcond = .false.
  anisotropic = .false.  
  chiral      = .false.
  kb     = 0._O
  string = 'MatProp'    
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) perfectcond
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable perfectcond;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) anisotropic
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anisotropic;')"
      stop
    end if    
    read (iInputNONAXSYM, *, iostat = ios) chiral
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable chiral;')"
      stop
    end if    
    read (iInputNONAXSYM, *, iostat = ios) kb
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable kb;')"
      stop
    end if              
  else
    print "(/,2x,'Group name MatProp not found;')"
    stop  
  end if   
  call check_MatPropNONAXSYM (perfectcond, anisotropic, chiral, kb) 
  if (chiral) call check_chirality (kb) 
!
  FileGeom = .false.  
  FileFEM  = ' '
  TypeGeom = 1  
  Nsurf = 3
  do i = 1, NsurfPD
    surf(i) = 1._O
  end do
  Nparam = 1
  anorm  = 1._O
  Rcirc  = 1._O
  miror  = .true.
  Nazimutsym = 0
  string = 'GeomProp'    
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) FileGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileGeom;')"
      stop
    end if      
    read (iInputNONAXSYM, *, iostat = ios) FileFEM
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable FileFEM;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) TypeGeom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable TypeGeom;')"
      stop
    end if          
    read (iInputNONAXSYM, *, iostat = ios) Nsurf
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nsurf;')"
      stop
    end if
    if (Nsurf > NsurfPD) then
      print "(/,2x,'Input error: Nsurf exceeds NsurfPD;')"                                    
      stop
    end if
    do i = 1, Nsurf
      read (iInputNONAXSYM, *, iostat = ios) surf(i)
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable surf;')"
        stop
      end if
    end do 
    read (iInputNONAXSYM, *, iostat = ios) Nparam
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nparam;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) anorm
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable anorm;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) Rcirc
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Rcirc;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) miror
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable miror;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) Nazimutsym
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable Nazimutsym;')"
      stop
    end if
  else
    print "(/,2x,'Group name GeomProp not found;')"
    stop  
  end if
  call check_inputNONAXSYM (FileGeom, FileFEM, miror, chiral, anisotropic,          &
       Nazimutsym)
  call check_geom3D (TypeGeom, Nsurf, Nparam, anisotropic, miror, Nazimutsym)    
  call check_anorm (anorm)
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
  Nface = 1
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
!
  DoConvTest   = .true.
  ExtThetaDom  = .true.
  string = 'ConvTest'    
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) DoConvTest
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable DoConvTest;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) ExtThetaDom
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable ExtThetaDom;')"
      stop
    end if                     
  else
    print "(/,2x,'Group name ConvTest not found;')"
    stop  
  end if     
!
  ind_refRelZ = (1.5_O,0._O)
  alphaPR = 0._O
  betaPR  = 0._O
  Nbeta   = 60
  if (anisotropic) then
    string = 'AnSVWF'    
    if (XFindPar (iInputNONAXSYM, string)) then
      read (iInputNONAXSYM, *, iostat = ios) ind_refRelZ
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable ind_refRelZ;')"
        stop
      end if    
      read (iInputNONAXSYM, *, iostat = ios) alphaPR
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable alphaPR;')"
        stop
      end if    
      read (iInputNONAXSYM, *, iostat = ios) betaPR
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable betaPR;')"
        stop
      end if            
      read (iInputNONAXSYM, *, iostat = ios) Nbeta
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nbeta;')"
        stop
      end if
    else
      print "(/,2x,'Group name AnSVWF not found;')"
      stop  
    end if      
  end if
  if (anisotropic) call check_ind_ref3 (ind_refRel, ind_refRelZ)  
  alphaPR = alphaPR * grd
  betaPR  = betaPR  * grd
  gammaPR = 0._O  
!  
  if (DoConvTest) then
    if (.not. FileGeom) then
      print "(/,2x,'Convergence Test for a Nonaxisymmetric Particle')"
      print "(  2x,'-----------------------------------------------')"      
    else
      print "(/,2x, a)",                                                            &
     'Convergence Test for a Nonaxisymmetric Particle over Nrank and Mrank'
      print "(  2x, a)",                                                            &
     '--------------------------------------------------------------------'             
    end if      
  else
    print "(/,2x,'T-Matrix Computation for a Nonaxisymmetric Particle')"
    print "(  2x,'---------------------------------------------------')"
  end if
!  
  x = k * Rcirc
  NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
  Nrank  = 16
  Mrank  = 14
  if (.not. DoConvTest) then
    string = 'NrankMrank'    
    if (XFindPar (iInputNONAXSYM, string)) then
      read (iInputNONAXSYM, *, iostat = ios) Nrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Nrank;')"
        stop
      end if
      read (iInputNONAXSYM, *, iostat = ios) Mrank
      if (ios /= 0) then
        print "(/,2x,'Error by reading the input variable Mrank;')"
        stop
      end if
    else
      print "(/,2x,'Group name NrankMrank not found;')"
      stop  
    end if
    print "(/,2x,'Nrank and Mrank input values:')"
    print "(  2x,'the input values of Nrank and Mrank are ', i3,' and ', i3, ',')",  &
              Nrank, Mrank
    print "(  2x, a, i3, a)",                                                        &
   'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
  else  
    print "(/,2x,'Nrank estimate:')"  
    print "(  2x, a, i3, a)" ,                                                       &
   'the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
    print "(/,2x,'- enter the estimated values of Nrank and Mrank,')"
    print "(  2x, a)",                                                               &
   '  where for almost spherical particles Mrank = Nrank - 2,...,Nrank,'
    print "(  2x, a)",                                                               &
   '  while for less spherical particles Mrank can be smaller than Nrank - 2;'
    call read_integer2 (Nrank, Mrank)                                                               
  end if 
  call check_MrankNrank (Mrank, Nrank)  
! 
  Nint1 = 12
  Nint2 = 12  
  if (.not. FileGeom) then
    if (DoConvTest) then
      more  = .true.
      icall = 0
      do while (more)
        icall = icall + 1     
        print "(/,2x,'- enter the estimated values of Nint1 and Nint2;')"
        call read_integer2 (Nint1, Nint2)
        if (Nint1 < 5) Nint1 = 5
        if (Nint2 < 5) Nint2 = 5
        if (.not. anisotropic) then          
          call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, miror,            &
               Nazimutsym, Nint1, Nint2)
          if (icall == 1) then
            print "(/,2x, a)",                                                      &
           '- enter true to perform an integration test over specific matrix '	                             
            print "(  2x,'elements and false otherwise;')"
            call read_logical (IntTest)
          end if           
          if (IntTest) then            
            call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, Mrank,       &
                 Mrank, Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb,            &
                 ind_refRel, miror, perfectcond, chiral)
            call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, 1, Mrank,    &
                 Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb, ind_refRel,       &
                 miror, perfectcond, chiral)
          end if
        else
          call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, .false., 0,       &
          Nint1, Nint2)
          if (icall == 1) then
            print "(/,2x, a)",                                                      &
           '- enter true to perform an integration test over specific matrix '	                      
            print "(  2x,'elements and false otherwise;')"
            call read_logical (IntTest)
          end if           
          if (IntTest) then                  
            call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,          &
                 ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, Mrank, Mrank,  &
                 Nrank, 5, Nint1, Nint2, Nbeta, Nparam)
            call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,          &
                 ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, 1, Mrank,      &
                 Nrank, 5, Nint1, Nint2, Nbeta, Nparam)     
          end if
        end if 
        print "(2x,'- enter true for a new input of number of integration points')"
        print "(2x,'or false to continue;')"        
        call read_logical (continuare)
        if (.not. continuare) more = .false.
      end do
    else
      string = 'Nint'    
      if (XFindPar (iInputNONAXSYM, string)) then
        read (iInputNONAXSYM, *, iostat = ios) Nint1
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nint1;')"
          stop
        end if
        read (iInputNONAXSYM, *, iostat = ios) Nint2
        if (ios /= 0) then
          print "(/,2x,'Error by reading the input variable Nint2;')"
          stop
        end if
      else
        print "(/,2x,'Group name Nint not found;')"
        stop  
      end if
      if (Nint1 < 5) Nint1 = 5
      if (Nint2 < 5) Nint2 = 5
      print "(/,2x,'Nint1 and Nint2 input values:')"
      print "(  2x,'the input values of Nint1 and Nint2 are ', i4,' and',i4,';')",  &
                Nint1, Nint2
      if (.not. anisotropic) then
        call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, miror, Nazimutsym,  &
             Nint1, Nint2)        
        print "(/,2x,'- enter true to perform an integration test over')"                  
        print "(  2x,'specific matrix elements and false otherwise;')"
        call read_logical (IntTest)                        
        if (IntTest) then       
          call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, Mrank, Mrank,  &
               Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb, ind_refRel, miror,  &
               perfectcond, chiral )
          call ConvergenceMatrixElem (3, 1, TypeGeom, Nsurf, Nparam, 1, Mrank,      &
               Nrank, 5, Nint1, Nint2, Nazimutsym, k, surf, kb, ind_refRel, miror,  &
               perfectcond, chiral )                                                         
        end if
      else
        call SurfaceElemNint (k, TypeGeom, Nparam, Nsurf, surf, .false., 0,         &
             Nint1, Nint2)
        print "(/,2x,'- enter true to perform an integration test over')"                  
        print "(  2x,'specific matrix elements and false otherwise;')"
        call read_logical (IntTest)                        
        if (IntTest) then                          
          call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,            &
               ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, Mrank, Mrank,    &
               Nrank, 5, Nint1, Nint2, Nbeta, Nparam)
          call ConvergenceMatrixElemAnis (TypeGeom, 3, 1, k, ind_refRel,            &
               ind_refRelZ, alphaPR, betaPR, gammaPR, Nsurf, surf, 1, Mrank,        &
               Nrank, 5, Nint1, Nint2, Nbeta, Nparam)       
        end if
      end if
    end if
  end if
!  
  if (DoConvTest) then   
    if (.not. FileGeom) then   
      print "(/,2x, a)",                                                            &
     '- enter the type of convergence test: 1 - Nint, 2 - Nrank and Mrank;'
      call read_integerbound (TypeConvTest, 1, 2)              
    else    
      TypeConvTest = 2 
    end if  
  else
    TypeConvTest = 0
  end if
!         
  epsNint  = 5.e-2_O
  epsNrank = 5.e-2_O
  epsMrank = 5.e-2_O
  dNint1 = 4
  dNint2 = 4
  string   = 'Errors'
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) epsNint
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNint;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) epsNrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsNrank;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) epsMrank
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable epsMrank;')"
      stop
    end if 
    read (iInputNONAXSYM, *, iostat = ios) dNint1
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint1;')"
      stop
    end if
    read (iInputNONAXSYM, *, iostat = ios) dNint2
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable dNint2;')"
      stop
    end if         
  else
    print "(/,2x,'Group name Errors not found;')"
    stop  
  end if  
!
  FileTmat = '../TMATFILES/T.dat'
  string   = 'Tmat' 
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) FileTmat
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
  if (XFindPar (iInputNONAXSYM, string)) then
    read (iInputNONAXSYM, *, iostat = ios) PrnProgress
    if (ios /= 0) then
      print "(/,2x,'Error by reading the input variable PrnProgress;')"
      stop
    end if             
  else
    print "(/,2x,'Group name PrintProgress not found;')"
    stop  
  end if     
  close (unit = iInputNONAXSYM)  
end subroutine readinputNONAXSYM                            
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
subroutine convergence_NintNONAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm, Nsurf,  &
           surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, dNint1,   &
           dNint2, miror, Nazimutsym, perfectcond, chiral, kb, epsNint, ExtThetaDom,&
           PrnProgress)
  use parameters
  implicit none
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, dNint1, &
                dNint2, Nazimutsym                  
  real(O)    :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD), &
                kb, epsNint
  complex(O) :: ind_ref
  logical    :: FileGeom, miror, perfectcond, chiral, ExtThetaDom, PrnProgress
!      
  integer    :: Nmax, Nteta, i, NthetaConv, iNint, NintAL   
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1) 
  call write_TypeConvHead (1)
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 7)
  do iNint = 1, 2
    NintAL = max(Nint1,Nint2)                            
    allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),         &
              weightsG(Nparam,NintAL*NintAL))
    allocate (Nintparam(Nparam))
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
         Nintparam, paramG1, paramG2, weightsG, miror, Nazimutsym) 
    call matrix_Q (FileGeom, TypeGeom, 3, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
         Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,    &
         weightsG, miror, Nazimutsym, perfectcond, chiral, kb, a, Nmax, Nmax)  
    if (PrnProgress) call write_progress (.false., 2+3*(iNint-1), 7)
    call matrix_Q (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area, &
         Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,    &
         weightsG, miror, Nazimutsym, perfectcond, chiral, kb, b, Nmax, Nmax)
    if (PrnProgress) call write_progress (.false., 3+3*(iNint-1), 7)
    call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+3*(iNint-1), 7)
    call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,   &
         Nmax, c)    
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
    call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,    &
         ExtThetaDom,.true., h, v)
    call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
    call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,            &
         alfap, k, snorm, Cext, Qext)     
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv)  
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)
    call write_DSCS (Nteta, ExtThetaDom, h, v)
    call write_Effic (Qscat, Qext)
    Nint1 = Nint1 + dNint1
    Nint2 = Nint2 + dNint2
    deallocate (paramG1, paramG2, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)
  deallocate (a, b, c, c1, h, v, oldh, oldv)    
end subroutine convergence_NintNONAXSYM 
!***********************************************************************************
subroutine convergence_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_ref, snorm,  &
           Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,    &
           miror, Nazimutsym, perfectcond, chiral, kb, epsNrank, epsMrank,          &
           ExtThetaDom, FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,      &
                   Nazimutsym                   
  real(O)       :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD),             &
                   area(NfacePD), kb, epsNrank, epsMrank
  complex(O)    :: ind_ref
  logical       :: FileGeom, miror, perfectcond, chiral, ExtThetaDom, PrnProgress 
  character(80) :: FileTmat
!      
  integer       :: Nmax, Nteta, i, j, NthetaConvN, NthetaConvM, NintAL
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  integer,allocatable    :: Nintparam(:) 
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:), a1(:,:), b1(:,:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)
  NintAL = max(Nint1,Nint2) 
  call write_TypeConvHead (4)  
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (Nmax, Nmax)     
  allocate (a1(2*Nmax,2*Nmax), b1(2*Nmax,2*Nmax))  
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))                        
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
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax)  
  call matrix_Q (FileGeom, TypeGeom, 1, 1, k, ind_ref, Nsurf, surf, rp, np, area,   &
       Nface, Mrank, Nrank, Nmax, NintAL, Nparam, Nintparam, paramG1, paramG2,      &
       weightsG, miror, Nazimutsym, perfectcond, chiral, kb, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)    
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax) 
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
  call write_FileTmat (Nmax, Nmax, b)
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)  
  if (.not. FileGeom) then             
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)                  
  else
    call write_2ConvParam (Nrank, Mrank)
  end if
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  deallocate (paramG1, paramG2, weightsG, Nintparam)
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)
  end do
  close (unit = iTmat)  
! --- (Nrank - 1) configuration ---
  if (PrnProgress) call write_progress_low  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax,2*Nmax)    
  call matrix_Nrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax) 
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a,2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)          
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  if (.not. FileGeom) then    
    call write_4ConvParam (Nint1, Nint2, Nrank - 1, Mrank)
  else    
    call write_2ConvParam (Nrank - 1, Mrank)
  end if
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)
! --- (Mrank - 1) configuration ---  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax)   
  call matrix_Mrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax) 
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)         
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  if (.not. FileGeom) then
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank - 1)
  else
    call write_2ConvParam (Nrank, Mrank - 1)
  end if
  call write_DSCS (Nteta, ExtThetaDom, h, v) 
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"                                
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., chiral)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., chiral)
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank    
  deallocate (a1, b1)  
  deallocate (a, b, c, c1, h, v, oldh, oldv, oldh0, oldv0)      
end subroutine convergence_Nrank_MrankNONAXSYM 
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
!                   ROUTINES FOR UNIAXIAL ANISOTROPIC PARTICLES                    *
!***********************************************************************************
subroutine convergence_NintAnis (FileGeom, TypeGeom, k, ind_ref, ind_refZ, alfaPR,  &
           betaPR, gamaPR, snorm, Nsurf, surf, rp, np, area, Nface, Nparam, Mrank,  &
           Nrank, Nbeta, Nint1, Nint2, dNint1, dNint2, epsNint, ExtThetaDom,        &
           PrnProgress)
  use parameters
  implicit none
  integer    :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nbeta, Nint1, Nint2,  &
                dNint1, dNint2
  real(O)    :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD), area(NfacePD), &
                epsNint, alfaPR, betaPR, gamaPR
  complex(O) :: ind_ref, ind_refZ  
  logical    :: FileGeom, ExtThetaDom, PrnProgress
!      
  integer    :: Nmax, Nteta, i, NthetaConv, iNint, NintAL  
  real(O)    :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,        &
                Cext, Qext
  integer,allocatable    :: Nintparam(:) 
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)         
  call write_TypeConvHead (1)
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do
  if (PrnProgress) call write_progress (.true., 1, 7)
  do iNint = 1, 2                       
    NintAL = max(Nint1,Nint2) 
    allocate (paramG1(Nparam,NintAL*NintAL), paramG2(Nparam,NintAL*NintAL),         &
              weightsG(Nparam,NintAL*NintAL))
    allocate (Nintparam(Nparam))
    call interpolation_list3D (TypeGeom, Nsurf, surf, Nint1, Nint2, NintAL, Nparam, &
         Nintparam, paramG1, paramG2, weightsG,.false., 1)   
    call matrix_Q_anis (FileGeom, TypeGeom, 3, 1, k, ind_ref, ind_refZ, alfaPR,     &
         betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax,      &
         Nbeta, NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, a, Nmax, Nmax)
    if (PrnProgress) call write_progress (.false., 2+3*(iNint-1), 7)                                                    
    call matrix_Q_anis (FileGeom, TypeGeom, 1, 1, k, ind_ref, ind_refZ, alfaPR,     &
         betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax,      &
         Nbeta, NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, b, Nmax, Nmax)
    if (PrnProgress) call write_progress (.false., 3+3*(iNint-1), 7)
    call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
    if (PrnProgress) call write_progress (.false., 4+3*(iNint-1), 7)
    call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,   &
         Nmax, c)    
    call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
    call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,    &
         ExtThetaDom,.true., h, v)
    call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
    call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,            &
         alfap, k, snorm, Cext, Qext)     
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsNint, NthetaConv) 
    call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)                  
    call write_DSCS (Nteta, ExtThetaDom, h, v)
    call write_Effic (Qscat, Qext)
    Nint1 = Nint1 + dNint1
    Nint2 = Nint2 + dNint2
    deallocate (paramG1, paramG2, weightsG, Nintparam)
  end do
  call write_NintConvRes (NthetaConv, Nteta, epsNint)                           
  deallocate (a, b, c, c1, h, v, oldh, oldv)    
end subroutine convergence_NintAnis
!***********************************************************************************
subroutine convergence_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_ref, ind_refZ,   &
           alfaPR, betaPR, gamaPR, snorm, Nsurf, surf, rp, np, area, Nface, Nparam, &
           Mrank, Nrank, Nbeta, Nint1, Nint2, epsNrank, epsMrank, ExtThetaDom,      &
           FileTmat, PrnProgress)
  use parameters
  implicit none
  integer       :: TypeGeom, Nsurf, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, Nbeta
  real(O)       :: k, snorm, surf(Nsurf), rp(3,NfacePD), np(3,NfacePD),             &
                   area(NfacePD), epsNrank, epsMrank, alfaPR, betaPR, gamaPR
  complex(O)    :: ind_ref, ind_refZ
  character(80) :: FileTmat
  logical       :: FileGeom, ExtThetaDom, PrnProgress
!        
  integer       :: Nmax, Nteta, i, NthetaConvN, NthetaConvM, NintAL
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  integer,allocatable    :: Nintparam(:)
  real(O),allocatable    :: paramG1(:,:), paramG2(:,:), weightsG(:,:), h(:), v(:),  &
                            oldh(:), oldv(:), oldh0(:), oldv0(:)
  complex(O),allocatable :: a(:,:), c(:), b(:,:), c1(:), a1(:,:), b1(:,:)   
!
  tetaGI = 0._O
  phiGI  = 0._O
  phiGS  = 0._O
  Nteta  = 10
  alfa   = Pi / 4._O
  beta   = Pi / 4._O
  gama   = 0._O
  alfap  = Pi / 4._O
  Nmax   = Nrank + Mrank * (2 * Nrank - Mrank + 1)   
  NintAL = max(Nint1,Nint2)
  call write_TypeConvHead (4)  
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call  write_HeadFileTmat (Nmax, Nmax)    
  allocate (a1(2*Nmax,2*Nmax), b1(2*Nmax,2*Nmax))  
  allocate (a(2*Nmax,2*Nmax), c(2*Nmax), b(2*Nmax,2*Nmax), c1(2*Nmax))
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta), oldh0(Nteta),             &
            oldv0(Nteta))                        
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
  call copy_matrix (2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax)
  call matrix_Q_anis (FileGeom, TypeGeom, 1, 1, k, ind_ref, ind_refZ, alfaPR,       &
       betaPR, gamaPR, Nsurf, surf, rp, np, area, Nface, Mrank, Nrank, Nmax, Nbeta, &
       NintAL, Nparam, Nintparam, paramG1, paramG2, weightsG, b, Nmax, Nmax)
  if (PrnProgress) call write_progress (.false., 3, 4)  
  call copy_matrix (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)
  if (PrnProgress) call write_progress (.false., 4, 4)
  call write_FileTmat (Nmax, Nmax, b)
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)             
  call write_4ConvParam (Nint1, Nint2, Nrank, Mrank)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  deallocate (paramG1, paramG2, weightsG, Nintparam)
  do i = 1, Nteta
    oldh(i)  = h(i)
    oldv(i)  = v(i)
    oldh0(i) = h(i)
    oldv0(i) = v(i)
  end do
  close (unit = iTmat)  
! --- (Nrank - 1) configuration ---
  if (PrnProgress) call write_progress_low  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax) 
  call matrix_Nrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1,  2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax)  
  call matrix_Nrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)         
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh, oldv, epsNrank, NthetaConvN)
  call write_4ConvParam (Nint1, Nint2, Nrank - 1, Mrank)
  call write_DSCS (Nteta, ExtThetaDom, h, v)
  call write_Effic (Qscat, Qext)
  call write_NrankConvRes (NthetaConvN, Nteta, epsNrank)                                          
! --- (Mrank - 1) configuration ---  
  call copy_matrix (2*Nmax, 2*Nmax, a1, 2*Nmax, 2*Nmax, a, 2*Nmax, 2*Nmax) 
  call matrix_Mrank_1_left (Mrank, Nrank, Nmax, a, Nmax, Nmax)  
  call copy_matrix (2*Nmax, 2*Nmax, b1, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax)  
  call matrix_Mrank_1_right (Mrank, Nrank, Nmax, b, Nmax, Nmax)
  call LU_SYSTEM (a, 2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, 2*Nmax)         
  call PWcoefficients_ab (tetaGI, phiGI, alfa, beta, gama, alfap, Mrank, Nrank,     &
       Nmax, c)    
  call product_matrix_vector (2*Nmax, 2*Nmax, b, 2*Nmax, 2*Nmax, c, c1)
  call DSCS (c1, Mrank, Nrank, Nmax, Nteta, phiGS, alfa, beta, gama, k, snorm,      &
       ExtThetaDom,.true., h, v)
  call CQscat (c1, Mrank, Nrank, Nmax, k, snorm, Cscat, Qscat)
  call CQext (c1, Mrank, Nrank, Nmax, tetaGI, phiGI, alfa, beta, gama,              &
       alfap, k, snorm, Cext, Qext)    
  call delta_DSCS (Nteta, h, v, oldh0, oldv0, epsMrank, NthetaConvM)
  call write_4ConvParam (Nint1, Nint2, Nrank, Mrank - 1)
  call write_DSCS (Nteta, ExtThetaDom, h, v) 
  call write_Effic (Qscat, Qext)
  call write_MrankConvRes (NthetaConvM, epsMrank)
  if (NthetaConvN >= int(0.8*Nteta) .and. NthetaConvM >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criteria for Nrank and Mrank are satisfied;')"
  else
    print "(/,2x,'Convergence criteria for Nrank and Mrank are not satisfied;')"
  end if          
  deallocate (a1, b1)  
  call write_InfoFileTmat (FileTmat, Mrank, Nrank, .false., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, Nrank, .false., .false., .false.)
  print "(/,2x,'The T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,';')", Nrank
  print "(  2x,'- number of azimuthal modes, Mrank = ',i3,';')", Mrank  
  deallocate (a, b, c, c1, h, v, oldh, oldv, oldh0, oldv0)      
end subroutine convergence_Nrank_MrankAnis
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
  if (DoConvTest) then
    open (unit = iOutput, file = FileOutput, status = "replace")
    call printinputNONAXSYM (FileGeom, TypeGeom, FileFEM, Nsurf, Nparam, Nazimutsym,&
         dNint1, dNint2, ind_refMed, wavelength, anorm, Rcirc, surf, kb, epsNint,   &
         epsNrank, epsMrank,  ind_refRel, ind_refRelZ, alphaPR, betaPR, anisotropic,&
         miror, perfectcond, chiral) 
    if (TypeConvTest == 1) then  
      if (.not. anisotropic) then             
        call convergence_NintNONAXSYM (FileGeom, TypeGeom, k, ind_refRel, snorm,    &
             Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2,  &
             dNint1, dNint2, miror, Nazimutsym, perfectcond, chiral, kb, epsNint,   &
             ExtThetaDom, PrnProgress)      
      else
        call convergence_NintAnis (FileGeom, TypeGeom, k, ind_refRel, ind_refRelZ,  &
             alphaPR, betaPR, gammaPR, snorm, Nsurf, surf, rp, np, area, Nface,     &
             Nparam, Mrank, Nrank, Nbeta, Nint1, Nint2, dNint1, dNint2, epsNint,    &
             ExtThetaDom, PrnProgress) 
      end if
    else if (TypeConvTest == 2) then      
      if (.not. anisotropic) then
        call convergence_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_refRel,    &
             snorm, Nsurf, surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1,  &
             Nint2, miror, Nazimutsym, perfectcond, chiral, kb, epsNrank, epsMrank, &
             ExtThetaDom, FileTmat, PrnProgress)             
      else
        call convergence_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_refRel,        &
             ind_refRelZ, alphaPR, betaPR, gammaPR, snorm, Nsurf, surf, rp, np,     &
             area, Nface, Nparam, Mrank, Nrank, Nbeta, Nint1, Nint2, epsNrank,      &
             epsMrank, ExtThetaDom, FileTmat, PrnProgress)  
        end if
    end if
    close (unit = iOutput)           
  else       
    if (.not. anisotropic) then                
      call TMatrix_Nrank_MrankNONAXSYM (FileGeom, TypeGeom, k, ind_refRel, Nsurf, &
           surf, rp, np, area, Nface, Nparam, Mrank, Nrank, Nint1, Nint2, miror,    &
           Nazimutsym, perfectcond, chiral, kb, FileTmat, PrnProgress,Nmax,b)                                                     
    else
      call TMatrix_Nrank_MrankAnis (FileGeom, TypeGeom, k, ind_refRel, ind_refRelZ, &
           alphaPR, betaPR, gammaPR, Nsurf, surf, rp, np, area, Nface, Nparam,      &
           Mrank, Nrank, Nbeta, Nint1, Nint2, FileTmat, PrnProgress,Nmax,b) 
    end if 
  end if    
end subroutine TNONAXSYM
