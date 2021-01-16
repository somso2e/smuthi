subroutine TLAY(kuser,ind_ref,surf,Npart,Nsurfmax,Nrankpmax,btot)
!------------------------------------------------------------------------------------
! 1. General Considerations                                                         !
! --------------------------                                                        !
! TLAY is a routine for computing the T matrix and the scattering characteristics   !
! of axisymmetric layered particles. The particle consists of Npart layers, each    !
! characterized by arbitrary but constant values of electric permittivity and       !
! magnetic permeability. Each layer is bounded by two surfaces. The enclosing       !
! (exterior) surface will be referred to as the layer surface, and obviously, the   !
! first layer surface is the largest surface.                                       ! 
!                                                                                   !
! The domain of analysis is divided into homogeneous layers. The electric and       !
! magnetic fields of a specific layer i can be approximated by using localized      !
! or distributed vector spherical wave functions, and Nrankp(i) will be referred    !
! to as the maximum expansion order of the layer i.                                 !
!                                                                                   !
! The transition matrix of a layered particle can be obtained by considering the    !
! null-field condition for the total electric field inside and outside each region. !
! Note that for axisymmetric layered particles, the scattering problem decouples    !
! over the azimuthal modes m and the T matrix can be computed separately for each m.!
!                                                                                   !
! Particle geometries currently supported include layered spheroids and cylinders.  !
! The following routines (from the file "GeomLib.f90") provide the required         !
! geometry parameters:                                                              !
! - the routine "interpolation_listLAY" computes the integration nodes and weights  !
!   for each layer generatrix,                                                      !
! - the routine "elem_geomLAY" provides the geometry parameters (position vector,   !
!   polar and azimuthal angles and normal unit vector in spherical coordinates)     !
!   at each integration node, and                                                   !
! - the routine "zDSLAY" generates the coordinates of the distributed sources.      !
! The generatrix of each layer is described with respect to a local coordinate      !
! system and the axial positions of the local coordinate systems are specified in   !
! the global coordinate system of the layered particle by the array OriginPart.     !
! The global coordinate system can be chosen as the local coordinate system of the  !
! first layer surface (host particle) by setting OriginPart(1) = 0.0. By convention,! 
! the COORDINATES OF THE DISTRIBUTED SOURCES corresponding to a specific layer are  !
! provided in the GLOBAL COORDINATE SYSTEM. The user can easily modify the above    !
! routines to generate particles with other geometries. Note that the list of       !
! parameters must be maintained and only geometries with an analytical description  !
! of the surface can be implemented.                                                !
!                                                                                   !
! 5. Additional Comments                                                            !
! -----------------------                                                           !
! For the expansion order test and distributed sources, two configurations of       !
! sources are considered. For both configurations, the coordinates of the           !
! distributed sources can be generated automatically or can be specified in the     !
! input file.                                                                       !
!                                                                                   !
! The convergence tests over Nint and Nrank can be switched off by setting the      !
! logical variable DoConvTest to false. In this case, the values of Nint and        !
! Nrank must be specified in the input file.                                        !
!                                                                                   !
! Nint and Nrank have been determined for the azimuthal orders m = 1 and m = - 1.   ! 
! For higher azimuthal orders and sources distributed in the complex plane, the     !
! chosen values of Nint and Nrank may lead to false convergence. There are two      !
! strategies for improving the numerical stability.                                 !
! 1. The number of integration points can be increased at each azimuthal mode       !
! calculation with the integration step dNintMrank. The parameter dNintMrank         !
! should be found experimentally by repeating the azimuthal order test for          !
! different input values until the differential scattering cross section converges. !
! A conservative estimate is: dNintMrank = (0.1,...,0.2) * Nint.                    !
! 2. For sources generated automatically, an interactive convergence test over      !
! Nint and Nrank can be performed at each azimuthal mode calculation. The general   !
! strategy is to reduce the number of discrete sources (Nrankp(i), i = 1,2,...,     !
! Npart) and to increase the number of integration points (Nint) for increasing     !
! values of m. However, in the current version of the program, this test is not     !
! active (see the comments in the subroutine "tmatrix_MrankDSLAY").             !
!                                                                                   !
! The number of azimuthal modes is the same for all layers, and for localized       !
! sources, the relation Mrank <= Nrankp(i), for all i = 1,2,...,Npart, should be    !
! satisfied. However, the Nrankp(i), i = 1,2,...,Npart, are determined before Mrank !
! is computed. Therefore, for the azimuthal mode calculation m, with |m| > 1,       !
! we set Nrankp(i) = |m|, if Nrankp(i) < |m| for some i. In this case, the maximum  !
! expansion order of a layer is m-dependent and                                     !
!                                                                                   !
!                    n = |m|,|m|+1,...,Nrankp(i,|m|).                               !
!                                                                                   !
! 6. Input Parameters                                                               !
! ---------------------                                                             !
! The parameters specified in the input file "/INPUTFILES/InputLAY.dat" are listed  !
! below                                                                             !
!                                                                                   !
! - wavelength (real) - wavelength of the incident light in vacuo.                  !
!                                                                                   !
! - ind_refMed (real) - refractive index of the ambient medium (nonabsorbing        !
!   dielectric medium).                                                             !
!                                                                                   !
! - TypeGeom (integer) - parameter specifying the type of the particle geometry.    !
!                                                                                   !
! - Npart (integer) - number of layers.                                             !
!                                                                                   !
! - anorm (real) - characteristic length of the layered particle which is used      !
!   to normalize the differential scattering cross sections.                        !
!                                                                                   !
! - DoConvTest (logical) - if DoConvTest = t, the interactive convergence tests     !
!   over Nint and Nrank are invoked. Estimates of Nrank for all layers are          !
!   provided by Wiscombe's truncation limit criterion. If DoConvTest = f, the       !
!   values of Nint and Nrank must be supplied in the input file.                    !
!                                                                                   !
! - DS (logical) - if DS = t, distributed sources are used for T-matrix             !
!   calculation, otherwise localized sources are employed.                          !
!                                                                                   !
! - autGenDS (logical) - if autGenDS = t, the routine "zDSLAY" generates the        !
!   coordinates of the distributed sources. Otherwise, the coordinates of the       !
!   distributed sources must be provided by the user. If DS = f, the logical        ! 
!   variable autGenDS is ignored.                                                   ! 
!                                                                                   !
! - Nint (integer) - global number of integration points for the layered            !
!   particle. This parameter is used if the convergence tests are not performed     !
!   (DoConvTest = f).                                                               !
!                                                                                   !
! The next parameters (specified in a sequence of group statements) characterize     !
! each layer.                                                                       !
! - ind_refPartRel (complex) - relative refractive index of the actual layer        !
!   with respect to the ambient medium. The imaginary part of the relative          !
!   refractive index must be zero for nonabsorbing particles and positive for       !
!   absorbing particles.                                                            !
!                                                                                   !
! - NsurfPart (integer) - number of surface parameters of the actual layer.         !
!                                                                                   !
! - NparamPart (integer) - number of smooth curves forming the generatrix curve of  !
!   the actual layer.                                                               !
!                                                                                   !
! - OriginPart (real) - axial position of the local coordinate system of the        !
!   actual layer with respect to the global coordinate system of the layered        !
!   particle. This parameter can be positive or negative.                           !
!                                                                                   !
! - surfPart (real array: surfPart(1), surfPart(2),...,surfPart(NsurfPart)) -       !
!   surface parameters specifying the shape of the actual layer. The dimension of   !
!   the array is NsurfPD. The integer parameter NsurfPD is specified in the         !
!   routine  "Parameters.f90" and has the value NsurfPD = 10. If                    !
!   NsurfPart > NsurfPD, the execution is automatically terminated.                 !
!                                                                                   !
!   The permissive values of the surface parameters (for the supported geometries)  !
!   are summarized below.                                                           !
!                                                                                   !
!      Particle      TypeGeom   Nsurf   Nparam                surf                  ! 
!     concentric        1         2        1        surf(1) - length of the semi-   ! 
!     spheroids                                               axis along the        ! 
!                                                             symmetry axis         !  
!                                                   surf(2) - length of the second  ! 
!                                                             semi-axis             !     
!                                                                                   ! 
!     concentric        2         2        3        surf(1) - half-length of        ! 
!     cylinders                                               the cylinder          !  
!                                                   surf(2) - cylinder radius       ! 
!                                                                                   !
! - lnormPart (real) - characteristic length of the actual layer (usually the       !
!   radius of the smallest circumscribing sphere) which is used to compute an       !
!   estimate of the maximum expansion order by using Wiscombe's truncation limit    !
!   criterion. Alternatively, lnormPart can be chosen as the equal-volume sphere    !
!   radius or the surface-equivalent-sphere radius.                                 !
!                                                                                   !
! - ComplexPlanePart (logical) - if ComplexPlanePart = t, the distributed sources   !
!   are situated in the complex plane. This parameter is used if distributed sources!
!   are required and the coordinates of the distributed sources are automatically   !
!   generated (DS = t and autGenDS = t).                                            !
!                                                                                   !
! - EpsZReImPart (real) - input parameter of the routine "zDSLAY" which controls    !
!   the distribution of the discrete sources. This parameter is used if distributed !
!   sources are required and the coordinates of the distributed sources are         !
!   automatically generated (DS = t and autGenDS = t).                              !
!                                                                                   !
! - NrankPart (integer) - maximum expansion order for the actual layer. This        !
!   parameter is used if the convergence tests are not performed, or the coordinates!
!   of the distributed sources are user-defined. More specifically, NrankPart is    !
!   used if                                                                         !
!   - (DS = f, DoConvTest = f),                                                     !
!   - (DS = t, autGenDS = t, DoConvTest = f), or                                    ! 
!   - (DS = t, autGenDS = f, DoConvTest = t/f).                                     !
!                                                                                   !
! - zRePart, zImPart (real arrays: zXPart(1), zXPart(2),...,zXPart(Nrankp),         !
!   X = Re, Im) - coordinates of the distributed sources for the actual layer and   !
!   the expansion order NrankPart. These parameters are used if the coordinates of  !
!   the distributed sources are user-defined (DS = t and autGenDS = f). The         !
!   dimension of the arrays zRePart and zImPart is NrankPD and the inequality       !
!   NrankPart <= NrankPD must hold. The integer parameter NrankPD is specified in   !
!   the routine "Parameters.f90" and has the value NrankPD = 200. If                !
!   NrankPart > NrankPD, the execution is automatically terminated. Note that the   !
!   coordinates of the distributed sources are defined with respect to the GLOBAL   !
!   COORDINATE SYSTEM of the layered particle.                                      !
!                                                                                   !
! - zRePart1, zImPart1 (real arrays: zXPart1(1), zXPart1(2),...,zXPart1(Nrankp-1),  !
!   X = Re, Im) - coordinates of the distributed sources for the actual layer and   !
!   the expansion order NrankPart - 1. These parameters are used if the             !
!   coordinates of the distributed sources are user-defined and the expansion order !
!   test is performed (DS = t and autGenDS = f and TypeConvTest = 2). The dimension !
!   of the arrays zRePart1 and zImPart1 is NrankPD. As before, zRePart1 and zImPart1!
!   are defined with respect to the GLOBAL COORDINATE SYSTEM of the layered         !
!   particle.                                                                       !
!                                                                                   !
! NOTE: THE INPUT ARRAYS zRe, zIm AND zRe1, zIm1 MUST BE SPECIFIED IF (DS = t AND   !
! autGenDS = f), AND MUST BE SEPARATED BY A BLANK LINE. IF THE EXPANSION ORDER TEST !
! IS NOT PERFORMED (TypeConvTest /= 2), THE INPUT ARRAYS zRePart1 AND zImPart1 CAN  !
! BE OMITTED.                                                                       !
!                                                                                   !
! - epsNint (real) - error tolerance for the integration test.                      !
!                                                                                   !
! - epsNrank (real) - error tolerance for the expansion order test.                 !
!                                                                                   !
! - epsMrank (real) - error tolerance for the azimuthal order test.                 !
!                                                                                   !
! - dNint (integer) - number of division points for the integration test. Note      !
!   that the scattering problem is solved for Nint and Nint + dNint.                ! 
!                                                                                   !
! - dNintMrank (integer) - number of division points for azimuthal mode             !
!   calculation. This parameter is used if distributed sources are employed         !
!   (DS = t) and THERE EXISTS AT LEST ONE LAYER with ComplexPlanePart = t). At      !
!   each azimuthal mode calculation, the number of integration points increases     !
!   with dNintMrank. For layers with sources distributed along the symmetry axis,   !
!   the number of integration points is also increased.                             ! 
!                                                                                   !
! - FileTmat (character(80)) - name of the file to which the T matrix is written.   !
!                                                                                   !
! - PrnProgress (logical) - if PrnProgress = t, the progress of calculation         !
!   is printed.                                                                     !
!                                                                                   !
! Note that all input parameters specifying lengths must be provided in the same    !
! units of lengths.                                                                 !
!                                                                                   !
! The statements with the group names OptRegProp, GeomRegProp, SourceRegPosAut,     !
! NrankReg and SourceRegPosInp characterize a specific layer. THESE STATEMENTS      !
! MUST BE REPEATED FOR ALL Npart LAYERS.                                            !
!------------------------------------------------------------------------------------
  use parameters
  use derived_parameters
  use allocation, only: Nsurf, Nparam, Nrankp, Nrankp1, zpart, zRe, zIm,      &
                        zRe1, zIm1, lnorm, ComplexPlane, EpsZReIm 
  implicit none 
  integer       :: TypeGeom, Npart, dNint, TypeConvTest, Nsurfmax, Nparammax,       &
                   Nrankpmax, Nrankpmax1, Nint, dNintMrank, ipart
  real(O)       :: k, ind_refMed, wavelength, anorm, snorm, epsNint, epsNrank,      &
                   epsMrank, kuser,                 &
                   surf(Npart, Nsurfmax)                                        
  logical       :: DoConvTest, DS, autGenDS, PrnProgress 
  character(80) :: FileTmat  
  complex(O) :: ref1, ref2, ind_ref(Npart)
  complex(O),intent(out) :: btot(2*Nrankpmax*(Nrankpmax+2),2*Nrankpmax*(Nrankpmax+2))
! -----------------------------------------------------------------------------------
!                            Read the input file                                    ! 
! -----------------------------------------------------------------------------------       
  call readinputLAY ( wavelength, ind_refMed, TypeGeom, Npart, anorm,               &
       DoConvTest, DS, autGenDS, Nint, epsNint, epsNrank, epsMrank, dNint,          &
       dNintMrank, FileTmat, PrnProgress, k, snorm, Nsurfmax, Nparammax,            &
       Nrankpmax, Nrankpmax1, TypeConvTest )
  k=kuser
  do ipart = 1, Npart
    Nrankp(ipart)=Nrankpmax
  end do
! -----------------------------------------------------------------------------------
!                                      Main                                         !
! -----------------------------------------------------------------------------------   
  if (.not. DS) then
    call tmatrix_MrankLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,       &
         Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, epsMrank,      &
         FileTmat, PrnProgress, btot)            
  else 
    call tmatrix_MrankDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
         Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nint,      &
         ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank, dNintMrank,        &
         FileTmat, PrnProgress)                
  end if                        
  deallocate (zRe, zIm, zRe1, zIm1)   
  deallocate (Nsurf, Nparam, Nrankp, Nrankp1, zpart, lnorm) 
  deallocate (ComplexPlane, EpsZReIm)   
end subroutine TLAY
!*****************************************************************************
subroutine tlayappendtotmat (b,btot,m,Nrank,Nmax,Nmaxmax)
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
end subroutine tlayappendtotmat

!***********************************************************************************
subroutine readinputLAY ( wavelength, ind_refMed, TypeGeom, Npart, anorm,           &
           DoConvTest, DS, autGenDS, Nint, epsNint, epsNrank, epsMrank, dNint,      &
           dNintMrank, FileTmat, PrnProgress, k, snorm, Nsurfmax, Nparammax,        &
           Nrankpmax, Nrankpmax1, TypeConvTest ) 
  use parameters
  use derived_parameters
  use allocation, only: Nsurf, Nparam, Nrankp, Nrankp1, zpart, zRe, zIm,            &
                        zRe1, zIm1, lnorm, EpsZReIm, ComplexPlane  
  implicit none 
  integer       :: TypeGeom, Npart, NsurfPart, NparamPart, dNint, TypeConvTest,     &
                   i, ipart, isurf, NrankPart, NrankW, Nsurfmax, Nparammax,         &
                   Nrankpmax, Nrankpmax1, Nint, dNintMrank, ios
  real(O)       :: k, ind_refMed, wavelength, anorm, surfPart(NsurfPD), xpart,      &
                   snorm, epsNint, epsNrank, epsMrank, zRePart(NrankPD),            &
                   zImPart(NrankPD), zRePart1(NrankPD), zImPart1(NrankPD),          &
                   OriginPart, x, lnormPart, EpsZReImPart                        
  complex(O)    :: ind_refPartRel
  logical       :: DoConvTest, DS, autGenDS, InputDS, ComplexPlanePart,             &
                   PrnProgress, XFindPar 
  character(80) :: FileTmat, string  
! -----------------------------------------------------------------------------------
!                           Read the input file FileInputLAY                        ! 
! -----------------------------------------------------------------------------------
  call DrvParameters  
  wavelength = 0.1_O * 2._O * Pi 
  ind_refMed = 1._O  
                                                         
  k = 2._O * Pi * ind_refMed / wavelength 
!
  TypeGeom = 1
  anorm = 1._O
  string   = 'GeomProp'

  call check_anorm (anorm)  
  xpart = k * anorm
  snorm = Pi * xpart * xpart 
!
  DoConvTest = .false.    
!                       
!
  DS = .false.
  autGenDS = .true.
!  
  InputDS = .false.
  if (DS .and. .not. autGenDS) InputDS = .true. 
!
  allocate (Nsurf(Npart), Nparam(Npart), zpart(Npart), lnorm(Npart))
  Nsurfmax  = 1
  Nparammax = 1
  do ipart = 1, Npart
    ind_refPartRel = (1.5_O,0._O)
    string = 'OptRegProp'
    call check_ind_ref1 (ipart, ind_refPartRel)
    NsurfPart  = 2
    NparamPart = 1
    OriginPart = 0._O
    do isurf = 1, NsurfPD
      surfPart(isurf) = 1._O
    end do
    lnormPart = 1._O    
    Nsurf(ipart) = NsurfPart
    Nparam(ipart)= NparamPart
    zpart(ipart) = OriginPart
    lnorm(ipart) = lnormPart
    if (Nsurf(ipart)  > Nsurfmax)  Nsurfmax  = Nsurf(ipart)
    if (Nparam(ipart) > Nparammax) Nparammax = Nparam(ipart)
  end do
!
!
  do ipart = 1, Npart  
    NsurfPart  = 2
    NparamPart = 1
    OriginPart = 0._O
    do isurf = 1, NsurfPD
      surfPart(isurf) = 1._O
    end do
    lnormPart = 1._O
  end do  
  call check_geomLAY (TypeGeom, Npart, Nsurf, Nparam)
!
  allocate (ComplexPlane(Npart), EpsZReIm(Npart))
!
  do ipart = 1, Npart
    ComplexPlane(ipart) = .false.
    EpsZReIm(ipart) = 0.95_O
    if (DS .and. autGenDS) then
      ComplexPlanePart = .false.
      EpsZReImPart = 0.95_O     
      ComplexPlane(ipart) = ComplexPlanePart
      EpsZReIm(ipart) = EpsZReImPart
    end if
  end do 
!  
  if (DoConvTest) then
    print "(/,2x,'Convergence Test for a Layered Particle')" 
    print "(  2x,'---------------------------------------')"         
  else
    print "(/,2x,'Convergence Test for a Layered Particle over Mrank')" 
    print "(  2x,'--------------------------------------------------')"    
  end if                          
  allocate (Nrankp(Npart), Nrankp1(Npart)) 
!
!  Nrankpmax  = 0
  Nrankpmax1 = 0
  do ipart = 1, Npart
    if (ipart == 1) then
      x = k * lnorm(ipart)
    else
      x = k * lnorm(ipart)! * abs(ind_ref(ipart))  or the real part of the ref. index
    end if
    NrankW = int(x + 4.05_O * x**0.33_O + 2._O)
    if (.not. DoConvTest .or. InputDS) then
      NrankPart = 17       
      Nrankp(ipart) = NrankPart
      if (ipart == 1) print "(/,2x,'Nrank input values:')"
      print "(2x,'the input value of Nrank for the layer ',i3,' is ', i4,',')",     &
              ipart, NrankPart
      print "(2x, a, i3, a)",                                                       &
     'while the estimated value of Nrank from Wiscombe''s criterion is ', NrankW, ';' 
    else
      if (ipart == 1) print "(/,2x,'Nrank estimates:')"  
      print "(2x,'the estimated value of Nrank from Wiscombe''s criterion')" 
      print "(2x,'for layer ',i2,' is ',i3,';')", ipart, NrankW 
      print "(2x,'- enter the estimated value of Nrank for layer ',i2)", ipart      
      call read_integer (Nrankp(ipart))
    end if
!    if (Nrankp(ipart) > Nrankpmax) Nrankpmax = Nrankp(ipart)
    Nrankp1(ipart) = Nrankp(ipart) - 1
    if (Nrankp1(ipart) > Nrankpmax1) Nrankpmax1 = Nrankp1(ipart)
  end do
!
!  
  if (DoConvTest) then
    print "(/,2x, a)",                                                              &
   '- enter the estimated values of Nint, where Nint = Ndgs * Nrankpmax,'        
    print "(  2x,'  Ndgs = 15,20,..., and Nrankpmax = ',i4,';')", Nrankpmax                 
    call read_integer (Nint) 
  else    
    Nint   = 100
    print "(/,2x,'Nint input value:')"
    print "(  2x, a, i4, a)",                                                       &
   'the input value of Nint is ', Nint, ', while Nint = Ndgs * Nrankpmax,'
    print "(  2x,'Ndgs = 15,20,..., and Nrankpmax = ',i4,';')", Nrankpmax                   
  end if
!
  if (DoConvTest) then     
    print "(/,2x, a)",                                                              &
   '- enter the type of convergence test: 1 - Nint, 2 - Nrank, 3 - Mrank;'
    call read_integerbound (TypeConvTest, 1, 3)     
  else    
    TypeConvTest = 3 
  end if
!
  allocate (zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax), zRe1(Npart,Nrankpmax1),     &
            zIm1(Npart,Nrankpmax1)) 
!
  do ipart = 1, Npart
    do i = 1, Nrankpmax     
      zRe(ipart,i) = 0._O
      zIm(ipart,i) = 0._O
    end do
    do i = 1, Nrankpmax1            
      zRe1(ipart,i) = 0._O
      zIm1(ipart,i) = 0._O
    end do     
    if (DS) then
      call check_MaxNrank (Nrankpmax)
      if (autGenDS) then
          print "( 2x, 'false')"  
!        call zDSLAY (TypeGeom, Npart, Nsurfmax, surf, Nrankpmax, Nrankp, zpart,     &
!             ComplexPlane, EpsZReIm, zRe, zIm)
!        if (TypeConvTest == 2) call zDSLAY (TypeGeom, Npart, Nsurfmax, surf,        &
!                                    Nrankpmax1, Nrankp1, zpart, ComplexPlane,       &
!                                    EpsZReIm, zRe1, zIm1)
      else
        do i = 1, NrankPD
          zRePart(i)  = 0._O
          zImPart(i)  = 0._O
          zRePart1(i) = 0._O
          zImPart1(i) = 0._O
        end do
        do i = 1, Nrankp(ipart)
          zRe(ipart,i) = zRePart(i)  
          zIm(ipart,i) = zImPart(i)               
        end do  
        if (TypeConvTest == 2) then             
          do i = 1, Nrankp1(ipart)
            zRe1(ipart,i) = zRePart1(i) 
            zIm1(ipart,i) = zImPart1(i) 
          end do          
        end if  
      end if                                
    end if
  end do         
!
  epsNint  = 5.e-2
  epsNrank = 5.e-2
  epsMrank = 5.e-2
  dNint = 4
  dNintMrank = 10
!
  FileTmat = '../TMATFILES/T.dat'
!
  PrnProgress = .true.    
end subroutine readinputLAY    
! **********************************************************************************
subroutine tmatrix_MrankLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,       &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zpart, Nint, epsMrank,      &
           FileTmat, PrnProgress, btot)          
  use parameters 
  implicit none
  integer       :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nint,           &
                   Nrankp(Npart), Nparam(Npart), isurf
  real(O)       :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsMrank                
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: PrnProgress
!             
  integer       :: Mstart, Mrank, Nmaxpmax, Nrank, Nmax, Nmaxmax, Nteta,            &
                   i, m, ipart, NthetaConv
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat, Qscat,     &
                   Cext, Qext
  integer,allocatable    :: Nintparam(:,:), Nmaxp(:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
  complex(O),intent(out) :: btot(2*Nrankpmax*(Nrankpmax+2),2*Nrankpmax*(Nrankpmax+2))
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
  Mrank  = Nrankpmax   
  Nrank = 0
  do ipart = 1, Npart        
    if(ipart < Npart) then
      Nrank = Nrank + 2 * Nrankp(ipart)
    else
      Nrank = Nrank + Nrankp(ipart)
    endif
  end do
  allocate (aa(2*Nrank,2*Nrank), bb(2*Nrank,2*Nrank))            
  allocate (a(2*Nrankpmax,2*Nrankpmax), b(2*Nrankpmax,2*Nrankpmax),                 &
            c(2*Nrankpmax,2*Nrankpmax), cv(2*Nrankpmax),cv1(2*Nrankpmax))       
  Nmaxmax = Nrankpmax + Mrank * (2 * Nrankpmax - Mrank + 1)
  allocate (cc(2*Nmaxmax)) 
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do     
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,     &
       Nintparam, paramG, weightsG)     
  Mrank = - 1
  do m = Mstart, Nrankpmax
    call write_1ConvParam (m)
    Mrank = Mrank + 1    
    do ipart = 1, Npart     
      if (Nrankp(ipart) < m) Nrankp(ipart) = m      
    end do
    allocate (Nmaxp(Npart))
    Nmaxpmax = 0
    Nmax = 0
    do ipart = 1, Npart
      if (m == 0) then
        Nmaxp(ipart) = Nrankp(ipart)
      else
        Nmaxp(ipart) = Nrankp(ipart) - iabs(m) + 1
      end if
      if (Nmaxp(ipart) > Nmaxpmax) Nmaxpmax = Nmaxp(ipart)
      if (ipart < Npart) then
        Nmax = Nmax + 2 * Nmaxp(ipart)
      else
        Nmax = Nmax + Nmaxp(ipart)
      end if
    end do
    if (PrnProgress) call write_progress_m (.true., m, 1, 6)  
    call matrix_Q31_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart, Nrankp,       &
         Nmaxpmax, Nmaxp, zpart, m, Nmax, Nint, Nparammax, Nparam, Nintparam,       &
         paramG, weightsG, aa, Nrank, Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 2, 6)  
    call inverse_matrix (aa, 2*Nrank, 2*Nrank, bb, 2*Nrank, 2*Nrank, 2*Nmax)
    if (PrnProgress) call write_progress_m (.false., m, 3, 6)
    call matrix_Q1_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,  &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, a, Nrankpmax, Nrankpmax)
    call extract_matrix3 (1, Nmaxpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)    
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,      &
         2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)     
    if (PrnProgress) call write_progress_m (.false., m, 4, 6)
    call matrix_Q1_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart, Nrankpmax,  &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, c, Nrankpmax, Nrankpmax)        
    call extract_matrix3 (2, Nmaxpmax, bb, Nrank, Nrank, b, Nrankpmax, Nrankpmax)    
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, c, 2*Nrankpmax,      &
         2*Nrankpmax, b, 2*Nrankpmax, 2*Nrankpmax)
    if (PrnProgress) call write_progress_m (.false., m, 5, 6)
    call sum_matrices (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax, 2*Nrankpmax,         &
         c, 2*Nrankpmax, 2*Nrankpmax)
    call incident_matrix_LAY (TypeGeom, k, Nsurfmax, surf, Npart, Nrankpmax,        &
         Nrankp, Nmaxpmax, Nmaxp, zpart, m, Nint, Nparammax, Nparam, Nintparam,     &
         paramG, weightsG, c, Nrankpmax, Nrankpmax)
    call product_matrices (2*Nmaxpmax, 2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,      &
         2*Nrankpmax, c, 2*Nrankpmax, 2*Nrankpmax)
    if (PrnProgress) call write_progress_m (.false., m, 6, 6) 
    call tlayappendtotmat(a, btot, m, Nrankpmax, Nmaxpmax, Nmaxmax)
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,            &
         Nrankpmax, Nmaxpmax, cv)
    call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,             &
         2*Nrankpmax, cv, cv1)
    call extend_vector_positive (cv1, cc, m, Mstart, Nrankpmax, Nmaxpmax, Nmaxmax)                                 
    if (m /= 0) then
      call matrix_m_negativ (Nmaxpmax, Nmaxpmax, a, Nrankpmax, Nrankpmax)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,         &
           Nrankpmax, Nmaxpmax, cv)
      call product_matrix_vector (2*Nmaxpmax, 2*Nmaxpmax, a, 2*Nrankpmax,           &
           2*Nrankpmax, cv, cv1)
      call extend_vector_negative (cv1, cc, m, Nrankpmax, Nmaxpmax, Nmaxmax)
    end if      
    deallocate (Nmaxp)        
    call DSCS (cc, Mrank, Nrankpmax, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,    &
         snorm,.false.,.true., h, v)    
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)    
    call write_DSCS (Nteta,.false., h, v)
!    if (NthetaConv >= int(0.8d0*Nteta)) exit
  end do   
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv, paramG, weightsG,     &
              Nintparam)
end subroutine tmatrix_MrankLAY    
! **********************************************************************************
subroutine tmatrix_MrankDSLAY (TypeGeom, k, ind_ref, snorm, Nsurfmax, surf,     &
           Nparammax, Nparam, Npart, Nrankpmax, Nrankp, zRe, zIm, zpart, Nint,      &
           ComplexPlane, EpsZReIm, autGenDS, epsMrank, epsNrank, dNintMrank,        &
           FileTmat, PrnProgress)    
  use parameters 
  implicit none
  integer       :: TypeGeom, Nsurfmax, Nparammax, Npart, Nrankpmax, Nint,           &
                   Nrankp(Npart), Nparam(Npart), dNintMrank
  real(O)       :: k, surf(Npart,Nsurfmax), zpart(Npart), snorm, epsMrank,          &
                   epsNrank, zRe(Npart,Nrankpmax), zIm(Npart,Nrankpmax),            &
                   EpsZReIm(Npart)
  complex(O)    :: ind_ref(Npart)
  character(80) :: FileTmat
  logical       :: autGenDS, PrnProgress, ComplexPlane(Npart)
!             
  integer       :: Mstart, Mrank, Nteta, Nrank, Nmaxmax, i, m, ipart, NthetaConv,   &
                   NrankAL, NrankpmaxAL, NrankG, NmaxG, NrankpmaxL, NrankpmaxL1
  real(O)       :: alfa, beta, gama, alfap, tetaGI, phiGI, phiGS, Cscat,            &
                   Qscat, Cext, Qext
  logical       :: ComplexDS, ChangeDS, more, continuare, ConvTest
  integer,allocatable    :: Nintparam(:,:), Nrankp1(:)
  real(O),allocatable    :: paramG(:,:,:), weightsG(:,:,:), h(:), v(:), oldh(:),    &
                            oldv(:), zReL(:,:), zImL(:,:), zReL1(:,:), zImL1(:,:)
  complex(O),allocatable :: aa(:,:), bb(:,:), a(:,:), b(:,:), c(:,:), cv(:),        &
                            cv1(:), cc(:)
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
  Mrank  = Nrankpmax 
  Nrank  = 0
  do ipart = 1, Npart - 1
    Nrank = Nrank + 2 * Nrankp(ipart)
  end do     
  Nrank   = Nrank + Nrankp(Npart)      
  NrankAL = Nrank
  NrankpmaxAL = Nrankpmax  
  NrankG  = Nrankpmax 
  Nmaxmax = NrankG + Mrank * (2 * NrankG - Mrank + 1)
  open (unit = iTmat, file = FileTmat, status = 'replace')
  call write_HeadFileTmat (NrankpmaxAL, NrankpmaxAL) 
  call write_TypeConvHead (3)
  call write_2ConvParamReg (Nint, Npart, Nrankp)                  
  allocate (aa(2*NrankAL,2*NrankAL), bb(2*NrankAL,2*NrankAL))            
  allocate (a(2*NrankpmaxAL,2*NrankpmaxAL), b(2*NrankpmaxAL,2*NrankpmaxAL),         &
            c(2*NrankpmaxAL,2*NrankpmaxAL)) 
  allocate (cv(2*NrankG), cv1(2*NrankG))      
  allocate (cc(2*Nmaxmax))   
  allocate (h(Nteta), v(Nteta), oldh(Nteta), oldv(Nteta))
  do i = 1, Nteta
    oldh(i) = 0._O
    oldv(i) = 0._O
  end do    
  allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),           &
            Nintparam(Npart,Nparammax))           
  call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint, Nparammax,     &
       Nintparam, paramG, weightsG)
  ComplexDS = .false.
  do ipart = 1, Npart           
    do i = 1, Nrankp(ipart)                    
      if (zIm(ipart,i) /= 0._O) ComplexDS = .true.
    end do                      
  end do 
  ChangeDS = .false.
! -----------------------------------------------------------------------------------  
! If the sources are distributed in the complex plane and they are automatically    !
! generated, an interactive convergence test over Nint and Nrank can be performed   !
! at each azimuthal mode calculation. To make this test active comment out the next !
! lines.                                                                            !  
! ----------------------------------------------------------------------------------- 
! if (ComplexDS .and. autGenDS) then
!   print "(/,2x, a)",                                                              &
!   '- enter true to change the values of Nint and Nrank at each azimuthal mode and'	                                  
!   print "( 2x, 'false otherwise')"        
!   call read_logical (ChangeDS)    
! end if                               
  Mrank = - 1
  do m = Mstart, NrankG
    call write_1ConvParam (m)
    Mrank = Mrank + 1
    if (m == 0) then
      NmaxG = NrankG
    else
      NmaxG = NrankG - m + 1
    end if
!   --- local convergence test for sources distributed in the complex plane ---    
    NrankpmaxL = Nrankpmax
    allocate (zReL(Npart,NrankpmaxL), zImL(Npart,NrankpmaxL))   
    do ipart = 1, Npart	    
      do i = 1, Nrankp(ipart)
        zReL(ipart,i) = zRe(ipart,i)
        zImL(ipart,i) = zIm(ipart,i)	    	
      end do        
    end do      
    if (ComplexDS .and. m > 1) then
      if (.not. ChangeDS) then
        deallocate (paramG, weightsG, Nintparam)
        Nint = Nint + dNintMrank
        allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),     &
                  Nintparam(Npart,Nparammax))           
        call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint,          &
             Nparammax, Nintparam, paramG, weightsG)     
      else
        print "(/,2x,'Azimuthal mode: m = ',i3)", m  
        more = .true.
        do while (more)
          print "(  2x, '- enter the number of integration points Nint')"
          call read_integer (Nint) 
          deallocate (paramG, weightsG, Nintparam)        
          allocate (paramG(Npart,Nparammax,Nint), weightsG(Npart,Nparammax,Nint),   &
                    Nintparam(Npart,Nparammax))           
          call interpolation_listLAY (TypeGeom, Npart, Nsurfmax, surf, Nint,        &
               Nparammax, Nintparam, paramG, weightsG)                 
          NrankpmaxL = 0  
          Nrank      = 0             
          do ipart = 1, Npart
            if (ComplexPlane(ipart)) then
              print "(2x,'- enter the estimated value of Nrank for region ',i2)",   &
	              ipart
              call read_integer (Nrankp(ipart))		
	    end if	 
            if (Nrankp(ipart)  > NrankpmaxL)  NrankpmaxL  = Nrankp(ipart)
            if (ipart < Npart) then
              Nrank = Nrank + 2 * Nrankp(ipart)
            else
              Nrank = Nrank + Nrankp(ipart)
            end if             	  
          end do 
          if (NrankpmaxL > Nrankpmax) then
            print "(/,2x,'Input error:')"
            print "(  2x, a)",                                                      &
           'the number of discrete sources for the first layer exceeds the number'	    
            print "(  2x, a)",                                                      &
           'of discrete sources corresponding to the initial configuration;'            
            stop   
	  end if                                     
          deallocate (zReL, zImL)	   	   
          allocate (zReL(Npart,NrankpmaxL), zImL(Npart,NrankpmaxL))               
          call zDSLAY (TypeGeom, Npart, Nsurfmax, surf, NrankpmaxL, Nrankp, zpart,  &
               ComplexPlane, EpsZReIm, zReL, zImL)               
          print "(/,2x, a)",                                                        &
         '- enter true to perform a convergence test over the new values of Nrank and'	                                  
          print "( 2x, 'false otherwise')"             
	  call read_logical (ConvTest)		
	  if (ConvTest) then
	    allocate (Nrankp1(Npart))        
	    NrankpmaxL1 = 0    	   
            do ipart = 1, Npart          
	      Nrankp1(ipart) = Nrankp(ipart) - 1	          
	      if (Nrankp1(ipart) > NrankpmaxL1) NrankpmaxL1 = Nrankp1(ipart) 	   
            end do  	     
	    allocate (zReL1(Npart,NrankpmaxL1), zImL1(Npart,NrankpmaxL1)) 
            call zDSLAY (TypeGeom, Npart, Nsurfmax, surf, NrankpmaxL1, Nrankp1,     &
                 zpart, ComplexPlane, EpsZReIm, zReL1, zImL1)  
            call convergence_NrankDSLAY_m (m, TypeGeom, k, ind_ref, snorm,          &
                 Nsurfmax, surf, Nparammax, Nparam, Npart, NrankpmaxL, Nrankp,      &
                 zReL, zImL, NrankpmaxL1, Nrankp1, zReL1, zImL1, zpart, NrankG,     &
                 NmaxG, Nint, Nintparam, paramG, weightsG, epsNrank)	    
            deallocate (zReL1, zImL1, Nrankp1)	     
          end if 	       	  	 	             	  	   	    	   	  	
	  if (ConvTest) then
            print "(2x, a)",                                                        &
           '- enter true for new input values of Nint and Nrank and false to continue;'	
            call read_logical (continuare)
            if (.not. continuare) more = .false.             	      
	  else
            more = .false.
	  end if  
        end do
      end if                                                                               
    end if
    if (Nrank > NrankAL) then
      NrankAL = Nrank
      deallocate (aa, bb)
      allocate (aa(2*NrankAL,2*NrankAL), bb(2*NrankAL,2*NrankAL))                                  
    end if    
!   --- end local convergence test ---     
    if (PrnProgress) call write_progress_m (.true., m, 1, 6)                    
    call matrix_Q31_DS_LAY (TypeGeom, k, ind_ref, Nsurfmax, surf, Npart,            &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, Nrank, Nint, Nparammax,          &
         Nparam, Nintparam, paramG, weightsG, aa, NrankAL, NrankAL)         
    if (PrnProgress) call write_progress_m (.false., m, 2, 6)                                            
    call inverse_matrix (aa, 2*NrankAL, 2*NrankAL, bb, 2*NrankAL, 2*NrankAL, 2*Nrank)
    if (PrnProgress) call write_progress_m (.false., m, 3, 6)        
    call matrix_Q1_DS_LAY (TypeGeom, 1, k, ind_ref, Nsurfmax, surf, Npart,          &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, NrankG, NmaxG, Nint, Nparammax,  &
         Nparam, Nintparam, paramG, weightsG, a, NrankpmaxAL, NrankpmaxAL)                               
    call extract_matrix3 (1, NrankpmaxL, bb, NrankAL, NrankAL, b,                   &    
         NrankpmaxAL, NrankpmaxAL)                
    call product_matrices (2*NmaxG, 2*NrankpmaxL, 2*NrankpmaxL, a, 2*NrankpmaxAL,   &
         2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)                       
    if (PrnProgress) call write_progress_m (.false., m, 4, 6)
    call matrix_Q1_DS_LAY (TypeGeom, 3, k, ind_ref, Nsurfmax, surf, Npart,          &
         NrankpmaxL, Nrankp, zReL, zImL, zpart, m, NrankG, NmaxG, Nint, Nparammax,  &
         Nparam, Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)                      
    call extract_matrix3 (2, NrankpmaxL, bb, NrankAL, NrankAL, b,                   &
         NrankpmaxAL, NrankpmaxAL)               
    call product_matrices (2*NmaxG, 2*NrankpmaxL, 2*NrankpmaxL, c, 2*NrankpmaxAL,   &
         2*NrankpmaxAL, b, 2*NrankpmaxAL, 2*NrankpmaxAL)             
    if (PrnProgress) call write_progress_m (.false., m, 5, 6)    
    call sum_matrices (2*NmaxG, 2*NrankpmaxL, a, 2*NrankpmaxAL, 2*NrankpmaxAL,      &
         c, 2*NrankpmaxAL, 2*NrankpmaxAL)         
    call incident_matrix_DS_LAY (TypeGeom, k, Nsurfmax, surf, Npart, NrankpmaxL,    &
         Nrankp, zReL, zImL, zpart, m, NrankG, NmaxG, Nint, Nparammax, Nparam,      &
         Nintparam, paramG, weightsG, c, NrankpmaxAL, NrankpmaxAL)                                         
    call product_matrices (2*NmaxG, 2*NrankpmaxL, 2*NmaxG, a, 2*NrankpmaxAL,        &
         2*NrankpmaxAL, c, 2*NrankpmaxAL, 2*NrankpmaxAL)         
    if (PrnProgress) call write_progress_m (.false., m, 6, 6)
    call write_FileTmat (NrankpmaxAL, NrankpmaxAL, a)                    
    call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, m,            &
         NrankG, NmaxG, cv)
    call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL, 2*NrankpmaxAL,  &
         cv, cv1)
    call extend_vector_positive (cv1, cc, m, Mstart, NrankG, NmaxG, Nmaxmax)                                 
    if (m /= 0) then
      call matrix_m_negativ (NmaxG, NmaxG, a, NrankpmaxAL, NrankpmaxAL)
      call PWcoefficients_ab_m (tetaGI, phiGI, alfa, beta, gama, alfap, -m,         &
           NrankG, NmaxG, cv)
      call product_matrix_vector (2*NmaxG, 2*NmaxG, a, 2*NrankpmaxAL,               &
           2*NrankpmaxAL, cv, cv1)
      call extend_vector_negative (cv1, cc, m, NrankG, NmaxG, Nmaxmax)
    end if      
    call DSCS (cc, Mrank, NrankG, Nmaxmax, Nteta, phiGS, alfa, beta, gama, k,       &
         snorm,.false.,.true., h, v)                                        
    call delta_DSCS (Nteta, h, v, oldh, oldv, epsMrank, NthetaConv)    
    call write_DSCS (Nteta,.false., h, v)
    deallocate (zReL, zImL)     
    if (NthetaConv >= int(0.8*Nteta)) exit
  end do 
  close (unit = iTmat)  
  call CQscat (cc, Mrank, NrankG, Nmaxmax, k, snorm, Cscat, Qscat)
  call CQext (cc, Mrank, NrankG, Nmaxmax, tetaGI, phiGI, alfa, beta, gama,          &
       alfap, k, snorm, Cext, Qext)      
  call write_Effic (Qscat, Qext)             
  call write_MrankConvRes (NthetaConv, epsMrank)
  if (NthetaConv >= int(0.8*Nteta)) then
    print "(/,2x,'Convergence criterion for Mrank is satisfied;')"                                                              
  else
    print "(/,2x,'Convergence criterion for Mrank is not satisfied;')"
  end if 
  call write_InfoFileTmat (FileTmat, Mrank, NrankG, .true., .false., .false.)
  call ScatCharact (k, FileTmat, Mrank, NrankG, .true., .false., .false.) 
  print "(/,2x,'T matrix is stored in ',a50)", FileTmat
  print "(  2x,'The dimensions of the T matrix are given by:')"
  print "(  2x,'- maximum expansion order,   Nrank = ',i3,',')", NrankG
  print "(  2x,'- number of azimuthal modes, Mrank = ',i2,';')", Mrank            
  deallocate (aa, bb, a, b, c, cv, cv1, cc, h, v, oldh, oldv, paramG, weightsG,     &
              Nintparam)
end subroutine tmatrix_MrankDSLAY
