C**********************************************************************C
C**********************************************************************C
C                                                                      C
C           IMPLEMENTATION OF LDPM FOR ABAQUS                          C
C                                                                      C
C           BY Erol Lale                                               C
C              Roozbeh Rezakhani                                       C
C              Gianluca Cusatis                                        C
C              Yuhui Lyu                                               C
C              Matthew Troemner                                        C
C              Elham Ramyar                                            C
C                                                                      C
C              30 August 2020                                          C
C                                                                      C
C**********************************************************************C
C**********************************************************************C





C======================================================================C
C     Notes                                                            C
C======================================================================C
C     Details regarding definition of double precision are available   C
C     here http://fortranwiki.org/fortran/show/Real+precision and in   C
C     Lemmon and Schafer, 2005, p. 20–22.                              C
C======================================================================C 





C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     LDPM MODULE                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      ! nomaxel =                                                      !
      ! FacetVar = Contains facet geometry                             !
      ! elecon = Element connectivity                                  !
      ! svars_selection = Contains svars selection                     !
      ! svars_freq (3) = Contains frequency to print svars data into file  !
      ! svars_nodal_stress =Flag for nodal stress computation 
      ! print_binary = Contains if you want to print data in binary or not
      !================================================================!


      module ModuleLDPM
      implicit none 
      
      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)

      ! Declare integers
      integer                                 :: nomaxel
      integer                                 :: svars_nodal_stress
      integer,  dimension(:,:),   allocatable :: elecon                
      integer,  dimension(:),     allocatable :: svars_selection 
      real(dp), dimension(:),     allocatable :: next_print      
      integer print_binary
      ! Declare double precision floats
      real(dp)                                :: b
      real(dp)                                :: svars_freq(3)
      real(dp), dimension(:,:,:), allocatable :: FacetVar           

      end module



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     USER ELEMENT SUBROUTINE                                          C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of the variables passed to VUEL                    !
      !================================================================!
      ! mcrd = dimension of the problem (3 for 3D)
      ! nnode = 4 (number of nodes per element)
      ! ndofel = 24 (number of dofs per element)
      ! nblock = number of block of elements (ABAQUS environment variable)
      ! rhs(nblock,ndofel) = right hand side vector in the global system of reference (OUTPUT)
      ! amass(nblock,ndofel,ndofel) = mass matrix (diagonal for explicit) in the global system of reference (OUTPUT)
      ! dtimeStable(nblock) = stable time step (INPUT)
      ! svars(nblock,nsvars) = state variables (INPUT/OUTPUT)
      ! nsvars = 12 * number of state variables for one facet (12 facets per tet) 
      ! energy(noblock,12) = various energy measures, see ABAQUS manual and code below (OUTPUT)
      ! props(nprops) = floating constant model parameters (INPUT)
      ! nprops = number of floating constant model parameters
      ! jprops = integer constant model parameters (INPUT)
      ! njprops = number of integer constant model parameters
      ! coords(nblock,nnnode,mcrd) = nodal coordinates (INPUT)
      ! u(nblock,ndofel) = displacement in the global system of reference (INPUT)
      ! du(nblock,ndofel) = displacement increment in the global system of reference (INPUT)
      ! v(nblock,ndofel) = velocity in the global system of reference (INPUT)
      ! a(nblock,ndofel) = accelleration in the global system of reference (INPUT)
      ! jtype = element type
      ! jElem(nblock) = Element ID in the global numbering
      ! time(2) = (1) current step time, (2) total time (INPUT)
      ! period = ?
      ! dtimeCur = current time increment (INPUT)
      ! dtimePrev = previous time increment (INPUT)
      ! kstep = current step ID
      ! kinc  = current time increment ID
      ! lflags = type of analysis (see ABAQUS manual and code below)
      ! dMassScaleFactor(nblock) = It is for mass scaling, ???? 
      ! predef(nblock, nnode, npredef, npred) = nodal predefined fields (for example T, h, p, ...)
      ! npredef = number of fields (first one is Temperature, then others)
      ! npred = 2, (1) new, (2) old
      ! jdltyp = load type (body forces, ...) ???
      ! adlmag = load magnitude ???          
      !           
      ! ndofn  ... number of degrees of freedom per node
      ! nshr   ... number of shear stress component
      ! ntens  ... total number of stress tensor components (=ndi+nshr)
      ! ninpt  ... number of integration points
      ! nsvint ... number of state variables per facet (stress, strain and other ones)
      !================================================================!


      SUBROUTINE VUEL(nblock,rhs,amass,dtimeStable,svars,nsvars,
     1                energy,
     2                nnode,ndofel,props,nprops,jprops,njprops,
     3                coords,mcrd,u,du,v,a,
     4                jtype,jElem,
     5                time,period,dtimeCur,dtimePrev,kstep,kinc,
     6                lflags,
     7                dMassScaleFactor,
     8                predef,npredef,
     9                jdltyp, adlmag)


      use ModuleLDPM
      implicit none  


      INTEGER(4) :: ORIGINAL_FPE_FLAGS, NEW_FPE_FLAGS
      integer(KIND=4),  parameter :: ndime                = 3
      integer(KIND=4),  parameter :: ndofn                = 6    
      integer(KIND=4),  parameter :: nshr                 = 3  
      integer(KIND=4),  parameter :: nsvint               = 25         
      integer(KIND=4),  parameter :: ntens                = 6
      integer(KIND=4),  parameter :: ndi                  = 3
      integer(KIND=4),  parameter :: nface                = 12

      integer(KIND=4),  parameter :: jMassCalc            = 1
      integer(KIND=4),  parameter :: jIntForceAndDtStable = 2  
      integer(KIND=4),  parameter :: jExternForce         = 3

      integer(KIND=4),  parameter :: iProcedure           = 1    
      integer(KIND=4),  parameter :: iNlgeom              = 2  
      integer(KIND=4),  parameter :: iOpCode              = 3   
      integer(KIND=4),  parameter :: nFlags               = 3 

      integer(KIND=4),  parameter :: jDynExplicit         = 17

      ! energy array index parameters
      integer(KIND=4),  parameter :: iElPd                = 1
      integer(KIND=4),  parameter :: iElCd                = 2
      integer(KIND=4),  parameter :: iElIe                = 3 
      integer(KIND=4),  parameter :: iElTs                = 4
      integer(KIND=4),  parameter :: iElDd                = 5
      integer(KIND=4),  parameter :: iElBv                = 6
      integer(KIND=4),  parameter :: iElDe                = 7
      integer(KIND=4),  parameter :: iElHe                = 8
      integer(KIND=4),  parameter :: iElKe                = 9
      integer(KIND=4),  parameter :: iElTh                = 10
      integer(KIND=4),  parameter :: iElDmd               = 11
      integer(KIND=4),  parameter :: nElEnergy            = 12

      ! predefined variables indices
      integer(KIND=4),  parameter :: iPredValueNew        = 1
      integer(KIND=4),  parameter :: iPredValueOld        = 2
      integer(KIND=4),  parameter :: nPred                = 2

      ! time indices
      integer(KIND=4),  parameter :: iStepTime            = 1
      integer(KIND=4),  parameter :: iTotalTime           = 2
      integer(KIND=4),  parameter :: nTime                = 2

  
      ! Declare integers
      integer(KIND=4)                         :: fMat
      integer(KIND=4)                         :: i
      integer(KIND=4)                         :: icol
      integer(KIND=4)                         :: icpu
      integer(KIND=4)                         :: idfile
      integer(KIND=4)                         :: idime
      integer(KIND=4)                         :: idofel
      integer(KIND=4)                         :: idofn
      integer(KIND=4)                         :: ielem
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: inode
      integer(KIND=4)                         :: ios
      integer(KIND=4)                         :: isvinc
      integer(KIND=4)                         :: jcol
      integer(KIND=4)                         :: jdltyp
      integer(KIND=4)                         :: jtype
      integer(KIND=4),  dimension(njprops)    :: jprops
      integer(KIND=4)                         :: jnode
      integer(KIND=4),  dimension(nblock)     :: jElem        
      integer(KIND=4)                         :: kblock
      integer(KIND=4)                         :: Kfindloc
      integer(KIND=4)                         :: KINC
      integer(KIND=4)                         :: KPROCESSNUM  
      integer(KIND=4)                         :: KSTEP 
      integer(KIND=4)                         :: LOP
      integer(KIND=4)                         :: LRESTART
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer(KIND=4),  dimension(nFlags)     :: lflags 
      integer(KIND=4)                         :: mcrd
      integer(KIND=4)                         :: nblock
      integer(KIND=4)                         :: ndofel
      integer(KIND=4)                         :: njprops
      integer(KIND=4)                         :: nnode
      integer(KIND=4)                         :: ninpt
      integer(KIND=4)                         :: npredef
      integer(KIND=4)                         :: nprops
      integer(KIND=4)                         :: NUMPROCESSES
      integer(KIND=4)                         :: nsvars
      integer(KIND=4)                         :: svarsidfile 


      ! Declare double precision floats
      real(dp), dimension(nblock, ndofel)     :: a
      real(dp), dimension(1,3)                :: accelVec
      real(dp), dimension(nblock)             :: adlmag     
      real(dp)                                :: al
      real(dp)                                :: aLe
      real(dp)                                :: aLen
      real(dp), dimension(nblock,ndofel,ndofel) :: amass  
      real(dp), dimension(ndofel)             :: BodyL      
      real(dp), dimension(mcrd)               :: coordgp   
      real(dp), dimension(nblock,nnode,mcrd)  :: coords           
      real(dp), dimension(mcrd)               :: cross
      real(dp), dimension(mcrd)               :: cross0
      real(dp)                                :: cur_step
      real(dp), dimension(ndofel)             :: DeltaU
      real(dp)                                :: density  
      real(dp), dimension(ndofel)             :: deigen      
      real(dp)                                :: diff
      real(dp), dimension(nblock)             :: dMassScaleFactor
      real(dp), dimension(ntens)              :: dstran
      real(dp)                                :: DTIME
      real(dp)                                :: dtimeCur
      real(dp)                                :: dtimePrev
      real(dp), dimension(nblock)             :: dtimeStable
      real(dp), dimension(nblock,ndofel)      :: du      
      real(dp)                                :: dx
      real(dp)                                :: dy
      real(dp)                                :: dz
      real(dp), dimension(mcrd,nnode)         :: elcoord0
      real(dp), dimension(mcrd,nnode)         :: elcoord1
      real(dp), dimension(mcrd,nnode)         :: elcoord
      real(dp), dimension(ndofel,ndofel)      :: elmass
      real(dp), dimension(ndofel,ndofel)      :: elmassInvStiff
      real(dp), dimension(ndofel,ndofel)      :: elstif
      real(dp)                                :: endTime
      real(dp), dimension(nblock,nElEnergy)   :: energy
      real(dp)                                :: epsV
      real(dp)                                :: factorStable  
      real(dp), dimension(nsvint)             :: facestatev
      real(dp)                                :: fcstrainEN
      real(dp), dimension(ndofn)              :: force  
      real(dp), dimension(ndofel)             :: forceVec
      real(dp)                                :: omax
      real(dp)                                :: period
      real(dp), dimension(nblock,nnode,npredef,nPred)  :: predef    
      real(dp), dimension(nprops)             :: props
      real(dp), dimension(nblock,ndofel)      :: rhs   
      real(dp), dimension(nsvint)             :: statevLocal
      real(dp), dimension(ntens)              :: stress
      real(dp), dimension(ntens)              :: stran      
      real(dp), dimension(nblock,nsvars)      :: svars  
      real(dp)                                :: tetVol0
      real(dp)                                :: tetVol      
      real(dp), dimension(nTime)              :: time
      real(dp)                                :: startTime
      real(dp)                                :: time_step
      real(dp), dimension(nblock,ndofel)      :: u
      real(dp), dimension(nblock,ndofel)      :: v
      real(dp), dimension(ndofel,ndofel)      :: veigen
      real(dp), dimension(mcrd)               :: va
      real(dp), dimension(mcrd)               :: vb
      real(dp), dimension(mcrd)               :: vc
      real(dp), dimension(mcrd)               :: va0
      real(dp), dimension(mcrd)               :: vb0
      real(dp), dimension(mcrd)               :: vc0
      real(dp), dimension(3)                  :: vnr
      real(dp), dimension(3)                  :: vmr
      real(dp), dimension(3)                  :: vlr


      ! Declare booleans
      logical(KIND=4) direxist
      logical(KIND=4) LargeDefFlag
      logical(KIND=4) RateEffectFlag


      ! Declare characters
      character(80)                           :: Chemin
      character(256)                          :: filename    
      character(80)                           :: JOBNAME
      character(256)                          :: OUTDIR      
      character(256)                          :: path
      character(1)                            :: path_separator       
      character(1)                            :: separator
      character(256)                          :: svarsfile
      character(256)                          :: svarsfilepath


      ! Assign values
      LargeDefFlag                     = .false.  
      RateEffectFlag                   = .false.
      factorStable                     = 0.900d0


      ! Extract job data
      CALL VGETJOBNAME( JOBNAME, LENJOBNAME ) ! Jobname
      CALL VGETOUTDIR( OUTDIR, LENOUTDIR )    ! Work directory

      CALL getcwd(path)

      separator=path(1:1)

      if (separator.eq.'/') then
         path_separator = '/'
      else
         path_separator = '\'
      endif


      ! Store working directory for error logs
      Chemin=trim(trim(OUTDIR)//path_separator//trim(JOBNAME))

      ! Store (and create) svar directory
      svarsfilepath=trim(OUTDIR(1:LENOUTDIR)//path_separator//trim(JOBNAME)//path_separator)
      inquire(directory=trim(svarsfilepath),exist=direxist)
      if (.not.direxist) then
         call system('mkdir '//trim(svarsfilepath))
      endif
  

C======================================================================C
C     MPI Files                                                        C
C======================================================================C

      ! Extract simulation parameters
      CALL VGETNUMCPUS( NUMPROCESSES ) 
      CALL VGETRANK( KPROCESSNUM ) 
      icpu=KPROCESSNUM

      ! If MPI, make multiple erl files    
      if(NUMPROCESSES.gt.0) then         
         idfile=105+10*icpu
         WRITE(filename, '(a,i3.3,a)') 
     1      trim(Chemin), icpu, ".erl"
         filename=trim(filename)     
         open(unit=idfile, FILE=filename, action='WRITE', access='APPEND', IOSTAT=ios)
      endif



C======================================================================C
C     Model Checks                                                     C
C======================================================================C

      density = Props(1) ! Used for time step calc

      ! Preliminaries    
      if(jtype.ne.4) then
         write(*,*)'Incorrect element type'
         !stop !call xit
      endif 

      if(nsvars.lt.nface*nsvint) then
         write(*,*)'Increase the number of SDVs to',nface*nsvint
         !stop !call xit
      endif 
      
      ! Find number of elements          
      if (time(iTotalTime).eq.0.0d0) then 
         if (allocated(FacetVar).eq.0.0d0) then
            nomaxel=jprops(1)                  
            allocate(FacetVar(30,12,nomaxel), elecon(nomaxel,5))
            call VUEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
         endif            

C      allocate an array next_print to record the next print time for each element
          if (allocated(next_print).eq.0) then
              nomaxel=jprops(1)                  
              allocate(next_print(nomaxel))
              next_print=0.0
          endif 
      endif
      
      ! Make svar files
      if (size(svars_selection).gt.0) then
        write(svarsfile, '(a,a,i3.3,a)') trim(svarsfilepath), 'svars-', icpu, ".txt"
        svarsidfile=19+icpu
        if (print_binary.eq.1) then
            open(unit=svarsidfile, FILE = svarsfile, action = "WRITE", access='APPEND', IOSTAT=ios, form='unformatted')
        else
            open(unit=svarsidfile, FILE = svarsfile, action = "WRITE", access='APPEND', IOSTAT=ios)
        endif
      endif
C======================================================================C
C     Select Procedures                                                C
C======================================================================C

      !Procedure type: Explicit dynamic 
      if (lflags(iProcedure).eq.jDynExplicit) then 

        ! Mass calculation
        if (lflags(iOpCode).eq.jMassCalc) then               


C======================================================================C   
C     CALCULATE MASS MATRIX                                            C
C======================================================================C      

           density = Props(1)
           do kblock=1,nblock              
               ielem=jelem(kblock) 
                do inode=1,nnode                
                    do idime=1,mcrd
                        elcoord(idime, inode)=coords(kblock,inode,idime)
                    end do
                end do   
                elmass=0.0d0              
                call KMASSLDPM(nnode,ndofn,ndofel,ninpt,mcrd,density,
     1                        elcoord,elecon(ielem,2:5),FacetVar(:,:,ielem),elmass)  
!               call KMASSCONTINUUM(nnode,ndofn,ndofel,4,mcrd,density,
!    1                             elcoord,nconnec,faceprop,elmass)   
                  
                amass(kblock,:,:)=elmass(:,:)              
           end do 
           

        else if (lflags(iOpCode).eq.jIntForceAndDtStable) then
    

C======================================================================C   
C     CALCULATE INTERNAL FORCES
C======================================================================C

            do kblock = 1, nblock                
               ielem=jelem(kblock) 
               ! Get displacement increment for ielem
               DeltaU(1:ndofel) = du(kblock,1:ndofel)
               do inode=1,nnode
                    icol=(inode-1)*ndofn
                    do idime=1,mcrd
                       elcoord0(idime, inode) = coords(kblock,inode,idime)
                       elcoord(idime, inode) = coords(kblock,inode,idime)+u(kblock,icol+idime)
                    end do                     
               end do    
               
               ale=svars(kblock,12*nsvint+1)
               if(time(iTotalTime).eq.0.0d0) then
                  do inode=1,nnode
                      do jnode=inode,nnode
                          if (inode.eq.jnode) cycle
                          dx=elcoord(1,inode)-elcoord(1,jnode)
                          dy=elcoord(2,inode)-elcoord(2,jnode)
                          dz=elcoord(3,inode)-elcoord(3,jnode)
                          al=DSQRT(DABS(dx*dx+dy*dy+dz*dz))               
                          if((al.lt.aLe).or.((aLe.lt.EPSILON(aLe)).and.(aLe.gt.-EPSILON(aLe)))) then
                              aLe=al                         
                          endif                       
                      end do                   
                  end do
                  svars(kblock,12*nsvint+1)=ale
               end if    


C======================================================================C
C     Element volumetric strain calculation
C======================================================================C

               do idime=1,mcrd
                  va0(idime) = elcoord0(idime,2) - elcoord0(idime,1)
                  vb0(idime) = elcoord0(idime,3) - elcoord0(idime,1)
                  vc0(idime) = elcoord0(idime,4) - elcoord0(idime,1)
                  va(idime) = elcoord(idime,2) - elcoord(idime,1)
                  vb(idime) = elcoord(idime,3) - elcoord(idime,1)
                  vc(idime) = elcoord(idime,4) - elcoord(idime,1)
               end do 
      
               cross0(1) = vb0(2) * vc0(3) - vb0(3) * vc0(2)
               cross0(2) = vb0(3) * vc0(1) - vb0(1) * vc0(3)
               cross0(3) = vb0(1) * vc0(2) - vb0(2) * vc0(1) 
               cross(1)  = vb(2) * vc(3) - vb(3) * vc(2)
               cross(2)  = vb(3) * vc(1) - vb(1) * vc(3)
               cross(3)  = vb(1) * vc(2) - vb(2) * vc(1)
               tetVol0 = dabs(va0(1)*cross0(1)+va0(2)*cross0(2)+va0(3)*cross0(3))/6.0d0
               tetVol = dabs(va(1)*cross(1)+va(2)*cross(2)+va(3)*cross(3))/6.0d0     
      
               if ((tetVol.lt.1.0d-10).or.(tetVol0.lt.1.0d-10)) then
                  write(*,*)'Warning: Tet vols in critical range; continuation may produce numerical issues.' 
                  write(*,*)'Setting tet volumes to 1.0E-11.'   
                  tetVol = 1.0d-11
                  tetVol0 = 1.0d-11
                  epsV = 0.0d0
               else if (((tetVol-tetVol0).lt.EPSILON(tetVol)).and.((tetVol-tetVol0).gt.-EPSILON(tetVol))) then
                  epsV = 0.0d0
               else
                  epsV = (tetVol-tetVol0)/(3.0d0*tetVol0)
               endif


C======================================================================C
C     Loop over facets
C======================================================================C

               Elstif=0.0d0
               fcstrainEN=0.0d0 
               forceVec=0.0d0
               aLe=0.0d0
               
               do iface=1,12 


C======================================================================C
C     Initial values of facets state variables
C======================================================================C

                  isvinc= (iface-1)*nsvint 
                  do i=1,nsvint                        
                      facestatev(i)=svars(kblock,isvinc+i)
                  end do   
                             
                  fMat=IDNINT(FacetVar(30,iface,ielem))
                 
                  call SUBLDPM (nprops,mcrd,nnode,ndofn,ndofel,ntens,
     1                          nsvint,elcoord0,elcoord,DeltaU,props,
     2                          FacetVar(:,iface,ielem),facestatev,
     3                          elecon(ielem,2:5),forceVec,Elstif,
     4                          epsV,fcstrainEN,alen,ielem,fMat,LargeDefFlag,
     5                          vnr,vmr,vlr,RateEffectFlag)
                   
                  if ((aLe.eq.0.0d0).or.(alen.lt.ale)) aLe=aLen
                       

C======================================================================C
C     Update state variables
C======================================================================C

              do i=1,nsvint                        
                 svars(kblock,i+isvinc)=facestatev(i)
              end do
              
              startTime = svars_freq(1)
              time_step = svars_freq(2)
              endTime = svars_freq(3)
              if (size(svars_selection).ne.0) then
                  if ((time(1).ge.next_print(ielem)).OR.(time(1).eq.startTime).OR.(time(1).eq.endTime)) then
                        inode=int(FacetVar(1,iface,ielem))
                        jnode=int(FacetVar(2,iface,ielem))
                        !inode=Kfindloc(elecon(ielem,2:5),inode,4)
                        !jnode=Kfindloc(elecon(ielem,2:5),jnode,4)
                        icol=(inode-1)*ndofn
                        jcol=(jnode-1)*ndofn
                        if (print_binary.eq.1) then
                            write(svarsidfile) time, real(ielem, 8), real(iface, 8)  ! time ielem iface
                            write(svarsidfile) u(kblock, icol+1:icol+6)  ! ui
                            write(svarsidfile) u(kblock, jcol+1:jcol+6)  ! uj
                        else
                            write(svarsidfile,*) time, ielem, iface ! time ielement iface 
                             write(svarsidfile,*) u(kblock, icol+1:icol+6)  ! ui
                            write(svarsidfile,*) u(kblock, jcol+1:jcol+6)  ! uj
                        endif
                        if (svars_nodal_stress.eq.1) then
                            if (print_binary.eq.1) then
                                write(svarsidfile) FacetVar(11,iface,ielem)  ! area
                                write(svarsidfile) vnr  ! vnr
                                write(svarsidfile) vmr  ! vmr
                                write(svarsidfile) vlr  ! vlr
                            else
                                write(svarsidfile,*) FacetVar(11,iface,ielem)  ! area
                                write(svarsidfile,*) vnr  ! vnr
                                write(svarsidfile,*) vmr  ! vmr
                                write(svarsidfile,*) vlr  ! vlr
                            endif
                        endif
C                       write(svarsidfile,*) facestatev
                        do i=1,size(svars_selection)
                            if (i.eq.size(svars_selection)) then
                                ! Write and change to new line in output svars file
                                if (print_binary.eq.1) then
                                    write(svarsidfile) facestatev(svars_selection(i))
                                else
                                    write(svarsidfile,'(1x,E32.16E3)') facestatev(svars_selection(i))
                                endif
                            else
                                ! Write and DON'T change to new line in output svars file (adding $ sign to the format)
                                if (print_binary.eq.1) then
                                    write(svarsidfile) facestatev(svars_selection(i))
                                else
                                    write(svarsidfile,'(1x,E32.16E3,$)') facestatev(svars_selection(i))
                                endif
                            endif
                        end do
                    endif                       
                endif
C======================================================================C

               end do !END of LOOP OVER FACETS

               ! Update the force vector        
               do idofn = 1,ndofel                           
                       rhs(kblock,idofn) = rhs(kblock,idofn) + forceVec(idofn)
               end do

               ! Internal energy calculation 
               energy(kblock,iElIe) = energy(kblock,iElIe) + fcstrainEN


C======================================================================C
C     Stable time calculations
C======================================================================C
        
               if (time(iTotalTime).eq.0.0d0) then
               call cal_stable_time(nnode,ndofn,ndofel,ninpt,mcrd,density,
     1                              elcoord0,elecon(ielem,2:5),FacetVar(:,:,ielem),elstif,omax)
             
               dtimeStable(kblock) = factorStable*omax
               end if
               
            if (size(svars_selection).ne.0) then
             if ((time(2).ge.next_print(ielem)).OR.(time(2).eq.startTime).OR.(time(2).eq.endTime)) then
                next_print(ielem) = next_print(ielem) + time_step
             endif
            endif               
            end do   ! END OF KBLOCK LOOP

        elseif(lflags(iOpCode) .eq. jExternForce) then   


C======================================================================C
C     Calculate external forces
C======================================================================C

      if ((jdltyp.eq.101).or.(jdltyp.eq.102).or.(jdltyp.eq.103)) then
      
      do kblock = 1, nblock
        ielem = jelem(kblock)         
        do inode=1,nnode
            do idime=1,mcrd
                elcoord(idime, inode)=coords(kblock,inode,idime)
            end do            
        end do   
        
        accelVec(1,:) = 0.0d0
        select case(jdltyp)            
            case(101)
                accelVec(1,1) = adlmag(kblock)
            case(102)
                accelVec(1,2) = adlmag(kblock)
            case(103)
                accelVec(1,3) = adlmag(kblock)
            end select            
            
        call BODYF_LDPM(nnode,ndofn,ndofel,ninpt,mcrd,density,elcoord,
     1            elecon(ielem,2:5),FacetVar(:,:,ielem),accelVec,BodyL) 
                    
          do idofel=1,ndofel
              rhs(kblock,idofel) = BodyL(idofel)
          end do
          
      end do   ! kblock loop
      
      endif
    
        
        end if  ! end of iOpCode if statement

      end if  ! end of iProcedure if statement
  

      close(NUMPROCESSES)
      if (NUMPROCESSES.gt.0) close(idfile)
950     format (A7,2X,I9,2X, I5, 2X, I5, 2X, I5, 2X, E16.8)
      
      return
      end



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     LDPM SUBROUTINNE                                                 C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


C======================================================================C
C     Notes                                                            C
C======================================================================C
C     This subroutine evaluates the strain-displacement                C
C     matrix for the 4-node tetrahedron                                C
C======================================================================C 


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !      
      !================================================================!


      subroutine SUBLDPM (nprops,ndime,nnode,ndofn,ndofel,ntens,nsvint,
     1                     coords0,coords,du,props,faceprop,facestatev,
     2                     nconnec,forceVec,Elstif,epsV,strainEN,alen,
     3                     jelem,fMat,LargeDefFlag,vnr,vmr,vlr,
     4                     RateEffectFlag)


      implicit none  


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare integers
      integer(KIND=4)                         :: fMat
      integer(KIND=4)                         :: i
      integer(KIND=4), dimension(3,3)         :: Iden
      integer(KIND=4)                         :: Icnt
      integer(KIND=4)                         :: icol
      integer(KIND=4)                         :: idof
      integer(KIND=4)                         :: idofn
      integer(KIND=4)                         :: inode 
      integer(KIND=4)                         :: irow
      integer(KIND=4)                         :: Jcnt
      integer(KIND=4)                         :: jcol
      integer(KIND=4)                         :: jdof
      integer(KIND=4)                         :: jdofn
      integer(KIND=4)                         :: jElem
      integer(KIND=4)                         :: jnode
      integer(KIND=4)                         :: jrow
      integer(KIND=4)                         :: Kfindloc
      integer(KIND=4), dimension(4)           :: nconnec
      integer(KIND=4)                         :: ndofn
      integer(KIND=4)                         :: nProps
      integer(KIND=4)                         :: ndime
      integer(KIND=4)                         :: ndofel
      integer(KIND=4)                         :: nnode
      integer(KIND=4)                         :: ntens
      integer(KIND=4)                         :: nsvint


      ! Declare double precision floats
      real(dp), dimension(3,6)                :: Ai
      real(dp), dimension(3,6)                :: Aj
      real(dp)                                :: AreLE_N
      real(dp)                                :: AreLE_T
      real(dp)                                :: alpha
      real(dp)                                :: alen
      real(dp)                                :: area     
      real(dp), dimension(6)                  :: BjkN
      real(dp), dimension(6)                  :: BikN
      real(dp), dimension(6)                  :: BjkM
      real(dp), dimension(6)                  :: BikM
      real(dp), dimension(6)                  :: BjkL
      real(dp), dimension(6)                  :: BikL
      real(dp), dimension(3,ndofn)            :: bmat
      real(dp), dimension(ndime, nnode)       :: coords
      real(dp), dimension(ndime, nnode)       :: coords0
      real(dp)                                :: cosine
      real(dp)                                :: depsN
      real(dp)                                :: depsM
      real(dp)                                :: depsL
      real(dp), dimension(3)                  :: dn
      real(dp), dimension(3)                  :: dm
      real(dp), dimension(3)                  :: dl
      real(dp)                                :: dtime
      real(dp), dimension(ndofel,1)           :: du
      real(dp), dimension(ndofn)              :: dui
      real(dp), dimension(ndofn)              :: duj
      real(dp)                                :: dx
      real(dp)                                :: dy
      real(dp)                                :: dz
      real(dp)                                :: epsN
      real(dp)                                :: epsM
      real(dp)                                :: epsL  
      real(dp)                                :: epsV    
      real(dp)                                :: E0
      real(dp), dimension(ndofel,ndofel)      :: Elstif
      real(dp)                                :: EN
      real(dp)                                :: ET
      real(dp), dimension(30)                 :: faceprop      
      real(dp), dimension(nsvint)             :: facestatev
      real(dp), dimension(ndime)              :: fcoord
      real(dp), dimension(ndofn)              :: forcei
      real(dp), dimension(ndofn)              :: forcej
      real(dp), dimension(ndofel)             :: forceVec
      real(dp), dimension(nprops)             :: props    
      real(dp)                                :: rbv
      real(dp), dimension(3,3)                :: RotMat   
      real(dp)                                :: sigN
      real(dp)                                :: sigM
      real(dp)                                :: sigL    
      real(dp)                                :: sigN0
      real(dp)                                :: sigM0
      real(dp)                                :: sigL0 
      real(dp)                                :: sine   
      real(dp)                                :: strainEN   
      real(dp), dimension(3)                  :: v          
      real(dp), dimension(3)                  :: vcross
      real(dp), dimension(3)                  :: vnr
      real(dp), dimension(3)                  :: vmr
      real(dp), dimension(3)                  :: vlr
      real(dp), dimension(3,3)                :: Vr
 

      ! Declare booleans
      logical(KIND=4) LargeDefFlag      
      logical(KIND=4) RateEffectFlag


      ! Assign values based on material
      SELECT CASE (fMat)
         CASE (0)
            E0 = props(2)
            alpha = props(3)
         CASE (1)
            E0 = props(2)
            alpha = props(3)
         CASE (2)
            E0 = props(25)
            alpha = props(26)
         CASE (3)
            E0 = props(48)
            alpha = props(49)
         CASE (4)
            E0 = props(71)
            alpha = props(72)
         CASE (5)
            E0 = props(94)
            alpha = props(95)
      END SELECT
   
      EN = E0;
      ET = alpha*E0;  

      inode = IDNINT(faceprop(1))
      jnode = IDNINT(faceprop(2))      


C======================================================================C
C     Get relative node numbers for abaqus element
C======================================================================C

      !inode = Kfindloc(nconnec,inode,4)
      !jnode = Kfindloc(nconnec,jnode,4)

      fcoord(:)=faceprop(3:5)
      area=faceprop(11)
      dn(:)=faceprop(12:14)
      dm(:)=faceprop(15:17)
      dl(:)=faceprop(18:20)


C======================================================================C
C     Length of the edge
C======================================================================C

      if (LargeDefFlag) then
         dx=coords(1,jnode)-coords(1,inode)
         dy=coords(2,jnode)-coords(2,inode)
         dz=coords(3,jnode)-coords(3,inode)
         alen=DSQRT(DABS(dx*dx+dy*dy+dz*dz))
         if (alen.lt.1.0d-10) then
            write(*,*)'Warning: alen in critical range; continuation may produce numerical issues.' 
            write(*,*)'Setting alen to 1.0E-13.'   
            alen = 1.0d-13
         endif
         vnr(1) = dx/alen
         vnr(2) = dy/alen
         vnr(3) = dz/alen
      else
         dx = coords0(1,jnode)-coords0(1,inode)
         dy = coords0(2,jnode)-coords0(2,inode)
         dz = coords0(3,jnode)-coords0(3,inode)
         alen = DSQRT(DABS(dx*dx+dy*dy+dz*dz))
         vnr(1) = dn(1)
         vnr(2) = dn(2)
         vnr(3) = dn(3)
      endif  


C======================================================================C
C     Updating the facet normal and tangential
C======================================================================C
      
      if (LargeDefFlag) then

         vcross(1) = dn(2)*vnr(3)-dn(3)*vnr(2)
         vcross(2) = dn(3)*vnr(1)-dn(1)*vnr(3)
         vcross(3) = dn(1)*vnr(2)-dn(2)*vnr(1)
         cosine = dn(1)*vnr(1)+dn(2)*vnr(2)+dn(3)*vnr(3)
         Vr(1,1) = 0
         Vr(1,2) = -vcross(3)
         Vr(1,3) = vcross(2) 
         Vr(2,1) = vcross(3)
         Vr(2,2) = 0
         Vr(2,3) = -vcross(1) 
         Vr(3,1) = -vcross(2)
         Vr(3,2) = vcross(1)
         Vr(3,3) = 0
         DO Icnt = 1,3
            DO Jcnt = 1,3
              Iden(Icnt,Jcnt) = 0
            END DO
            Iden(Icnt,Icnt) = 1
         END DO

         if ((cosine.lt.-0.99d0).and.(cosine.gt.-1.01d0)) then
            write(*,*)'Warning: Cosine in critical range; continuation may produce numerical issues.' 
            write(*,*)'Flipping normal.'   
            vmr(1) = -dm(1)
            vmr(2) = -dm(2)
            vmr(3) = -dm(3)
            vlr(1) = -dl(1)
            vlr(2) = -dl(2)
            vlr(3) = -dl(3)    
         else   
            RotMat = Iden + Vr + matmul(Vr,Vr)/(1+cosine)
            vmr(1) = DOT_PRODUCT(RotMat(1,:),dm(:))
            vmr(2) = DOT_PRODUCT(RotMat(2,:),dm(:))
            vmr(3) = DOT_PRODUCT(RotMat(3,:),dm(:))
            vlr(1) = DOT_PRODUCT(RotMat(1,:),dl(:))
            vlr(2) = DOT_PRODUCT(RotMat(2,:),dl(:))
            vlr(3) = DOT_PRODUCT(RotMat(3,:),dl(:)) 
         endif
      else
         vmr(1) = dm(1)
         vmr(2) = dm(2)
         vmr(3) = dm(3)
         vlr(1) = dl(1)
         vlr(2) = dl(2)
         vlr(3) = dl(3)
      endif   


C======================================================================C
C     Displacement increment of the nodes
C======================================================================C

      idof = (inode-1)*ndofn
      dui = du(idof+1:idof+ndofn,1)
      jdof = (jnode-1)*ndofn
      duj = du(jdof+1:jdof+ndofn,1) 
      
      call Cal_AMAT(ndime,ndofn,coords0(:,inode),fcoord,Ai)
      call Cal_AMAT(ndime,ndofn,coords0(:,jnode),fcoord,Aj)
          
      BjkN = MATMUL(vnr,Aj)/alen
      BikN = MATMUL(vnr,Ai)/alen
      
      BjkM = MATMUL(vmr,Aj)/alen
      BikM = MATMUL(vmr,Ai)/alen
     
      BjkL = MATMUL(vlr,Aj)/alen
      BikL = MATMUL(vlr,Ai)/alen


C======================================================================C
C     Calculate strain increment for a facet
C======================================================================C
           
      depsN=0.0d0; depsM=0.0d0; depsL=0.0d0     

      do idof=1,ndofn
         depsN=depsN+BjkN(idof)*duj(idof)-BikN(idof)*dui(idof)
         depsM=depsM+BjkM(idof)*duj(idof)-BikM(idof)*dui(idof)
         depsL=depsL+BjkL(idof)*duj(idof)-BikL(idof)*dui(idof)
      end do      

      epsN=facestatev(4)+depsN
      epsM=facestatev(5)+depsM 
      epsL=facestatev(6)+depsL


C======================================================================C
C     Call Material Subroutine (depending on material flag)
C======================================================================C

      sigN0=facestatev(1); sigM0=facestatev(2); sigL0=facestatev(3)

      SELECT CASE (fMat)
         CASE (0)
            call umatLDPM(23, props(1:23), rbv, facestatev, dtime, E0, 
     1                 alpha, EN, ET, aLen, epsV, depsN, depsM, depsL, RateEffectFlag)
         CASE (1)
            call umatLDPM(23, props(1:23), rbv, facestatev, dtime, E0, 
     1                 alpha, EN, ET, aLen, epsV, depsN, depsM, depsL, RateEffectFlag)
         CASE (2)
            call umatLDPM(23, props(24:46), rbv, facestatev, dtime, E0, 
     1                 alpha, EN, ET, aLen, epsV, depsN, depsM, depsL, RateEffectFlag)
         CASE (3)
            call umatLDPM(23, props(47:69), rbv, facestatev, dtime, E0, 
     1                 alpha, EN, ET, aLen, epsV, depsN, depsM, depsL, RateEffectFlag)
         CASE (4)
            call umatLDPM(23, props(70:92), rbv, facestatev, dtime, E0, 
     1                 alpha, EN, ET, aLen, epsV, depsN, depsM, depsL, RateEffectFlag)
         CASE (5)
            call umatLDPM(23, props(93:115), rbv, facestatev, dtime, E0, 
     1                 alpha, EN, ET, aLen, epsV, depsN, depsM, depsL, RateEffectFlag)
      END SELECT

      sigN=facestatev(1); sigM=facestatev(2); sigL=facestatev(3)


C======================================================================C
C     Strain energy
C======================================================================C

      strainEN=strainEN+0.5d0*alen*area*((sigN+sigN0)*depsN+
     1             (sigM+sigM0)*depsM+(sigL+sigL0)*depsL)    


C======================================================================C
C     Calculate Force Vector
C======================================================================C

      forcei = 0.0d0; forcej = 0.0d0

      do idof=1,ndofn
         forcei(idof)=forcei(idof)+sigN*BikN(idof)+
     1             sigM*BikM(idof)+sigL*BikL(idof)
         forcej(idof)=forcej(idof)+sigN*BjkN(idof)+
     1             sigM*BjkM(idof)+sigL*BjkL(idof)
      end do

      forcei=-forcei*area*alen
      forcej=forcej*area*alen
            
      do idof=1,ndofn
         icol=(inode-1)*ndofn
         jcol=(jnode-1)*ndofn
         forceVec(icol+idof)=forceVec(icol+idof)+forcei(idof)
         forceVec(jcol+idof)=forceVec(jcol+idof)+forcej(idof)
      end do    


C======================================================================C
C     Calculate Stifness Matrix
C======================================================================C      

      irow = (inode-1)*ndofn
      icol = (jnode-1)*ndofn   
      AreLE_N = area*alen*EN
      AreLE_T = area*alen*ET

      do idofn=1,ndofn          
         do jdofn=1,ndofn
         !Kii
         Elstif(irow+idofn,irow+jdofn) = Elstif(irow+idofn,irow+jdofn)+
     1            AreLE_N*BikN(idofn)*BikN(jdofn)+AreLE_T*
     2            (BikM(idofn)*BikM(jdofn)+BikL(idofn)*BikL(jdofn))

         !Kjj
         Elstif(icol+idofn,icol+jdofn) = Elstif(icol+idofn,icol+jdofn)+
     1            AreLE_N*BjkN(idofn)*BjkN(jdofn)+AreLE_T*
     2            (BjkM(idofn)*BjkM(jdofn)+BjkL(idofn)*BjkL(jdofn))

         !Kij
         Elstif(irow+idofn,icol+jdofn) = Elstif(irow+idofn,icol+jdofn)-
     1            AreLE_N*BikN(idofn)*BjkN(jdofn)-AreLE_T*
     2            (BikM(idofn)*BjkM(jdofn)+BikL(idofn)*BjkL(jdofn))

         !Kji
         Elstif(icol+idofn,irow+jdofn) = Elstif(icol+idofn,irow+jdofn)-
     1            AreLE_N*BjkN(idofn)*BikN(jdofn)-AreLE_T*
     2            (BjkM(idofn)*BikM(jdofn)+BjkL(idofn)*BikL(jdofn))

         end do
      end do      

      return                                      
      end subroutine  



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     A MATRIX CALCULATION IN LDPM                                     C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!


      subroutine Cal_AMAT(ndime,ndofn,coord,fcoord,Amat)
      

      implicit none


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare integers
      integer(KIND=4)                         :: ndime
      integer(KIND=4)                         :: ndofn
      

      ! Declare double precision floats     
      real(dp), dimension(ndime,ndofn)        :: Amat
      real(dp), dimension(ndime)              :: coord
      real(dp), dimension(ndime)              :: fcoord


      Amat = 0.0D0
      Amat(1,1) = 1.0d0
      Amat(1,5) = fcoord(3)-coord(3)
      Amat(1,6) = coord(2)-fcoord(2)
      Amat(2,2) = 1.0d0
      Amat(2,4) = -Amat(1,5)
      Amat(2,6) = fcoord(1)-coord(1)
      Amat(3,3) = 1.0d0
      Amat(3,4) = -Amat(1,6)
      Amat(3,5) = -Amat(2,6)

    
      return
      end subroutine



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     LDPM MASS MATRIX                                                 C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


C======================================================================C
C     Notes                                                            C
C======================================================================C
C     Xi and Xj are tet node coordinates related to eah facet          C
C     coordinates of three facet nodes:                                C
C     Xc node inside tet, XXa node on the edge, XXb node on tet face   C          
C     FacetVar(:,iface,jelem)=faceprop                                 C 
C======================================================================C 


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !      
      !================================================================!


      subroutine KMASSLDPM(nnode,ndofn,ndofel,ninpt,mcrd,density,
     1                     coords,nconnec,faceprop,amatrx) 


      implicit none  


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare integers
      integer(KIND=4)                         :: i
      integer(KIND=4)                         :: idofel
      integer(KIND=4)                         :: idofn
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: inode
      integer(KIND=4)                         :: ir
      integer(KIND=4)                         :: j
      integer(KIND=4)                         :: jdofel
      integer(KIND=4)                         :: jdofn
      integer(KIND=4)                         :: jnode
      integer(KIND=4)                         :: jr
      integer(KIND=4)                         :: Kfindloc
      integer(KIND=4), dimension(4)           :: nconnec
      integer(KIND=4)                         :: ndofel
      integer(KIND=4)                         :: nnode
      integer(KIND=4)                         :: ndofn
      integer(KIND=4)                         :: ninpt
      integer(KIND=4)                         :: mcrd


      ! Declare double precision floats     
      real(dp), dimension(ndofel,ndofel)      :: amatrx
      real(dp), dimension(mcrd,nnode)         :: coords
      real(dp)                                :: density
      real(dp), dimension(30,12)              :: faceprop
      real(dp), dimension(3)                  :: fcoord
      real(dp), dimension(ndofn,ndofn)        :: pMi
      real(dp), dimension(ndofn,ndofn)        :: pMj
      real(dp)                                :: vol
      real(dp)                                :: volhes
      real(dp), dimension(3)                  :: Xi
      real(dp), dimension(3)                  :: Xj
      real(dp), dimension(3)                  :: Xc
      real(dp), dimension(3)                  :: XXa
      real(dp), dimension(3)                  :: XXb
     

         
C======================================================================C
C     Initialize Mass Matrix
C======================================================================C   

      do idofel = 1, ndofel       
         do jdofel = 1, ndofel
            amatrx(jdofel,idofel) = 0.0d0
         end do
      end do      


C======================================================================C
C     Loop Over Facets
C======================================================================C      
      
      do iface=1,12           
     
         inode = IDNINT(faceprop(1,iface))
         jnode = IDNINT(faceprop(2,iface))
        
         ! Get relative node numbers for abaqus element
         !inode = Kfindloc(nconnec,inode,4)
         !jnode = Kfindloc(nconnec,jnode,4)
         fcoord(:) = faceprop(3:5,iface)
         Xc = faceprop(21:23,iface)
         XXa = faceprop(27:29,iface)
         XXb = faceprop(24:26,iface)


C======================================================================C
C     Compute Mass Matrix
C======================================================================C   

         Xi(1)=coords(1,inode)
         Xi(2)=coords(2,inode)
         Xi(3)=coords(3,inode)

         call CalcMI(ndofn,density,Xi,Xc,XXa,XXb,pMi,inode)

         volhes = volhes+vol

         Xj(1) = coords(1,jnode)
         Xj(2) = coords(2,jnode)
         Xj(3) = coords(3,jnode)

         call CalcMI(ndofn,density,Xj,Xc,XXa,XXb,pMj,jnode)        
  
         ir = (inode-1)*ndofn   
         jr = (jnode-1)*ndofn           
          
         do idofn=1,ndofn          
            do jdofn=1,ndofn
            !Mii
            amatrx(ir+idofn,ir+jdofn) = amatrx(ir+idofn,ir+jdofn) + pMi(idofn,jdofn)
            !Mjj
            amatrx(jr+idofn,jr+jdofn) = amatrx(jr+idofn,jr+jdofn) + pMj(idofn,jdofn)
            end do
         end do

      end do !END of LOOP OVER FACETS       
     

C======================================================================C
C     Evaluate a lumped mass matrix using the row sum method
C======================================================================C         
    
      do i = 1 , ndofel
         do j = 1 , ndofel
            if (i.ne.j) then
               amatrx(i,j) = 0.0d0;
            end if
         end do
      end do

      return
      end   



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     CALC MI SUBROUTINE                                               C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!


      subroutine CalcMI(ndofn,density,Xi,Xc,Xa,Xb,pMi,inode)


      implicit none


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare integers
      integer(KIND=4)                         :: inode
      integer(KIND=4)                         :: ndofn


      ! Declare double precision floats           
      real(dp), dimension(3,3)                :: coef 
      real(dp)                                :: density
      real(dp), dimension(ndofn,ndofn)        :: pMi
      real(dp)                                :: Sx
      real(dp)                                :: Sy
      real(dp)                                :: Sz
      real(dp)                                :: tetvol  
      real(dp)                                :: vIx
      real(dp)                                :: vIy
      real(dp)                                :: vIz
      real(dp)                                :: vIxy
      real(dp)                                :: vIxz
      real(dp)                                :: vIyz
      real(dp)                                :: vol
      real(dp)                                :: volutetra
      real(dp), dimension(3)                  :: Xi
      real(dp), dimension(3)                  :: Xa
      real(dp), dimension(3)                  :: Xb
      real(dp), dimension(3)                  :: Xc
      real(dp), dimension(3)                  :: Xg      
      real(dp), dimension(3)                  :: X
      real(dp), dimension(3)                  :: Y
      real(dp), dimension(3)                  :: Z    


      coef(1,1) = 0.1d0;  coef(1,2) = 0.05d0; coef(1,3) = 0.05d0      
      coef(2,1) = 0.05d0; coef(2,2) = 0.1d0;  coef(2,3) = 0.05d0      
      coef(3,1) = 0.05d0; coef(3,2) = 0.05d0; coef(3,3) = 0.1d0  

      tetvol = volutetra(Xi,Xa,Xb,Xc)

      vol = density*tetvol 

      Xg = (Xi+Xa+Xb+Xc)/4.0d0

      X(1) = Xa(1)-Xi(1); X(2)=Xb(1)-Xi(1); X(3)=Xc(1)-Xi(1)
      Y(1) = Xa(2)-Xi(2); Y(2)=Xb(2)-Xi(2); Y(3)=Xc(2)-Xi(2)
      Z(1) = Xa(3)-Xi(3); Z(2)=Xb(3)-Xi(3); Z(3)=Xc(3)-Xi(3)

      Sx = (Xg(1)-Xi(1))*vol
      Sy = (Xg(2)-Xi(2))*vol
      Sz = (Xg(3)-Xi(3))*vol

      vIx = DOT_PRODUCT(X,MATMUL(coef,X))*vol
      vIy = DOT_PRODUCT(Y,MATMUL(coef,Y))*vol
      vIz = DOT_PRODUCT(Z,MATMUL(coef,Z))*vol
      vIxy = DOT_PRODUCT(X,MATMUL(coef,Y))*vol
      vIxz = DOT_PRODUCT(X,MATMUL(coef,Z))*vol
      vIyz = DOT_PRODUCT(Y,MATMUL(coef,Z))*vol

      pMi(1,1)=vol; pMi(1,5)= Sz; pMi(1,6)=-Sy
      pMi(2,2)=vol; pMi(2,4)=-Sz; pMi(2,6)= Sx
      pMi(3,3)=vol; pMi(3,4)= Sy; pMi(3,5)=-Sx
      pMi(4,2)=-Sz; pMi(4,3)= Sy; pMi(4,4)=vIy+vIz; pMi(4,5)=-vIxy;   pMi(4,6)=-vIxz 
      pMi(5,1)= Sz; pMi(5,3)=-Sx; pMi(5,4)=-vIxy;   pMi(5,5)=vIx+vIz; pMi(5,6)=-vIyz 
      pMi(6,1)=-Sy; pMi(6,2)= Sx; pMi(6,4)=-vIxz;   pMi(6,5)=-vIyz;   pMi(6,6)=vIx+vIy 


      return
      end subroutine



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     CALC GI SUBROUTINE                                               C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!


      subroutine CalcGI(ndofn,density,Xi,Xc,Xa,Xb,pGi,inode)


      implicit none 


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare integers
      integer(KIND=4)                         :: inode
      integer(KIND=4)                         :: ndofn


      ! Declare double precision floats     
      real(dp)                                :: density
      real(dp), dimension(3,ndofn)            :: pGi
      real(dp)                                :: tetvol
      real(dp)                                :: vol
      real(dp)                                :: Sx
      real(dp)                                :: Sy
      real(dp)                                :: Sz
      real(dp)                                :: volutetra
      real(dp), dimension(3)                  :: Xi
      real(dp), dimension(3)                  :: Xa
      real(dp), dimension(3)                  :: Xb
      real(dp), dimension(3)                  :: Xc
      real(dp), dimension(3)                  :: Xg      
      real(dp), dimension(3)                  :: X
      real(dp), dimension(3)                  :: Y
      real(dp), dimension(3)                  :: Z  

      tetvol=volutetra(Xi,Xa,Xb,Xc)

      vol=tetvol*density

      Xg=(Xi+Xa+Xb+Xc)/4.0d0

      X(1)=Xa(1)-Xi(1); X(2)=Xb(1)-Xi(1); X(3)=Xc(1)-Xi(1)
      Y(1)=Xa(2)-Xi(2); Y(2)=Xb(2)-Xi(2); Y(3)=Xc(2)-Xi(2)
      Z(1)=Xa(3)-Xi(3); Z(2)=Xb(3)-Xi(3); Z(3)=Xc(3)-Xi(3)

      Sx=(Xg(1)-Xi(1))*vol; Sy=(Xg(2)-Xi(2))*vol; Sz=(Xg(3)-Xi(3))*vol

      pGi(1,1)=vol; pGi(1,5)= Sz; pGi(1,6)=-Sy
      pGi(2,2)=vol; pGi(2,4)=-Sz; pGi(2,6)= Sx
      pGi(3,3)=vol; pGi(3,4)= Sy; pGi(3,5)=-Sx;


      return
      end subroutine



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FUNCTION FOR VOLUME OF TET CALCULATION                           C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!      


      FUNCTION volutetra(p1,p2,p3,p4) result(tetvol) 


      implicit none


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare double precision floats     
      real(dp), dimension(3)                  :: p1
      real(dp), dimension(3)                  :: p2
      real(dp), dimension(3)                  :: p3
      real(dp), dimension(3)                  :: p4            
      real(dp)                                :: tetVol    


      tetvol =0.0d0
      tetvol = p2(1)*(p3(2)*p4(3)-p4(2)*p3(3))-p3(1)*(p2(2)*p4(3)-p4(2)*p2(3))+p4(1)*(p2(2)*p3(3)-p3(2)*p2(3))
      tetvol = tetvol-(p3(1)*(p4(2)*p1(3)-p1(2)*p4(3))-p4(1)*(p3(2)*p1(3)-p1(2)*p3(3))+p1(1)*(p3(2)*p4(3)-p4(2)*p3(3)))
      tetvol = tetvol+p4(1)*(p1(2)*p2(3)-p2(2)*p1(3))-p1(1)*(p4(2)*p2(3)-p2(2)*p4(3))+p2(1)*(p4(2)*p1(3)-p1(2)*p4(3))
      tetvol = tetvol-(p1(1)*(p2(2)*p3(3)-p3(2)*p2(3))-p2(1)*(p1(2)*p3(3)-p3(2)*p1(3))+p3(1)*(p1(2)*p2(3)-p2(2)*p1(3)))
      tetvol = DABS(tetvol)/6.0d0


      END FUNCTION volutetra



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE FOR BODY FORCE LDPM                                   C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


C Xi and Xj are tet node coordinates related to eah facet
C coordinates of three facet nodes:
C Xc node inside tet, XXa node on the edge, XXb node on tet face            
C FacetVar(:,iface,jelem)=faceprop    


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!  


      subroutine BODYF_LDPM(nnode,ndofn,ndofel,ninpt,mcrd,density,
     1                     coords,nconnec,faceprop,accelVec,BodyL) 


      implicit none 


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare integers
      integer(KIND=4)                         :: idofel
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: inode
      integer(KIND=4)                         :: ir
      integer(KIND=4)                         :: jnode
      integer(KIND=4)                         :: jr
      integer(KIND=4)                         :: Kfindloc      
      integer(KIND=4), dimension(4)           :: nconnec
      integer(KIND=4)                         :: nnode
      integer(KIND=4)                         :: ndofn
      integer(KIND=4)                         :: ndofel
      integer(KIND=4)                         :: ninpt
      integer(KIND=4)                         :: mcrd


      ! Declare double precision floats     
      real(dp), dimension(1,3)                :: accelVec
      real(dp), dimension(ndofel)             :: BodyL
      real(dp), dimension(1,ndofn)            :: BforceI
      real(dp), dimension(1,ndofn)            :: BforceJ
      real(dp), dimension(mcrd,nnode)         :: coords
      real(dp)                                :: density
      real(dp), dimension(30,12)              :: faceprop
      real(dp), dimension(3)                  :: fcoord
      real(dp), dimension(3,ndofn)            :: pGi
      real(dp), dimension(3,ndofn)            :: pGj
      real(dp), dimension(3)                  :: Xi
      real(dp), dimension(3)                  :: Xj
      real(dp), dimension(3)                  :: Xc
      real(dp), dimension(3)                  :: XXa
      real(dp), dimension(3)                  :: XXb


C======================================================================C
C     Initialize the BodyL
C======================================================================C   

      do idofel = 1,ndofel
         BodyL(idofel) = 0.0d0
      end do

C======================================================================C
C     Loop Over Facet
C======================================================================C               

      do iface=1,12       
    
         inode = IDNINT(faceprop(1,iface))
         jnode = IDNINT(faceprop(2,iface))
         
         ! Get relative node numbers for abaqus element
         !inode = Kfindloc(nconnec,inode,4)
         !jnode = Kfindloc(nconnec,jnode,4)  

         Xc  = faceprop(21:23,iface)
         XXa = faceprop(27:29,iface)
         XXb = faceprop(24:26,iface)

         Xi(1) = coords(1,inode)
         Xi(2) = coords(2,inode)
         Xi(3) = coords(3,inode)

         call CalcGI(ndofn,density,Xi,Xc,XXa,XXb,pGi,inode)      

         BforceI = MATMUL(accelVec,pGi)
      
         Xj(1) = coords(1,jnode)
         Xj(2) = coords(2,jnode)
         Xj(3) = coords(3,jnode)

         call CalcGI(ndofn,density,Xj,Xc,XXa,XXb,pGj,jnode)    

         BforceJ = MATMUL(accelVec,pGj)
       
         ir = (inode-1)*ndofn
         jr = (jnode-1)*ndofn
C
            
         BodyL(ir+1:ir+ndofn) = BodyL(ir+1:ir+ndofn) + BforceI(1,:)
         BodyL(jr+1:jr+ndofn) = BodyL(jr+1:jr+ndofn) + BforceJ(1,:)    
C
           
      end do !END of LOOP OVER FACETS    
C
      return
      end   



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE ELEMENT CONNECTIVITY                                  C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!  


      subroutine get_elem_con 


      use ModuleLDPM
      implicit none 


      ! Declare integers
      integer(KIND=4)                         :: io
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer(KIND=4)                         :: nelcnt


      ! Declare characters
      character(256)                          :: filename
      character(256)                          :: input
      character(256)                          :: jobname
      character(256)                          :: outdir
      character(256)                          :: separator
      character(1)                            :: path
      character(1)                            :: path_separator


C======================================================================C
C     Read element connectivity for facet information                  C
C     (https://www.eng-tips.com/viewthread.cfm?qid=269244)             C
C======================================================================C                   

              
      !Open INP file for reading
      call vgetjobname(jobname,LENJOBNAME)
      call vgetoutdir(outdir,LENOUTDIR)

      separator = outdir(1:1)

      if (separator.eq.'/') then
         path_separator = '/'
      else
         path_separator = '\'
      endif


      filename=outdir(1:LENOUTDIR)//path_separator//
     *jobname(1:LENJOBNAME)//'.inp' 

      open(unit=17,file=filename(1:LENOUTDIR+                          
     *LENJOBNAME+5),status='unknown')
     
      nelcnt = 0

  
      !Skip down to *Element in INP File
      do while (index(input,'*Element')==0)              
         read(17,'(a)',ERR=200) input              
      end do
                              
      !Read in Element Nodal Connectivity      
      do while(.true.)
         read(17,*,ERR=200)input                                                
         if(index(input,'*')==0)then
            backspace(17)
            nelcnt=nelcnt+1      
            read(17,*,ERR=200) elecon(nelcnt,1),elecon(nelcnt,2),elecon(nelcnt,3),elecon(nelcnt,4),elecon(nelcnt,5)
            elecon(nelcnt,1)=nelcnt
         else
            exit
         endif
      end do  
          
200   close(17)!closes file    


      return
      end subroutine get_elem_con



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE FACET INFORMATION                                     C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!  


      subroutine get_facet_info_old 


      use ModuleLDPM
      implicit none  


      ! Declare integers
      integer(KIND=4)                         :: io
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: iwhere
      integer(KIND=4)                         :: itet
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer(KIND=4)                         :: nelcnt


      ! Declare characters
      character(256)                          :: filename
      character(256)                          :: input
      character(256)                          :: jobname
      character(256)                          :: line
      character(256)                          :: outdir
      character(1)                            :: path
      character(1)                            :: path_separator
      character(256)                          :: separator

                
C======================================================================C
C     Read facet information from facet file                           C
C======================================================================C
   
      !Open facedata file for reading.
      call vgetjobname(jobname,lenjobname)
      call vgetoutdir(outdir,lenoutdir)


      CALL getcwd(path)

      separator = outdir(1:1)
      if (separator.eq.'/') then
         path_separator = '/'
      else
         path_separator = '\'
      endif


      filename=outdir(1:lenoutdir)//path_separator//
     *jobname(1:lenjobname)//'_facet.dat' 
      open(unit=18,file=filename(1:lenoutdir+                          
     *lenjobname+11),status='old')
c            
      nelcnt=0
10    read (18, '(A)', end = 20) line
      iwhere = index (line, 'Tet index:')
      if (iwhere.ne.0) then
          read(line(11:20),'(I8)') itet !int(line(11:20))         
          nelcnt=nelcnt+1      
  
          do iface=1,12
              read(18,*,ERR=200) FacetVar(1:30,iface,nelcnt)          
          end do       
      endif
      goto 10
20    continue      
200   close(18)!closes file 
   
      return
      end subroutine get_facet_info_old



      
   
      subroutine get_facet_info 


      use ModuleLDPM
      implicit none  


      ! Declare integers
      integer(KIND=4)                         :: io
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: iwhere
      integer(KIND=4)                         :: itet
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer(KIND=4)                         :: nelcnt
      integer(KIND=4)                         :: mf
      integer(KIND=4)                         :: ipos
      integer(KIND=4)                         :: iver
      integer(KIND=4)                         :: Kfindloc_new
      
      
      ! Declare floats
      real(dp)                                :: vol
      real(dp)                                :: area
      !real(dp), dimension(30,12,2348)         :: facetvar
      real(dp), dimension(3,36*nomaxel)       :: vertices
      real(dp), dimension(3)                  :: verts
      real(dp), dimension(3)                  :: center 
      real(dp), dimension(3)                  :: n
      real(dp), dimension(3)                  :: m
      real(dp), dimension(3)                  :: l

      ! Declare characters
      character(256)                          :: filename
      character(256)                          :: input
      character(256)                          :: jobname
      character(256)                          :: line
      character(256)                          :: outdir
      character(1)                            :: path
      character(1)                            :: path_separator
      character(256)                          :: separator

                
C======================================================================C
C     Read facet information from facet file                           C
C======================================================================C
   
      !Open facedata file for reading.
      call vgetjobname(jobname,lenjobname)
      call vgetoutdir(outdir,lenoutdir)


      CALL getcwd(path)

      separator = outdir(1:1)
      if (separator.eq.'/') then
         path_separator = '/'
      else
         path_separator = '\'
      endif

      
      filename=outdir(1:lenoutdir)//path_separator//
     *jobname(1:lenjobname)//'-data-facets.dat' 
      open(unit=18,file=filename(1:lenoutdir+                          
     *lenjobname+17),status='old')
      
      filename=outdir(1:lenoutdir)//path_separator//
     *jobname(1:lenjobname)//'-data-facetsVertices.dat' 
      open(unit=19,file=filename(1:lenoutdir+                          
     *lenjobname+25),status='old')
      
c      
c     Read vertices data
c
      iver=0
10    read (19, '(A)', end = 20) line
c     
      if (index (line, '//').eq.0) then
          iver=iver+1          
          read(line,*,ERR=200) vertices(1:3,iver)
      endif
      goto 10
20    continue
200   close(19)!closes file 
     
     
c      
c     Read facet data and combine with vertices data
c            
      iface=0
11    read (18, '(A)', end = 21) line
c     
      if (index (line, '//').eq.0) then
          iface=iface+1          
          read(line,*,ERR=201) nelcnt, verts, vol, area, center, 
     *                                                n, m, l, mf 
           
           FacetVar(1,iface,nelcnt+1) =Kfindloc_new(iface,1)
           FacetVar(2,iface,nelcnt+1) =Kfindloc_new(iface,2)
           FacetVar(3:5,iface,nelcnt+1) =center
           FacetVar(10,iface,nelcnt+1) = vol
           FacetVar(11,iface,nelcnt+1) = area
           FacetVar(12:14,iface,nelcnt+1)=n
           FacetVar(15:17,iface,nelcnt+1)=m
           FacetVar(18:20,iface,nelcnt+1)=l
           FacetVar(21:23,iface,nelcnt+1)=
     *                vertices(:,verts(1)+1)
           FacetVar(27:29,iface,nelcnt+1)=
     *                vertices(:,verts(2)+1)
           FacetVar(24:26,iface,nelcnt+1)=
     *                vertices(:,verts(3)+1)
           FacetVar(30,iface,nelcnt+1)=mf		   
          if (iface.eq.12) iface=0               
      endif
      goto 11
21    continue      
201   close(18)!closes file 
      
      return
      end subroutine get_facet_info
      
      
      
      
      
      function Kfindloc_new(iface,j) result(ipos)


      implicit none


      ! Declare integers      
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: j
      integer(KIND=4)                         :: ipos
      integer(KIND=4), dimension(2)           :: mpos
      
      


C       select case(iface)
C       case(1)
C           mpos=(/1, 2/) 
C       case(2)   
C           mpos=(/1, 2/)
C       case(3)
C           mpos=(/1, 3/) 
C       case(4)   
C           mpos=(/1, 3/)
C       case(5)
C           mpos=(/2, 3/) 
C       case(6)   
C           mpos=(/2, 3/)
C       case(7)
C           mpos=(/2, 4/) 
C       case(8)   
C           mpos=(/2, 4/)
C       case(9)
C           mpos=(/3, 4/) 
C       case(10)   
C           mpos=(/3, 4/)
C       case(11)
C           mpos=(/1, 4/) 
C       case(12)   
C           mpos=(/1, 4/) 
C       end select
	  
	  
	   select case(iface)
      case(1)
          mpos=(/1, 2/) 
      case(2)   
          mpos=(/1, 2/)
      case(3)
          mpos=(/1, 3/) 
      case(4)   
          mpos=(/1, 3/)
	  case(5)
          mpos=(/1, 4/) 
      case(6)   
          mpos=(/1, 4/) 
      case(7)
          mpos=(/2, 3/) 
      case(8)   
          mpos=(/2, 3/)
      case(9)
          mpos=(/2, 4/) 
      case(10)   
          mpos=(/2, 4/)
      case(11)
          mpos=(/3, 4/) 
      case(12)   
          mpos=(/3, 4/)     
      end select
      
      ipos=mpos(j)     
          
      end function
      
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE FACET INFORMATION                                     C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================!  


      subroutine get_svars_config


      use ModuleLDPM
      implicit none  
                

      ! Declare integers
      integer(KIND=4)                         :: io
      integer(KIND=4)                         :: iface
      integer(KIND=4)                         :: iwhere
      integer(KIND=4)                         :: itet
      integer(KIND=4)                         :: LENJOBNAME
      integer(KIND=4)                         :: LENOUTDIR
      integer(KIND=4)                         :: nelcnt
      integer(KIND=4)                         :: num_svars


      ! Declare characters
      character(256)                          :: filename
      character(256)                          :: input
      character(256)                          :: jobname
      character(256)                          :: line
      character(256)                          :: outdir
      character(1)                            :: path
      character(1)                            :: path_separator
      character(256)                          :: separator


C======================================================================C
C     Read svars selection from local file                             C
C======================================================================C

      CALL getcwd(path)
      call vgetoutdir(outdir,LENOUTDIR)

      separator=outdir(1:1)
      if (separator.eq.'/') then
         path_separator = '/'
      else
         path_separator = '\'
      endif
      

      filename=trim(outdir(1:LENOUTDIR)//path_separator//'svars_config.dat')
      open(unit=138,file=filename,status='old')
        
10      read (138, '(A)', end = 20) line
      iwhere = index (line, 'Number of svars:')
      if (iwhere.ne.0) then
          read(line(17:20),'(I3)') num_svars
          allocate(svars_selection(num_svars))
          if (num_svars.gt.0) then
            read(138,*,ERR=200) svars_selection
          else
            read(138,*,ERR=200) line
          endif
          read(138,*,ERR=200) svars_freq
          read(138,*,ERR=200) svars_nodal_stress
          read(138,*,ERR=200) print_binary
      endif
      goto 10
20    continue       
200   close(138) !closes file 
      return
      end subroutine get_svars_config



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FUNCTION KFINDLOC                                                C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================! 


      function Kfindloc(Iarray,ival,nel) result(ipos)


      implicit none


      ! Declare integers
      integer(KIND=4)                         :: nel
      integer(KIND=4)                         :: i
      integer(KIND=4)                         :: ipos
      integer(KIND=4)                         :: ival
      integer(KIND=4), dimension(nel)         :: Iarray


      do i=1,nel
          if (Iarray(i).eq.ival) then
              ipos=i
              return
          endif
      end do
      end function



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE LDPM CONSTITUTIVE EQUATIONS                           C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================! 

     
      subroutine umatLDPM(nprops, props, rbv,stv, dt, E0, alpha, EN,ET, 
     1                 aLength,epsV,depsN, depsM, depsL, RateEffectFlag)

     
c     LDPM MATERIAL MODEL
c      
c     LDPM Facet Constitutive Law - Cusatis/Rezakhani Oct 2019
c
c     stv state variable
c     stv(1)  Normal N stress
c     stv(2)  Shear M stress
c     stv(3)  Shear L stress
c     stv(4)  Normal N strain
c     stv(5)  Shear M strain
c     stv(6)  Shear L strain
c     stv(7)  Max normal N strain
c     stv(8)  Max shear T strain
c     stv(9)  Tensile strength
c     stv(10) Post-peak slope in tension
c     stv(11) Shear L crack opening
c     stv(12) Minimum Normal Strain
c     stv(13) Normal N crack opening
C     stv(14) Shear M crack opening
C     stv(15) Total crack opening
C     stv(16) Volumetric Strain
C     stv(17) Dissipated energy density rate
C     stv(18) Dissipated energy density
C     stv(19) Dissipated energy density rate in tension
C     stv(20) Dissipated energy density in tension
 

      implicit none  


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)
     

      ! Declare parameters
      real(dp), parameter :: PI = 3.141592653589793d0
      real(dp), parameter :: eps = 1.0d-15


      ! Declare integers
      integer(KIND=4)                         :: icod  
      integer(KIND=4)                         :: nprops  
      integer(KIND=4),  dimension(25)         :: VectorSVIndex  


      ! Declare double precision floats
      real(dp)                                :: aKc
      real(dp)                                :: aKc1
      real(dp)                                :: aKcc
      real(dp)                                :: aKt
      real(dp)                                :: aLength
      real(dp)                                :: alpha
      real(dp)                                :: ateta
      real(dp)                                :: beta  
      real(dp)                                :: bound_f
      real(dp)                                :: bound_f_tmp
      real(dp)                                :: bound_N
      real(dp)                                :: bound_N_tmp
      real(dp)                                :: bound_T
      real(dp)                                :: chi
      real(dp)                                :: chLen 
      real(dp)                                :: csi
      real(dp)                                :: csiO
      real(dp)                                :: csiMax
      real(dp)                                :: cteta
      real(dp)                                :: cteta2
      real(dp)                                :: c0
      real(dp)                                :: c1
      real(dp)                                :: Dcsi
      real(dp)                                :: DDE
      real(dp)                                :: DensRatio
      real(dp)                                :: depsL
      real(dp)                                :: depsM
      real(dp)                                :: depsN
      real(dp)                                :: depsIL
      real(dp)                                :: depsIM
      real(dp)                                :: depsIN
      real(dp)                                :: dk1
      real(dp)                                :: dk2
      real(dp)                                :: dk3
      real(dp)                                :: dmu
      real(dp)                                :: dns  
      real(dp)                                :: dt
      real(dp)                                :: dtDyn
      real(dp)                                :: ec
      real(dp)                                :: ec0
      real(dp)                                :: E0
      real(dp)                                :: EN
      real(dp)                                :: ENc
      real(dp)                                :: Eu
      real(dp)                                :: ep0
      real(dp)                                :: epsEq
      real(dp)                                :: epsD
      real(dp)                                :: epsN
      real(dp)                                :: epsT
      real(dp)                                :: epsV
      real(dp)                                :: epsV0
      real(dp)                                :: ET  
      real(dp)                                :: ETc
      real(dp)                                :: fc  
      real(dp)                                :: Fdyn
      real(dp)                                :: fmu_0  
      real(dp)                                :: fmu_inf
      real(dp)                                :: fr
      real(dp)                                :: fs
      real(dp)                                :: ft  
      real(dp)                                :: H
      real(dp)                                :: Hr  
      real(dp)                                :: Hs
      real(dp), dimension(45)                 :: hxv
      real(dp)                                :: phiD
      real(dp), dimension(nprops)             :: Props
      real(dp)                                :: ptf
      real(dp)                                :: rat
      real(dp)                                :: rat2
      real(dp), dimension(8)                  :: rbv
      real(dp)                                :: RinHardMod 
      real(dp)                                :: s0
      real(dp)                                :: sen_c  
      real(dp)                                :: sqalpha
      real(dp)                                :: sf0  
      real(dp)                                :: sigL0
      real(dp)                                :: sigM0
      real(dp)                                :: sigN0
      real(dp)                                :: sigL
      real(dp)                                :: sigM
      real(dp)                                :: sigN
      real(dp)                                :: sigT
      real(dp)                                :: sigMe
      real(dp)                                :: sigLe
      real(dp)                                :: sigTe
      real(dp)                                :: sc
      real(dp)                                :: sc0
      real(dp)                                :: ss
      real(dp)                                :: ssD
      real(dp)                                :: ssDdam
      real(dp)                                :: st
      real(dp)                                :: steta
      real(dp)                                :: steta2
      real(dp)                                :: Stot
      real(dp)                                :: Stot0
      real(dp), dimension(25)                 :: stv
      real(dp)                                :: teta
      real(dp)                                :: tmp
      real(dp)                                :: tmp1
      real(dp)                                :: tmp2
      real(dp)                                :: tsrn_e  
      real(dp)                                :: unkt
      real(dp)                                :: unks
      real(dp)                                :: unkc
      real(dp)                                :: x


      ! Declare booleans
      logical(KIND=4) EAF
      logical(KIND=4) RateEffectFlag


      !Material parameters
      dns           = Props(1)     !Density
      E0            = Props(2)     !Young Modulus ! GC: This should be EN=E0
      alpha         = Props(3)     !Poisson ratio ! GC: This should be alpha
      ft            = Props(4)     !Tensile Strength
      chLen         = Props(5)     !Tensile characteristic length
      fr            = Props(6)     !Shear strength ratio
      sen_c         = Props(7)     !Softening exponent
      fc            = Props(8)     !Compressive Yield Strength
      RinHardMod    = Props(9)     !Initial hardening modulus ratio
      tsrn_e        = Props(10)    !Transitional Strain ratio
      dk1           = Props(11)    !Deviatoric strain threshold ratio
      dk2           = Props(12)    !Deviatoric damage parameter
      fmu_0         = Props(13)    !Initial friction
      fmu_inf       = Props(14)    !Asymptotic friction
      sf0           = Props(15)    !Transitional stress 
      DensRatio     = Props(16)    !Densification ratio  
      beta          = Props(17)    !Volumetric deviatoric coupling
      unkt          = Props(18)    !Tensile unloading parameter
      unks          = Props(19)    !Shear unloading parameter
      unkc          = Props(20)    !Compressive unloading parameter
      Hr            = Props(21)    !Shear softening modulus ratio
      dk3           = Props(22)    !FinalHardeningModulusRatio
      EAF           = Props(23)    !ElasticAnalysisFlag
       
      VectorSVIndex=-1
      
      !Facet failure 
      stv(16) = epsV
      if (stv(16)>0.2d0) then
          stv(1) = 0.0D0
          stv(2) = 0.0d0
          stv(3) = 0.0d0
          stv(4) = stv(4)+depsN
          if(stv(4).ge.0.0d0) then
              stv(13) = stv(13)+ depsN*aLength   ! wN
              stv(14) = stv(14)+ depsM*aLength   ! wM
              stv(11) = stv(11)+ depsL*aLength   ! wL
              stv(15) = DSQRT(DABS(stv(13)*stv(13) + stv(14)*stv(14) + 
     1                   stv(11)*stv(11)))               ! wtotal
          endif
          return
      endif

      sigN = 0.0d0
      sigM = 0.0d0
      sigL = 0.0d0

      icod = 0
      if (EAF) icod = 1; ! Elastic Behavior
      
      ! Parameters Initializing    
      sqalpha=DSQRT(alpha)   
      ENc = EN
      ETc = ET    


C======================================================================C
C     Strain Increment                                                 C
C======================================================================C

      !Old Stresses 
      sigN0 = stv(1)
      sigM0 = stv(2)
      sigL0 = stv(3)
      if (icod.eq.1) then     ! Elastic behavior
            ! Compute New Stress Vector               
            stv(1) = stv(1) + EN * DepsN
            stv(2) = stv(2) + ET * DepsM
            stv(3) = stv(3) + ET * DepsL
            stv(4) = stv(4) + DepsN
            stv(5) = stv(5) + DepsM 
            stv(6) = stv(6) + DepsL
      else 
        fs = fr*ft
        ss = fs
        sc = fc       

        ! Dissipation Energy Objectivity during Fracture Processes
        if ((stv(9).eq.0.0d0).or.(stv(10).eq.0.0d0)) then
            chi = 0.99d0

            !tensile strength and post-peak slope
            if (aLength < chi * chLen) then
                st = ft
                if (aLength.lt.1.0d-10) then
                   write(*,*)'Warning: aLength in critical range; continuation may produce numerical issues.'     
                endif
                if (aLength.eq.chLen) then
                  write(*,*)'Warning: aLength equals chLen; continuation may produce numerical issues.'
                endif
                aKt = 2.0d0 * E0 / (-1.0d0 + chLen / aLength);
            else                 
                Write(*,*)'Warning: Edge length > tensile characteristic length'
                aKt = 2.0d0 * E0 / (-1.0d0 + 1.0d0 / chi);
                st = ft * DSQRT(DABS(chi * chLen / aLength));
            endif
            stv(9) = st
            stv(10) = aKt
         else 
            st = stv(9)
            aKt = stv(10)
         endif

        ! Shear strength and post-peak slope
        Hs = Hr * E0
        ss = fs

        ! Old effective Strain
        csiO = DSQRT(DABS(stv(4) * stv(4) + alpha * 
     1              (stv(5) * stv(5) + stv(6) * stv(6))))


C======================================================================C
C     New Strains                                                      C
C======================================================================C


        stv(4) = stv(4) + DepsN
        stv(5) = stv(5) + DepsM
        stv(6) = stv(6) + DepsL
        epsN = stv(4)  ! Normal Strain
        epsT = DSQRT(DABS(stv(5) * stv(5) + stv(6) * stv(6)))  ! Total Shear Strain
        csi  = DSQRT(DABS(epsN * epsN + alpha * epsT * epsT)) ! New Effective Strain
        if ((epsN.lt.EPSILON(epsN)).and.(epsN.gt.-EPSILON(epsN))) epsN = 0.0d0
        if ((epsV.lt.EPSILON(epsV)).and.(epsV.gt.-EPSILON(epsV))) epsV = 0.0d0
        if (((epsV.lt.EPSILON(epsV)).and.(epsV.gt.-EPSILON(epsV))).and.((epsN.lt.EPSILON(epsN)).and.(epsN.gt.-EPSILON(epsN)))) then
           epsD = 0.0d0
        else
           epsD = epsN - epsV
        endif

        ! Effective Strain Increment
        if ((csi.lt.EPSILON(csi)).and.(csi.gt.-EPSILON(csi))) csi = 0.0d0
        if ((csiO.lt.EPSILON(csiO)).and.(csiO.gt.-EPSILON(csiO))) csiO = 0.0D0
        if (((csi.lt.EPSILON(csi)).and.(csi.gt.-EPSILON(csi))).and.((csiO.lt.EPSILON(csiO)).and.(csiO.gt.-EPSILON(csiO)))) then
           Dcsi = 0.0d0
        else
           Dcsi = csi - csiO
        endif

        ! Coupling Variable
        teta = 0.0d0
        tmp = epsT * sqalpha  
        

        if ((epsN.lt.EPSILON(epsN)).and.(epsN.gt.-EPSILON(epsN))) then
            teta = DATAN(0.0d0)
        else if ((tmp.lt.(-EPSILON(tmp))).or.(tmp.gt.(EPSILON(tmp)))) then
            teta = DATAN(epsN / tmp) 
        else if (epsN > 0) then
            teta = PI / 2.0d0
        else if  (epsN < 0) then
            teta = -PI / 2.0d0
        endif

        ateta = teta * 2.0d0 / PI


        tmp = stv(1)*stv(1) + (stv(2) * stv(2) + stv(3) * stv(3))/alpha
        Stot0 = DSQRT(DABS(tmp))  !effective old stress

        if (epsN < stv(12)) stv(12) = epsN  ! Min Normal Strain
        !  epsNmin = stv(12) GC: check this


        !  Cohesive Damage Irreversibility
        if (epsN > stv(7)) stv(7) = epsN  ! Max Normal Strain
        if (epsT > stv(8)) stv(8) = epsT  ! Max Shear Strain
        csiMax = DSQRT(DABS(stv(7) * stv(7) + alpha * stv(8) * stv(8)))  ! Max Effective Strain


        ! Strain Rate Effect
        if (RateEffectFlag) then
            dtDyn = ptf * dt
            if ((dtDyn.lt.EPSILON(dtDyn)).and.(dtDyn.gt.-EPSILON(dtDyn))) then
                Fdyn = 1.0d0
            else 
                 x = aLength * DABS(Dcsi) / (c0 * dtDyn)
                 Fdyn = 1.0d0 + c1 * DLOG(DABS(x + DSQRT(DABS(1.0d0 + x * x))))
            endif
         else
            Fdyn = 1.0d0
         endif                  
        
C======================================================================C
C     Fracturing                                                       C
C======================================================================C


            if (epsN > 0.0d0) then ! Fracturing Behavior                

                steta  = DSIN(teta)
                steta2 = steta * steta
                cteta  = DCOS(teta)
                cteta2 = cteta * cteta

                if (st.lt.1.0d-10) then 
                  write(*,*)'Warning: st in critical range; continuation may produce numerical issues.',st         
                endif
                rat = ss / sqalpha / st
                rat2 = rat*rat !GC double check           
                if (cteta2 > eps) then
                    s0 = st*(-steta + DSQRT(DABS(steta2+4.0d0*cteta2/rat2)))
     1                      / (2.0d0 * cteta2 / rat2)
                else
                    s0 = st
                endif

                ep0 = s0 / E0
                if ((ateta.lt.0.0d0).and.(sen_c.lt.1.0d0)) then
                  write(*,*)'Warning: ateta and sec_c in ranges to produce exponent failure',ateta,sen_c
                endif
                H = Hs / alpha + (aKt - Hs / alpha) * ateta**sen_c
                tmp = csiMax - Fdyn * ep0       
              
                if (tmp.lt.-EPSILON(tmp)) then
                    bound_f_tmp = Fdyn * s0
                else
                    bound_f_tmp = Fdyn * s0 * DEXP(-H * tmp / s0)
                endif

                tmp1 = unkt * (csiMax - bound_f_tmp / E0)
                if (csi.gt.tmp1) then
                    bound_f = bound_f_tmp
                else
                    bound_f = 0.0d0
                endif                
     
            ! New effective Stress
                Stot = MIN(bound_f , MAX(0.0d0, Stot0 + E0 * Dcsi))
            ! Facet Failure Definition


            ! New Stresses
                if ((csi.gt.(-EPSILON(csi))).and.(csi.lt.(EPSILON(csi)))) then
                    sigN = 0.0d0
                    sigM = 0.0d0
                    sigL = 0.0d0
                else 
                    sigN = Stot * stv(4) / csi
                    sigM = Stot * alpha * stv(5) / csi
                    sigL = Stot * alpha * stv(6) / csi
                endif                
   
            else  ! Pore Collapse / Frictional Behavior                
 
  
C======================================================================C
C     Normal Stress                                                    C
C======================================================================C


                ec = sc / E0
                ec0 = tsrn_e * ec                
    
                if (sigN0.lt.-sc) then
                    ENc = DensRatio * E0;
                    !ETc = alpha * Eu
                endif
                phiD  = 1.0d0
                epsV0 = 0.1d0*ec
                aKc   = RinHardMod*E0
                aKc1  = E0*dk3
                if (epsV.lt.0.0d0) then
                    tmp = dk2*(-DABS(epsD) / (epsV-epsV0) - dk1)
                else 
                    tmp = dk2*(DABS(epsD) / (epsV0) - dk1)
                endif
                if (tmp >= 0.0d0 )  phiD = 1.0d0 / ( 1.0d0 + tmp )        
                !
                aKcc =  (aKc-aKc1) * phiD + aKc1
                !
                epsEq = epsV + beta * epsD
                tmp = epsEq + ec
                tmp1 = epsEq + ec0
                !
                if ((tmp < 0.0d0).and.(tmp1 > 0.0d0)) then
                    bound_N_tmp = -sc + aKcc * tmp
                else if (tmp1 < 0.0d0 ) then
                    sc0 = sc + aKcc * ( ec0 - ec )
                    bound_N_tmp = -sc0 * DEXP( -aKcc * tmp1 / sc0 )
                else 
                    bound_N_tmp = -sc
                endif                
    
          !    tmp1 = unkc * (- epsNmin + bound_N_tmp / ENc)
          !    if (- epsN > tmp1)
               bound_N = bound_N_tmp;
          !  else
          !    bound_N = 0.0;
          !
                sigN = MAX(bound_N, MIN(0.0d0, sigN0 + ENc * DepsN))                
           
                ! Damage of the Cohesive Component 
                ssD  = Fdyn * ss;
                tmp2 = csiMax / sqalpha - ssD / ET
                ssDdam = 0.0d0
                if (tmp2 < 0.0d0) then
                    ssDdam = ssD
                else
                    ssDdam = ssD * dexp( -Hs * tmp2 / ss )
                endif
    
                ! Shear Boundary
                dmu = fmu_0 - fmu_inf                   
                bound_T = ssDdam + dmu * sf0 - fmu_inf * sigN -
     1                      dmu * sf0 * DEXP(sigN / sf0)                

                ! New Tangential Stresses
                sigMe = sigM0 + ETc * DepsM
                sigLe = sigL0 + ETc * DepsL
                sigTe = DSQRT(DABS(sigMe * sigMe + sigLe * sigLe))
                sigT  = MIN(bound_T , MAX(0.0d0, sigTe))
                !
                if (sigTe.gt.EPSILON(sigTe)) then
                    sigM = sigT * sigMe / sigTe
                    sigL = sigT * sigLe / sigTe
                else 
                    sigM = 0.0d0
                    sigL = 0.0d0
                endif                

        endif  
 
        stv(1) = sigN
        stv(2) = sigM
        stv(3) = sigL     
      
      endif      ! **********  END OF ICOD IF STATEMENT


C======================================================================C
C     Increment of Inelastic Behavior                                  C
C======================================================================C


      DepsIN = DepsN - (stv(1) - sigN0) / ENc
      DepsIM = DepsM - (stv(2) - sigM0) / ETc
      DepsIL = DepsL - (stv(3) - sigL0) / ETc

      ! increment of dissipated energy density
      DDE = 0.5d0 * (stv(1) + sigN0) * DepsIN
      DDE = DDE + 0.5d0 * (stv(2) + sigM0) * DepsIM
      DDE = DDE + 0.5d0 * (stv(3) + sigL0) * DepsIL
      stv(17)  = DDE
      stv(18)  = stv(18) + stv(17) ! dissipated energy per unit volume
      if (dt > 0.0d0) stv(17)=stv(17)/dt  ! dissipated power per unit volume
      !
      if (stv(4) >= 0.0d0) then
        stv(13) = stv(13) + aLength * DepsIN ! wN
        stv(14) = stv(14) + aLength * DepsIM ! wM
        stv(11) = stv(11) + aLength * DepsIL ! wL
        stv(15) = DSQRT(DABS(stv(13)*stv(13) + stv(14)*stv(14) + 
     1                                         stv(11)*stv(11))) ! wtotal
        stv(19) = DDE
        stv(20) = stv(20) + stv(19) ! dissipated energy during fracture per unit volume
        if (dt > 0.0d0) stv(19) = stv(19)/dt  ! dissipated power during fracture per unit volume
      endif
      return
      end



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE UMAT                                                  C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================! 


!      subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt,
!     1 drplde,drpldt,stran,dstran,time,dtime,temp2,dtemp,predef,dpred,
!     2 cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
!     3 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

!      implicit none
       

      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
!      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare integers
!      integer(KIND=4), dimension(4)           :: jstep
!      integer(KIND=4)                         :: kinc
!      integer(KIND=4)                         :: kspt
!      integer(KIND=4)                         :: ndi  
!      integer(KIND=4)                         :: nelem
!      integer(KIND=4)                         :: nprops
!      integer(KIND=4)                         :: nshr
!      integer(KIND=4)                         :: ntens
!      integer(KIND=4)                         :: nstatv
!      integer(KIND=4)                         :: noel
!      integer(KIND=4)                         :: noffset
!      integer(KIND=4)                         :: npt
!      integer(KIND=4)                         :: layer       


      ! Declare double precision floats
!      real(dp)                                :: celent
!      real(dp), dimension(3)                  :: coords
!      real(dp), dimension(ntens,ntens)        :: ddsdde
!      real(dp), dimension(ntens)              :: ddsddt
!      real(dp), dimension(3,3)                :: dfgrd0
!      real(dp), dimension(3,3)                :: dfgrd1
!      real(dp), dimension(1)                  :: dpred      
!      real(dp), dimension(3,3)                :: drot
!      real(dp), dimension(ntens)              :: drplde
!      real(dp)                                :: drpldt
!      real(dp), dimension(ntens)              :: dstran
!      real(dp)                                :: dtemp
!      real(dp)                                :: dtime
!      real(dp)                                :: pnewdt
!      real(dp), dimension(1)                  :: predef
!      real(dp), dimension(nprops)             :: props
!      real(dp)                                :: rpl
!      real(dp)                                :: scd
!      real(dp)                                :: sse
!      real(dp)                                :: spd
!      real(dp), dimension(ntens)              :: stress
!      real(dp), dimension(nstatv)             :: statev  
!      real(dp), dimension(ntens)              :: stran
!      real(dp)                                :: temp2
!      real(dp), dimension(2)                  :: time


      ! Declare characters
!      character(8)                            :: cmname

!      ddsdde=0.0d0      
!      noffset=noel-nelem
     
!      return
!      end


C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE VUMAT                                                 C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================! 

!      subroutine vumat(
C     Read only (unmodifiable)variables -
!     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
!     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
!     3  props, density, strainInc, relSpinInc,
!     4  tempOld, stretchOld, defgradOld, fieldOld,
!     5  stressOld, stateOld, enerInternOld, enerInelasOld,
!     6  tempNew, stretchNew, defgradNew, fieldNew,
C     Write only (modifiable) variables -
!     7  stressNew, stateNew, enerInternNew, enerInelasNew 
C     Read only extra arguments -
!     8    nElement, nMatPoint, nLayer, nSecPoint )
C
!      implicit double precision (A-H,O-Z)   
      
C
!      dimension props(nprops), density(nblock), coordMp(nblock,*),
!     1  charLength(nblock), strainInc(nblock,ndir+nshr),
!     2  relSpinInc(nblock,nshr), tempOld(nblock),
!     3  stretchOld(nblock,ndir+nshr),
!     4  defgradOld(nblock,ndir+nshr+nshr),
!     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
!     6  stateOld(nblock,nstatev), enerInternOld(nblock),
!     7  enerInelasOld(nblock), tempNew(nblock),
!     8  stretchNew(nblock,ndir+nshr),
!     8  defgradNew(nblock,ndir+nshr+nshr),
!     9  fieldNew(nblock,nfieldv),
!     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
!     2  enerInternNew(nblock), enerInelasNew(nblock)
!      dimension nElement(nblock)
C
!      character*80 cmname
C

!      do kblock = 1,nblock
!          noel=nElement(kblock)          
!          noffset=noel-nelem
          !stressNew(kblock,:)=UserVar(noffset,1:ntens,npt)
          !stran=UserVar(noffset,ntens+1:2*ntens,npt)
C         stateNew(kblock,1:nstatev)=UserVar(noffset,1:nstatev,npt)
!      end do

!      return
!      end
C*****************************************************************            
C
C*****************************************************************
!      subroutine kstatevar(igauss,nsvint,statev,statev_ip,icopy)
C
C Transfer data to/from element-level state variable array from/to
C material-point level state variable array.
C
!      implicit double precision (A-H,O-Z)   
      

!      dimension statev(*),statev_ip(*)

!      isvinc=(igauss-1)*nsvint ! integration point increment

!      if (icopy.eq.1.0d0) then ! Prepare arrays for entry into umat
!       do i=1,nsvint
!        statev_ip(i)=statev(i+isvinc)
!       enddo
!      else ! Update element state variables upon return from umat
!       do i=1,nsvint
!        statev(i+isvinc)=statev_ip(i)
!       enddo
!      end if

!      return
!      end



C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     SUBROUTINE VUEXTERNALDB                                          C
C                                                                      C
C**********************************************************************C
C**********************************************************************C


      !================================================================!
      ! Explanation of variables                                       !
      !================================================================!
      !                                                                !
      !                                                                !
      !================================================================! 

      subroutine VUEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)


      use ModuleLDPM
      implicit none  


      ! Declare integers
      integer(KIND=4)                         :: KINC
      integer(KIND=4)                         :: KSTEP
      integer(KIND=4)                         :: LRESTART
      integer(KIND=4)                         :: LOP


      ! Declare double precision floats
      real(dp), dimension(2)                  :: time
      real(dp)                                :: DTIME


      ! Declare characters    
      character(256)                          :: OUTDIR
      character(256)                          :: Chemin
      character(256)                          :: JOBNAME        


      call get_elem_con
      call get_facet_info
      call get_svars_config


      RETURN
      END
C*****************************************************************            
C
C*****************************************************************
!      subroutine KMASSCONTINUUM(nnode,ndofn,ndofel,ninpt,mcrd,density,
!     1                          coords,nconnec,faceprop,amatrx) 
      ! 
!      implicit double precision (A-H,O-Z)   
c
c     dimension amatrx(ndofel, ndofel),coords(mcrd, nnode),
c     1          faceprop(20,12),fcoord(3),nconnec(4),
c     2          pMi(ndofn,ndofn),pMj(ndofn,ndofn)      
             
!      dimension amatrx(ndofel, ndofel),posgp(ninpt,mcrd),
!     1          weigp(ninpt),shapef(nnode),gpcod(mcrd),
!     2          dNdxi(nnode,mcrd),dNdx(nnode,mcrd),
!     3          coords(mcrd, nnode), alen(3), pmi(ndofn,ndofn)         
c
c     INITIALIZE THE MASS MATRÝX
c
!      do idofel = 1, NDOFEL       
!         do jdofel = 1, NDOFEL
!            amatrx(jdofel,idofel) = 0.0d0
!         enddo
!      enddo      
C    
C     Defining Gauss points coordinates and Gauss Gauss weights         
!      call gaussQtet(ninpt,ngauss,mcrd,posgp,weigp)
C
C     LOOP OVER INTEGRATION POINTS
C      
!      do igauss=1,ninpt  
C
C compute shape functions
C       
!      call SFR2(shapef,dNdxi,posgp(igauss,:),mcrd,nnode)
C                 
c compute jacobian and derivatives of shape functions at the end of increment
!      call CalcJacob(dNdx,dNdxi,det,coords,gpcod,
!     1               igauss,mcrd,nnode,shapef,ninpt,alen)
C          
!      dvol=det*weigp(igauss)           
C 
C Compute Mass matrix
C
!      dmass=density*dvol
      !
!      dx=posgp(igauss,1)*alen(1)
!      dy=posgp(igauss,2)*alen(2)
!      dz=posgp(igauss,3)*alen(3)
      !
!      pmi=0.0d0
!      pmi(1,1)=dmass
!      pmi(2,2)=dmass
!      pmi(3,3)=dmass
!      if (ndofn.gt.3) then
      !
!      rx2=dy*dy+dz*dz
!      ry2=dx*dx+dz*dz
!      rz2=dx*dx+dy*dy
!      rxy=dx*dy
!      rxz=dx*dz
!      ryz=dy*dz
!      Sx=dx*dmass*0; Sy=dy*dmass*0; Sz=dz*dmass*0           
      !
!      pmi(1,5)=Sz
!      pmi(1,6)=-Sy
!      pmi(2,4)=-Sz
!      pmi(2,6)=Sx
!      pmi(3,4)=Sy
!      pmi(3,5)=-Sx
      
!      pmi(5,1)=Sz
!      pmi(6,1)=-Sy
!      pmi(4,2)=-Sz
!      pmi(6,2)=Sx
!      pmi(4,3)=Sy
!      pmi(5,3)=-Sx
      !
!      pmi(4,4)=dmass*rx2
!      pmi(5,5)=dmass*ry2
!      pmi(6,6)=dmass*rz2
!      pmi(4,5)=-dmass*rxy
!      pmi(4,6)=-dmass*rxz
!      pmi(5,4)=-dmass*rxy
!      pmi(5,6)=-dmass*ryz
!      pmi(6,4)=-dmass*rxz
!      pmi(6,5)=-dmass*ryz
!      end if           
      !      
c    --------------------------------------------------------------------------------
!      do inode=1,nnode
!        irow=(inode-1)*ndofn          
!        do jnode=1,nnode
!           icol=(jnode-1)*ndofn              
c
c Assemble mass matrix
c          
!           do idofn=1,ndofn
!              ii=irow+idofn                         
!              do jdofn=1,ndofn
!                 jj=icol+jdofn   
!                 amatrx(ii,jj)=amatrx(ii,jj)+pmi(idofn,jdofn)*shapef(inode)*shapef(jnode)
!              end do
!           end do     
!        end do ! end of jnode   
!      end do ! end of inode
c    --------------------------------------------------------------------------------
      
!      end do !END of LOOP OVER INTEGRATION POINTS
       
       !write(*,*) "massmatrix"
       !do ii=1,ndofel
       !    write(*,"(24f20.8)") (amatrx(ii,jj),jj=1,ndofel)
       !end do
       
c
c Evaluate a lumped mass matrix using the row sum method
c      
!        do i = 1 , ndofel
!          do j = 1 , ndofel
!              if (i.ne.j) then
!                  amatrx(i,i) = amatrx(i,i) + amatrx(i,j)*0
!                  amatrx(i,j) = 0.0d0
!              end if
!          end do
!      end do
c       
      !write(*,*) "element"
      !do i=1,ndofel
      !    write(*,*) amatrx(i,i)
      !enddo 
!      return
!      end 
C*****************************************************************            
C
C*****************************************************************
!      subroutine gaussQtet(ninpt,ngauss,ndime,posgp,weigp) 
!
!
!*** THIS SUBROUTINE SETS UP THE GAUSS-LEGENDRE INTEGRATION CONSTANTS
!
!
!      implicit double precision (A-H,O-Z)   
!      dimension posgp(ninpt,ndime),weigp(ninpt)               
!      if ( ninpt == 1 ) then      

!      posgp(1,1) =   0.25D+00  
!      posgp(1,2) =   0.25D+00
!      posgp(1,3) =   0.25D+00

!      weigp(1) = 0.166666666666667D+00      

!      else if ( ninpt == 4 ) then

!      posgp(1,1:3) =   (/0.58541020D0,  0.13819660D0,  0.13819660D0/)   
!      posgp(2,1:3) =   (/0.13819660D0,  0.58541020D0,  0.13819660D0/)
!      posgp(3,1:3) =   (/0.13819660D0,  0.13819660D0,  0.58541020D0/)
!      posgp(4,1:3) =   (/0.13819660D0,  0.13819660D0,  0.13819660D0/)

!      weigp(1) = 0.041666666666667D+00;weigp(2) = 0.041666666666667D+00
!      weigp(3) = 0.041666666666667D+00;weigp(4) = 0.041666666666667D+00

!      else if ( ninpt == 5 ) then

!      posgp(1,1:3) =  (/0.250D0,  0.250D0,  0.250D0/)   
!      posgp(2,1:3) =  (/0.50D0,0.166666666666667D0,0.166666666666667D0/)
!      posgp(3,1:3) =  (/0.166666666666667D0,0.50D0,0.166666666666667D0/)
!      posgp(4,1:3) =  (/0.166666666666667D0,0.166666666666667D0,0.50D0/)
!      posgp(5,1:3) =  (/0.166666666666667D0,0.166666666666667D0,
!     1                0.166666666666667D0/)

!      weigp(1) = -0.133333333333333D+00;       weigp(2) = 0.075D+00
!      weigp(3) = 0.075D+00;                    weigp(4) = 0.075D+00
!      weigp(5) = 0.075D+00
     
!      end if      
c                
!      return
!      end subroutine  
C*****************************************************************            
C
C*****************************************************************
!!      subroutine SFR2(shapef,dNdxi,coord,ndime,nnode)
!!
!!    THIS SUBROUTINE EVALUATES SHAPE FUNCTIONS AND THEIR DERIVATIVES
!!    FOR LINEAR,QUADRATIC LAGRANGIAN AND SERENDIPITY
!!    ISOPARAMETRIC 2-D ELEMENTS
!!
!      implicit double precision (A-H,O-Z)   
!      dimension dNdxi(nnode,ndime),shapef(nnode),coord(ndime)
!      data R1,R2,ZERO,P25/1.0D0, 2.0D0, 0.D0, 0.25D0/             
!      xi=coord(1)
!      eta=coord(2)
!      zeta=coord(3)
!      if(nnode.EQ.4) then      
!!
!!*** SHAPE FUNCTIONS FOR 4 NODED TETRAHEDRAL ELEMENT
!!
!         shapef(1)=R1-xi-eta-zeta
!         shapef(2)=xi
!         shapef(3)=eta
!         shapef(4)=zeta
!!
!!*** SHAPE FUNCTION DERIVATIVES
!!
!         dNdxi(1,1)= -R1
!         dNdxi(2,1)= R1
!         dNdxi(3,1)= ZERO
!         dNdxi(4,1)= ZERO
!            
!         dNdxi(1,2)= -R1
!         dNdxi(2,2)= ZERO
!         dNdxi(3,2)= R1
!         dNdxi(4,2)= ZERO
!            
!         dNdxi(1,3)= -R1
!         dNdxi(2,3)= ZERO
!         dNdxi(3,3)= ZERO
!         dNdxi(4,3)= R1
!      
!      elseif(nnode.eq.10) then
!!
!!*** SHAPE FUNCTIONS FOR 10 NODED TETRAHEDRAL ELEMENT
!!
!         xi=coord(1)
!         eta=coord(2)
!         zeta=coord(3)
!      
!         shapef(1) = (R1-xi-eta-zeta)*(R1-R2*xi-R2*eta-R2*zeta);
!         shapef(2) = xi*(R2*xi-R1);
!         shapef(3) = eta*(R2*eta-R1);
!         shapef(4) = zeta*(R2*zeta-R1);
!         shapef(5) = 4.0d0*xi*(R1-xi-eta-zeta);
!         shapef(6) = 4.0d0*xi*eta;
!         shapef(7) = 4.0d0*eta*(R1-xi-eta-zeta);
!         shapef(8) = 4.0d0*zeta*(R1-xi-eta-zeta);
!         shapef(9) = 4.0d0*xi*zeta;
!         shapef(10) = 4.0d0*eta*zeta; 
!!
!!*** SHAPE FUNCTION DERIVATIVES
!!
!         dNdxi(1,1)= -3.0d0+4.0d0*xi+4.0d0*eta+4.0d0*zeta
!         dNdxi(2,1)= 4.0d0*xi - R1
!         dNdxi(3,1)= ZERO
!         dNdxi(4,1)= ZERO
!         dNdxi(5,1)= 4.0d0 - 8.0d0*xi - 4.0d0*eta - 4.0d0*zeta
!         dNdxi(6,1)= 4.0d0*eta
!         dNdxi(7,1)= -4.0d0*eta
!         dNdxi(8,1)= -4.0d0*zeta
!         dNdxi(9,1)= 4.0d0*zeta
!         dNdxi(10,1)= ZERO
!            
!         dNdxi(1,2)= -3.0d0 + 4.0d0*xi + 4.0d0*eta + 4.0d0*zeta
!         dNdxi(2,2)= ZERO
!         dNdxi(3,2)= 4.0d0*eta-R1
!         dNdxi(4,2)= ZERO
!         dNdxi(5,2)= -4.0d0*xi
!         dNdxi(6,2)= 4.0d0*xi
!         dNdxi(7,2)= 4.0d0 - 4.0d0*xi - 8.0d0*eta - 4.0d0*zeta
!         dNdxi(8,2)= -4.0d0*zeta
!         dNdxi(9,2)= ZERO
!         dNdxi(10,2)= 4.0d0*zeta
!            
!         dNdxi(1,3)= -3.0d0 + 4.0d0*xi + 4.0d0*eta + 4.0d0*zeta
!         dNdxi(2,3)= ZERO
!         dNdxi(3,3)= ZERO
!         dNdxi(4,3)= 4.0d0*zeta - R1
!         dNdxi(5,3)= -4.0d0*xi
!         dNdxi(6,3)= ZERO
!         dNdxi(7,3)= -4.0d0*eta
!         dNdxi(8,3)= 4.0d0 - 4.0d0*xi - 4.0d0*eta - 8.0d0*zeta
!         dNdxi(9,3)= 4.0d0*xi
!         dNdxi(10,3)= 4.0d0*eta
!          
!      else
!        write(*,*)"ERROR: element is not supported by user subroutine"
!        stop      
!      endif
!      return
!      end subroutine SFR2      
C*****************************************************************            
C
C*****************************************************************      
!     subroutine CalcJacob(dNdx,dNdxi,det,elcod,gpcod,
!     1            igauss,ndime,nnode,shapef,ninpt,alen) 
!
!*** THIS SUBROUTINE EVALUATES THE JACOBIAN MATRIX AND THE CARTESIAN
!    SHAPE FUNCTION DERIVATIVES       
!      
!      implicit double precision (A-H,O-Z)       
!     dimension dNdxi(nnode,ndime),elcod(ndime,nnode),shapef(nnode),
!    1            dNdx(nnode,ndime),gpcod(ndime),alen(3),
!     2                  dxidx(ndime,ndime),dxdxi(ndime,ndime)
!      logical Ok_flag
 !       DATA R0/0.0D0/
!      parameter (eps=1.0d-10)       
                 
!
!*** CALCULATE COORDINATES OF SAMPLING POINT
!
!      do  idime=1,ndime
!          gpcod(idime)=R0
!          do  inode=1,nnode
!              gpcod(idime)=gpcod(idime)+
!     1                   elcod(idime,inode)*shapef(inode)
!          end do
!      end do 
         
!
!*** CREATE JACOBIAN MATRIX dxdxi
!
!      do  idime=1,ndime
!          do  jdime=1,ndime
!              dxdxi(idime,jdime)=R0
!              do  inode=1,nnode
!                  dxdxi(idime,jdime)=dxdxi(idime,jdime)+
!     1                   dNdxi(inode,jdime)*elcod(idime,inode)
!              end do
!          end do
!      end do
      !write(*,*) "dxdxi",dxdxi      
!
!*** CALCULATE DETERMINANT AND INVERSE OF JACOBIAN MATRIX
!     
!!!!!!!!!! DIKKAT :: 3x3 matris tersi alan bir subroutine lazým buraya
!      if (ndime<3) then
!          det=dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)     
!          dxidx(1,1)=dxdxi(2,2)/det
!          dxidx(2,2)=dxdxi(1,1)/det
!          dxidx(1,2)=-dxdxi(1,2)/det
!          dxidx(2,1)=-dxdxi(2,1)/det
!      else
!          det=dxdxi(1,1)*(dxdxi(2,2)*dxdxi(3,3)-dxdxi(3,2)*dxdxi(2,3)) 
!     1        - dxdxi(1,2)*(dxdxi(2,1)*dxdxi(3,3)-dxdxi(2,3)*dxdxi(3,1))
!     2        + dxdxi(1,3)*(dxdxi(2,1)*dxdxi(3,2)-dxdxi(3,1)*dxdxi(2,2))
!          alen(1)=norm2(dxdxi(1,:))
!          alen(2)=norm2(dxdxi(2,:))
!          alen(3)=norm2(dxdxi(3,:))
!        if (dabs(DET) .LE. EPS) then          
!         write(*,*) "jacobian problem", eps,det         
!         dxidx = 0.0D0
!         write(*,1003) lmn,igauss
!1003     Format(//' *** Error in subroutinejacobian.f90 ***' 
!     1       //   ' Unable to invert jacobian matrix for element ',i4/ 
!     2         ' at integration point ',i4)
!                 stop          
!      end if
!      dxidx(1,1) = +(dxdxi(2,2)*dxdxi(3,3)-dxdxi(2,3)*dxdxi(3,2))/DET
!      dxidx(2,1) = -(dxdxi(2,1)*dxdxi(3,3)-dxdxi(2,3)*dxdxi(3,1))/DET
!      dxidx(3,1) = +(dxdxi(2,1)*dxdxi(3,2)-dxdxi(2,2)*dxdxi(3,1))/DET
!      dxidx(1,2) = -(dxdxi(1,2)*dxdxi(3,3)-dxdxi(1,3)*dxdxi(3,2))/DET
!      dxidx(2,2) = +(dxdxi(1,1)*dxdxi(3,3)-dxdxi(1,3)*dxdxi(3,1))/DET
!      dxidx(3,2) = -(dxdxi(1,1)*dxdxi(3,2)-dxdxi(1,2)*dxdxi(3,1))/DET
!      dxidx(1,3) = +(dxdxi(1,2)*dxdxi(2,3)-dxdxi(1,3)*dxdxi(2,2))/DET
!      dxidx(2,3) = -(dxdxi(1,1)*dxdxi(2,3)-dxdxi(1,3)*dxdxi(2,1))/DET
!      dxidx(3,3) = +(dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1))/DET   
!      endif      
!     
!*** CALCULATE CARTESIAN DERIVATIVES
!
!      dNdx=R0       
!      do idime=1,ndime
!          do jdime=1,ndime
!              do inode=1,nnode       
!                  dNdx(inode,idime)=dNdx(inode,idime)+
!     1            dxidx(jdime,idime)*dNdxi(inode,jdime)
!              end do
!          end do
!      end do
!   !   write(*,*) "dNdx",dNdx
! 600  FORMAT(//,36H PROGRAM HALTED IN SUBROUTINE JACOB2,/,11X
!     1       ,22H ZERO OR NEGATIVE AREA,/,10X,16H ELEMENT NUMBER ,I5)
!      return
!      end subroutine 

C      
C*****************************************************************
C*****************************************************************
C
      subroutine cal_stable_time(nnode,ndofn,ndofel,ninpt,mcrd,density,
     1                  elcoord,elecon,facetvar,elstif,omax)


      implicit none  


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


      ! Declare parameters
      real(dp), parameter :: eps = 1.0d-10


      ! Declare integers
      integer(KIND=4)                         :: i
      integer(KIND=4)                         :: idofn
      integer(KIND=4)                         :: ielem
      integer(KIND=4)                         :: it_max
      integer(KIND=4)                         :: it_num
      integer(KIND=4)                         :: j
      integer(KIND=4)                         :: k
      integer(KIND=4)                         :: nnode  
      integer(KIND=4)                         :: ndofn  
      integer(KIND=4)                         :: ndofel
      integer(KIND=4)                         :: ninpt
      integer(KIND=4)                         :: mcrd  
      integer(KIND=4)                         :: one = 1


      ! Declare double precision floats
      real(dp)                                :: a
      real(dp)                                :: density
      real(dp), dimension(mcrd,nnode)         :: elcoord
      real(dp), dimension(4)                  :: elecon
      real(dp), dimension(30,12)              :: FacetVar
      real(dp), dimension(ndofel,ndofel)      :: elmass
      real(dp), dimension(ndofel,ndofel)      :: Elstif
      real(dp), dimension(ndofel,ndofel)      :: elmassInvStiff
      real(dp)                                :: irot_num
      real(dp), dimension(ndofel,ndofel)      :: veigen
      real(dp), dimension(ndofel)             :: deigen
      real(dp), dimension(ndofel,ndofel)      :: CML
      real(dp), dimension(ndofel,ndofel)      :: CMLinv
      real(dp), dimension(ndofel,ndofel)      :: CMLinvT
      real(dp)                                :: omax


      elmass = 0.0d0      

      call KMASSLDPM(nnode,ndofn,ndofel,ninpt,mcrd,density,
     1               elcoord,elecon,FacetVar,elmass)  
      
        
      !----------------------------------------------------------------------------------------!
      !----------------------------------------------------------------------------------------!
      ! Choleski factorization 

      CML = 0.0d0
      
      CML(1,1) = DSQRT(DABS(elmass(1,1)))

      do i = 2,ndofel
            if (elmass(1,1).lt.EPSILON(elmass(1,1))) then
                  write(*,*)'Warning: Element mass is numerically equivalent to zero; failure is likely. Elmass = ',Elmass(1,1)
            endif
            CML(i,1) = elmass(1,i)/elmass(1,1)
      enddo

      do j = 2 , ndofel
            do i = 1 , ndofel
                  if(i.lt.j) then
                        CML(i,j) = 0.0d0
                  else
                        if(i.eq.j) then
                              a = 0.0d0
                              do k = 1 , i-one
                                   a = a + CML(i,k)**2.0d0
                              enddo
                              CML(i,j) = DSQRT(DABS(elmass(i,i) - a))
                        else
                              a = 0.0d0
                              do k = 1 , i-one
                                   a = a + CML(j,k) * CML(i,k)
                              enddo
                              CML(i,j) = (elmass(j,i) - a)/CML(j,j)
                        endif
                  endif
            enddo
      enddo      
      
      !write(*,*) "CML =" 
      !do i=1,ndofel
      !    write(*,'(24E18.10)') (CML(i,j),j=1,ndofel)
      !end do  
      
      CMLinv  = 0.0d0
      do idofn = 1,ndofel 
          CMLinv(idofn,idofn) = 1.0d0/CML(idofn,idofn)          
      end do      
      CMLinvT = 0.0d0
      CMLinvT = transpose(CMLinv)
      elmassInvStiff = matmul(CMLinv,matmul(elstif,CMLinvT))
      
      !write(*,*) "massInvstiff =" 
      !do i=1,ndofel
      !    write(*,'(24E18.10)') (elmassInvStiff(i,j),j=1,ndofel)
      !end do  
      
      !----------------------------------------------------------------------------------------! 
      !----------------------------------------------------------------------------------------! 
      it_max = 1000000
      call jacobi_eigenvalue (ndofel,elmassInvStiff,it_max,veigen,deigen,it_num,irot_num,ielem)
      omax = 2.0d0/DSQRT(DABS(deigen(ndofel)))
      !write(*,*) "time step =", omax      
      return
      end subroutine
       
C
C***************************************************************** 
C***************************************************************** 

      subroutine jacobi_eigenvalue(n,a,it_max,v,d,it_num,rot_num,ielem)

C*****************************************************************
C*****************************************************************
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
!
!    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
!
!    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
!
        implicit none   

      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)


        integer ( kind = 4 ) n
        real(dp)                              :: a(n,n)
        real(dp)                              :: bw(n)
        real(dp)                              :: c
        real(dp) :: d(n)
        real(dp) :: g
        real(dp) :: gapq
        real(dp) :: h
        integer ( kind = 4 ) i
        integer ( kind = 4 ) it_max
        integer ( kind = 4 ) it_num
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ) l
        integer ( kind = 4 ) m
        integer ( kind = 4 ) p
        integer ( kind = 4 ) q
        integer ( kind = 4 ) rot_num
        real(dp) :: s
        real(dp) :: t
        real(dp) :: tau
        real(dp) :: term
        real(dp) :: termp
        real(dp) :: termq
        real(dp) :: theta
        real(dp) :: thresh
        real(dp) :: v(n,n)
        real(dp) :: w(n)
        real(dp) :: zw(n)
        integer (kind = 4 ) ielem
        
        do j = 1, n
          do i = 1, n
            v(i,j) = 0.0D+00
          end do
          v(j,j) = 1.0D+00
        end do

        do i = 1, n
          d(i) = a(i,i)
        end do

        bw(1:n) = d(1:n)
        zw(1:n) = 0.0D+00
        it_num = 0
        rot_num = 0

        do while ( it_num < it_max )

          it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
          thresh = 0.0D+00
          do j = 1, n
            do i = 1, j - 1
              thresh = thresh + a(i,j) ** 2.0d0
            end do
          end do

          thresh = DSQRT(DABS(thresh)) / real ( 4 * n, dp )

          if ( thresh == 0.0D+00 ) then
!           write(*,*) "Threshhold is zero in eigen value calculation!"
            exit 
          end if

          do p = 1, n
            do q = p + 1, n

              gapq = 10.0D+00 * DABS( a(p,q) )
              termp = gapq + DABS( d(p) )
              termq = gapq + DABS( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
             if ( 4 < it_num.and.termp == DABS(d(p)).and.termq == DABS(d(q))) then
                a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
!
                   else if ( thresh <= DABS( a(p,q) ) ) then
                      h = d(q) - d(p)
                      term = DABS( h ) + gapq

                   if ( term == DABS( h ) ) then
                      t = a(p,q) / h
                   else
                      theta = 0.5D+00 * h / a(p,q)
                      t = 1.0D+00 / ( DABS( theta ) + DSQRT(DABS( 1.0D+00 + theta * theta )))
                   if ( theta < 0.0D+00 ) then 
                          t = -t
                   end if
                   end if

                c = 1.0D+00 / DSQRT(DABS ( 1.0D+00 + t * t ))
                s = t * c
                tau = s / ( 1.0D+00 + c )
                h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
!
                zw(p) = zw(p) - h                  
                zw(q) = zw(q) + h
                d(p) = d(p) - h
                d(q) = d(q) + h
                a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
                do j = 1, p - 1
                  g = a(j,p)
                  h = a(j,q)
                  a(j,p) = g - s * ( h + g * tau )
                  a(j,q) = h + s * ( g - h * tau )
                end do
C
                do j = p + 1, q - 1
                  g = a(p,j)
                  h = a(j,q)
                  a(p,j) = g - s * ( h + g * tau )
                  a(j,q) = h + s * ( g - h * tau )
                end do
C
                do j = q + 1, n
                  g = a(p,j)
                  h = a(q,j)
                  a(p,j) = g - s * ( h + g * tau )
                  a(q,j) = h + s * ( g - h * tau )
                end do
!
!  Accumulate information in the eigenvector matrix.
!
                do j = 1, n
                  g = v(j,p)
                  h = v(j,q)
                  v(j,p) = g - s * ( h + g * tau )
                  v(j,q) = h + s * ( g - h * tau )
                end do

                rot_num = rot_num + 1

           end if

      end do
      end do

      bw(1:n) = bw(1:n) + zw(1:n)
      d(1:n) = bw(1:n)
      zw(1:n) = 0.0D+00

      end do
!
!  Restore upper triangle of input matrix.
!
      do j = 1, n
         do i = 1, j - 1
            a(i,j) = a(j,i)
         end do
      end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
      do k = 1, n - 1
         m = k
         do l = k + 1, n
            if ( d(l) < d(m) ) then
              m = l
            end if
         end do
         if ( m /= k ) then
            t    = d(m)
            d(m) = d(k)
            d(k) = t
            w(1:n)   = v(1:n,m)
            v(1:n,m) = v(1:n,k)
            v(1:n,k) = w(1:n)
         end if
      end do
      return
      end
C***************************************************************** 
C*****************************************************************
C User subroutine VFRIC
      subroutine vfric (fTangential,statev, kStep, kInc,
     1             nContact, nFacNod, nSlvNod, nMstNod, nFricDir, 
     2            nDir, nStateVar, nProps, nTemp, nPred, 
     3            numDefTfv, jSlvUid, jMstUid, jConSlvid, jConMstid, 
     4            timStep, timGlb,dTimPrev, surfInt, surfSlv, 
     5            surfMast, lContType, dSlipFric, fStickForce, 
     6            fTangPrev, fNormal, frictionWork,shape,
     7            coordSlv, coordMst, dircosSl, dircosN, props, 
     8            areaSlv, tempSlv, preDefSlv, tempMst, preDefMst )

      implicit none   


      ! Initialize double precision (Lemmon and Schafer, 2005, p. 20–22)
      integer,  parameter :: dp = selected_real_kind(15, 307)      


      ! Declare parameters
      real(dp), parameter :: eps = 1.0d-10


      ! Declare integers
      integer(KIND=4), dimension(nSlvNod)     :: jSlvUid
      integer(KIND=4), dimension(nMstNod)     :: jMstUid
      integer(KIND=4), dimension(nContact)    :: jConSlvid
      integer(KIND=4), dimension(nFacNod,nContact)  :: jConMstid
      integer(KIND=4)                         :: kcon
      integer(KIND=4)                         :: kStep
      integer(KIND=4)                         :: kInc
      integer(KIND=4)                         :: nContact
      integer(KIND=4)                         :: nFacNod
      integer(KIND=4)                         :: nSlvNod
      integer(KIND=4)                         :: nMstNod
      integer(KIND=4)                         :: nFricDir
      integer(KIND=4)                         :: ndir
      integer(KIND=4)                         :: nStateVar
      integer(KIND=4)                         :: nProps
      integer(KIND=4)                         :: nTemp
      integer(KIND=4)                         :: nPred
      integer(KIND=4)                         :: numDefTfv
      integer(KIND=4)                         :: lContType


      ! Declare double precision floats
      real(dp), dimension(nSlvNod)            :: areaSlv
      real(dp), dimension(nDir,nMstNod)       :: coordSlv
      real(dp), dimension(nDir,nMstNod)       :: coordMst
      real(dp), dimension(nDir,nContact)      :: dirCosSl
      real(dp), dimension(nDir,nContact)      :: dircosN
      real(dp), dimension(nDir,nContact)      :: dSlipFric
      real(dp), dimension(nContact)           :: fNormal
      real(dp), dimension(nFricDir,nContact)  :: fTangential
      real(dp), dimension(nDir,nContact)      :: fTangPrev
      real(dp), dimension(nContact)           :: fStickForce
      real(dp), dimension(nContact,nPred)     :: preDefSlv
      real(dp), dimension(numDefTfv,nPred)    :: preDefMst
      real(dp), dimension(nProps)             :: props
      real(dp), dimension(nStateVar,nSlvNod)  :: statev
      real(dp), dimension(nContact)           :: tempSlv
      real(dp), dimension(numDefTfv)          :: tempMst
      real(dp)                                :: timStep
      real(dp)                                :: fn
      real(dp)                                :: fs
      real(dp)                                :: ft
      real(dp)                                :: timGlb
      real(dp)                                :: dtimPrev
      real(dp)                                :: frictionWork
      real(dp), dimension(nFacNod,nContact)   :: shape

    
      ! Declare characters
      character*80 surfInt, surfSlv, surfMast


      dimension slip(nContact)
      dimension dSlipNslv(nSlvNod)
      dimension jCon(nContact)
      real(dp) :: slip
      real(dp) :: dSlipNslv
      real(dp) :: jCon
      real(dp) :: mus
      real(dp) :: mud
      real(dp) :: so
      real(dp) :: mu

      mus = props(1)
      mud = props(2)
      so = props(3)
     
      
      ! Assign values
      dSlipNslv(:)=0.0d0
      slip(:)=0.0d0
      jCon(:)=0
C

      do kcon = 1, nContact
        jCon(kcon) = jSlvUid(JConSlvid(kcon))
        Slip(kcon) = statev(1,jConSlvid(kcon))+ dSlipFric(1,kcon)
        statev(1,jConSlvid(kcon)) = Slip(kcon)
      end do

C define fTangential for slip from previous increment at nodes that
C are in contact during this increment
      do kcon = 1, nContact
        fn = fNormal(kcon)
        fs = fStickForce(kcon)
        mu = mud+(mus-mud)*so/(so+Slip(kcon))
        ft = MIN(mu*fn,fs)
        fTangential(1,kcon)= -ft
      end do
      return
      end

