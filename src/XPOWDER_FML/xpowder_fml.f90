 !!----
 !!---- PROGRAM:  XPOWDER_FML (Calculation_of_Powder_Patterns)
 !!----
 !!---- Simple program for calculating powder patterns by reading a CIF 
 !!----
 !!---- If a CIF file is read there is no control on powder pattern parameters. A standard
 !!---- pattern corresponding to CuKa X-ray radiation is calculated. The parameters of this
 !!---- powder pattern are described by the components of the object PPC (powder pattern constants)
 !!---- Default Powder Diffraction Pattern
 !!----    stlmax=0.6; PPC%Title="Default Powder Pattern"
 !!----    PPC%U=0.0002; PPC%V=-0.0002; PPC%W=0.012; PPC%LAMBDA=1.54056; PPC%X=0.0015
 !!----    PPC%Thmin=1.00; PPC%step=0.05;  PPC%Thmax= int(2.0*asind(stlmax*1.54056)); PPC%job=0
 !!----    PPC%Ls=1900.0; nf=30; PPC%bkg=0.0
 !!----    powfile="powder_pattern.dat"
 !!----
 !!----
 !!---- The program uses CrysFML and a module called Gen_Powder_Pattern where the subroutine for
 !!---- calculating the powder diffraction pattern is stored
 !!----
 Module Gen_Powder_Pattern
    !---- Use Modules ----!
    use CFML_GlobalDeps,           only: to_Deg
    use CFML_Math_General,         only: asind,locate
    use CFML_Reflections_Utilities,only: Reflection_List_Type
    use CFML_Diffraction_Patterns, only: Diffraction_Pattern_Type, Allocate_Diffraction_Pattern
    use CFML_PowderProfiles_CW,    only: PseudoVoigt

    !---- Variables ----!
    implicit none

    private

    public  :: calc_powder_pattern
    public  :: tch
!srio    private :: TCH

    Type, public :: PowPat_CW_Conditions
       character(len=140) :: title
       integer :: job
       real    :: Lambda, U, V, W, X, Ls
       real    :: Thmin, Thmax, step
       real    :: scalef,bkg
    End Type PowPat_CW_Conditions

 Contains
    !!----
    !!---- Pure Subroutine TCH(Hg,Hl,Fwhm,Eta)
    !!----
    !!---- Calculation of eta and FWHM of the pV-function for the
    !!---- T-C-H representation.
    !!----
    !!
    Pure Subroutine TCH(Hg,Hl,Fwhm,Eta)
       !---- Arguments ----!
       real, intent(in)  :: hg
       real, intent(in)  :: hl
       real, intent(out) :: fwhm
       real, intent(out) :: eta

       !---- Variables ----!
       real, parameter :: o1= 2.69269, o2=2.42843, o3=4.47163, o4= 0.07842
       real, parameter :: e1= 1.36603, e2=0.47719, e3=0.11116
       real            :: ctl, tlr

       ! There is no exception handling because it is supposed to be
       ! perfomed before calling TCH
       ctl=hg**5.0+o1*hg**4.0*hl+o2*hg**3.0*hl**2.0+o3*hg**2.0*hl**3.0+  &
           o4*hg*hl**4.0+hl**5.0
       fwhm=ctl**0.2
       tlr = hl/fwhm
       eta = max(1.0e-06, e1*tlr-e2*tlr*tlr+e3*tlr**3.0)  !eta

       Return
    End Subroutine TCH

    Subroutine Calc_Powder_Pattern(Ppc,Hkl,Pat)
       !---- Argument ----!
       Type(PowPat_CW_Conditions),     intent(in)  :: PPC
       Type(Reflection_List_Type),     intent(in)  :: hkl
       Type(Diffraction_Pattern_Type), intent(out) :: Pat

       !--- Local Variables ----!
       integer :: i,j,npts,i1,i2
       real    :: step,Intens,Bragg,Hl,Hg, ss,cs,tt,th1,th2,LorentzF, Y,eta,fwhm,chw

       npts=(PPC%Thmax-PPC%Thmin)/PPC%step + 1.02
       call Allocate_Diffraction_Pattern(Pat,npts)
       Pat%Title=adjustl(Trim(PPC%title))
       i=len_trim(Pat%Title)
       write(unit=Pat%Title(i+2:),fmt="(a,f7.4,f7.1)") " => lambda,Ls: ", &
                  PPC%Lambda,PPC%Ls

       Pat%scat_var="2-Theta"
       Pat%instr="Calculated Pattern"
       Pat%xmin= PPC%Thmin
       Pat%xmax= PPC%Thmax
       Pat%ymin= 0.0
       Pat%ymax=0.0
       Pat%scal=1.0
       Pat%monitor=0.0
       Pat%step=PPC%step
       Pat%Tsamp=300.0
       Pat%Tset=300.0
       Pat%npts=npts
       Pat%ct_step=.true.
       Pat%conv=(/PPC%Lambda,PPC%Lambda,0.0/)
       chw=15.0
       do i=1,npts
          Pat%x(i)=Pat%xmin+real(i-1)*Pat%step
       end do

       Y= to_deg*PPC%Lambda/PPC%Ls
       do i=1,hkl%nref
          ss=PPC%Lambda*hkl%ref(i)%S
          cs=sqrt(abs(1.0-ss*ss))
          tt=ss/cs
          LorentzF=0.5/(ss*ss*cs)
          Bragg=2.0*asind(ss)
          HG=sqrt(tt*(PPC%U*tt+PPC%V)+PPC%W)
          HL=PPC%X*tt + Y/cs
          call TCH(hg,hl,fwhm,eta)
          Select Case(nint(eta*10.0))
             Case(:2)
                chw=25.0
             case(3:5)
                chw=45.0
             case(6:7)
                chw=60.0
             case(8:)
                chw=90.0
          End Select

          th1=Bragg-chw*fwhm
          th2=Bragg+chw*fwhm
          i1=Locate(Pat%x,npts,th1)
          i2=Locate(Pat%x,npts,th2)
          i1=max(i1,1)
          i2=min(i2,npts)
          Intens= LorentzF *hkl%ref(i)%mult * hkl%ref(i)%Fc**2 * PPC%Scalef
          do j=i1,i2
             Pat%ycalc(j)=Pat%ycalc(j)+ PseudoVoigt( Pat%x(j)-Bragg, (/fwhm,eta /) ) * Intens 
          end do
       end do

       return
    End Subroutine Calc_Powder_Pattern


  End Module Gen_Powder_Pattern

  !!----
  !!----  Program Calculation_of_Powder_Patterns
  !!----
  !!----
  !!---- Update: June - 2009
  !!
  Program XPOWDER_FML
     !---- Use Modules ----!
     use CFML_Math_General,              only: asind,sind
     use CFML_Atom_TypeDef,              only: Atom_List_Type, Allocate_Atom_List,Write_Atom_List
     use CFML_Crystal_Metrics,           only: Crystal_Cell_type, set_Crystal_Cell,write_crystal_cell
     use CFML_string_utilities,          only: u_case
     use CFML_Reflections_Utilities,     only: Reflection_List_Type,Hkl_uni,get_maxnumref
     Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Set_SpaceGroup, &
                                               Write_SpaceGroup
     Use CFML_Structure_Factors,         only: Write_Structure_Factors, Structure_Factors,&
                                               Init_Structure_Factors,ERR_SFac,ERR_SFac_Mess
     use CFML_Diffraction_Patterns,      only: Diffraction_Pattern_Type, &
                                               Allocate_Diffraction_Pattern
     use CFML_IO_Formats,                only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type

     use Gen_Powder_Pattern
     
     use CFML_GlobalDeps,                only: to_Deg

     !---- Variables ----!
     implicit none

     integer                :: i,j,k,l,m,n,maxnumref,ier,mult,nf
     integer                :: lun=1,lp=2
     real, dimension(3)     :: ad,ang,x,fr
     character(len=1)       :: ans,outa
     integer, dimension(3)  :: ncel,h
     real                   :: stlmax,tini,tfin,sn,sf2,tim,ftim,box
     character(len=132)     :: line,powfile,filcod
     character(len=3)       :: mode
     character(len=8)       :: units,radiation
     character(len=4),dimension(:),allocatable :: ch
     

     real    :: Intens,Bragg,Hl,Hg, ss,cs,tt,LorentzF, Y,eta,fwhm
     real    :: tmpf

     Type(Crystal_Cell_type)        :: cell
     Type(Space_Group_Type)         :: SpG
     Type(Atom_List_Type)           :: A
     Type(Reflection_List_Type)     :: hkl
     Type(Diffraction_Pattern_Type) :: Pat
     Type(PowPat_CW_Conditions)     :: PPC
     Type(file_list_type)           :: fich_cfl
     integer                        :: narg, tmpi
     Logical                        :: esta, arggiven=.false., fail


     !---- Arguments on the command line ----!
     narg=command_argument_count()

     if (narg > 0) then
        call get_command_argument(1,filcod)
        arggiven=.true.
     end if


     write(unit=*,fmt="(/,/,6(a,/))")                                                     &
          "              ------ PROGRAM POWDER PATTERN CALCULATION  ------"             , &
          "                     ---- Version 0.1 April-2009----"                        , &
          "    **********************************************************************"  , &
          "    * Calculates powder diffraction pattern from a *.CIF file            *"  , &
          "    **********************************************************************"  , &
          "                          (JRC- April 2009 )"
     write(unit=*,fmt=*) " "

     if (.not. arggiven) then
        write(unit=*,fmt="(a)", advance='no') " => cif file: "
        read(unit=*,fmt="(a)") filcod
        if(len_trim(filcod) == 0) stop
     end if


     inquire(file=trim(filcod),exist=esta)
     if (esta) then
        Mode="CIF"
        call Readn_set_Xtal_Structure(trim(filcod),Cell,SpG,A,Mode="CIF")
     end if

    if (err_form) then
       write(unit=*,fmt="(a)") trim(err_form_mess)
    else
       !open(unit=lun,file=trim(filcod)//".par", status="replace",action="write")
       open(unit=lun,file="xpowder_fml.par", status="replace",action="write")
       write(unit=lun,fmt="(/,/,6(a,/))")                                                 &
          "              ------ PROGRAM: POWDER PATTERN CALCULATION  ------"             , &
          "                     ---- Version 0.1 April-2009----"                        , &
          "    ***********************************************************************"  , &
          "    * Calculates powder diffraction patterns from a *.CIF file            *"  , &
          "    ***********************************************************************"  , &
          "                            (JRC- April 2009 )"
       ! Calculate a default Powder Diffraction Pattern
       stlmax=0.6; PPC%Title="Default Powder Pattern"
       PPC%U=0.0002; PPC%V=-0.0002; PPC%W=0.012; PPC%LAMBDA=1.54056; PPC%X=0.0015
       PPC%Thmin=1.00; PPC%step=0.05;  PPC%Thmax= int(2.0*asind(stlmax*1.54056)); PPC%job=0
       PPC%Ls=1900.0;  nf=30; PPC%bkg=0.0
       powfile="xpowder_fml.out"
       !powfile = trim(filcod)//".dat"
       units=" seconds"
       tim=0.0

       ! Write initial structure information in the .powder file
       call Write_Crystal_Cell(Cell,lun)
       call Write_SpaceGroup(SpG,lun)
       call Write_Atom_List(A,lun=lun)


       !
       ! define parameters for powder pattern
       !
       
       ! title
       write(*,fmt="(a)",advance="no") "New Title: ["//trim(PPC%TITLE)//"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) PPC%TITLE=line
       
       
       ! lambda
       write(*,fmt="(a,f12.6,a)",advance="no") "New Lambda: [",PPC%LAMBDA,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%LAMBDA=tmpf 
       end if

       !
       ! job
       !
       write(*,fmt="(a,I3,a)",advance="no") "New Job 0=XRays 1=Neutrons: [",PPC%Job,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpi
         if (ier == 0) PPC%job=tmpi
       end if   
       
       
       ! u v w x
       write(*,fmt="(a,f12.6,a)",advance="no") "New U: [",PPC%U,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%U=tmpf 
       end if
       write(*,fmt="(a,f12.6,a)",advance="no") "New V: [",PPC%V,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%V=tmpf 
       end if
       write(*,fmt="(a,f12.6,a)",advance="no") "New W [",PPC%W,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%W=tmpf 
       end if
       write(*,fmt="(a,f12.6,a)",advance="no") "New X: [",PPC%X,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%X=tmpf 
       end if
       
       
       ! bkg
       !write(*,fmt="(a,f12.6,a)",advance="no") "New BCK: [",PPC%BKG,"]: "
       !read(*,fmt="(a)") line
       !if (len_trim(line) > 0) then
       !  read(line,fmt=*,iostat=ier) tmpf
       !  if (ier == 0) PPC%BKG=tmpf 
       !end if
       
       !
       ! Ls
       !
       write(*,fmt="(a,f12.6,a)",advance="no") "New Ls: [",PPC%Ls,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%Ls=tmpf 
       end if   
       
       !
       ! thmin, step, thMax
       !
       write(*,fmt="(a,f12.6,a)",advance="no") "New ThMin: [",PPC%ThMin,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%ThMin=tmpf 
       end if
       write(*,fmt="(a,f12.6,a)",advance="no") "New Step: [",PPC%Step,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%Step=tmpf 
       end if
       write(*,fmt="(a,f12.6,a)",advance="no") "New ThMax: [",PPC%ThMax,"]: "
       read(*,fmt="(a)") line
       if (len_trim(line) > 0) then
         read(line,fmt=*,iostat=ier) tmpf
         if (ier == 0) PPC%ThMax=tmpf 
       end if         

                
         
         
         !write(*,fmt="(a,f12.6)") "Default Lambda: ",PPC%LAMBDA
         !write(*,fmt="(a)",advance="no") "New Lambda: "
         !if (len_trim(line) > 0) then 
         !  !print *,"Changing line..."
         !  read(line,fmt=*) PPC%LAMBDA
         !  !write(*,fmt="(a,f12.6)") "New Lambda: ",PPC%LAMBDA
         !end if


    ! Calculate sinTheta/Lambda max from 2Thetamax     PPC%Thmax= int(2.0*asind(stlmax*1.56))
    stlmax=sind(min((PPC%Thmax+10.0)*0.5,90.0))/PPC%lambda

    end if !if error

    if (PPC%job == 0) then      !X-rays
       write(unit=lun,fmt="(/,a)")  " => CALCULATION OF X-RAY POWDER DIFFRACTION PATTERN "
       PPC%title=Trim(PPC%title)//"; X-RAYS: "
    else
       write(unit=lun,fmt="(/,a)")  " => CALCULATION OF NEUTRON POWDER DIFFRACTION PATTERN"
       PPC%title=Trim(PPC%title)//"; NEUTRONS: "
    end if
    write(unit=lun,fmt="(  a,4f10.5)")  " => Resolution parameters UVWX: ",PPC%U,PPC%V,PPC%W,PPC%X
    write(unit=lun,fmt="(  a, f10.5)")  " => Lambda: ",PPC%lambda
    write(unit=lun,fmt="(  a, f10.5,a)")" => Background level: ",PPC%bkg," counts"
    write(unit=lun,fmt="(  a,2f10.2)")  " => Lorentzian size: ",PPC%Ls
    write(unit=lun,fmt="(  a,3f10.5)")  " => 2Theta range and step: ",PPC%Thmin,PPC%step,PPC%Thmax
    write(unit=lun,fmt="(  a,3f10.5)")  " => Maximum sin(Theta)/Lambda (for generating reflections): ",stlmax

    ! Now calculate a powder diffraction pattern
    ! First generate reflections and calculate structure factors
    Mult=2*SpG%NumOps
       MaxNumRef = get_maxnumref(stlmax,Cell%CellVol,mult=Mult)
       call cpu_time(tini)
       call Hkl_Uni(Cell,SpG,.true.,0.0,stlmax,"s",MaxNumRef,hkl)
       call cpu_time(tfin)
       tim=tim+ tfin-tini
       write(*,"(a,i8)") "  => Total number of generated reflections is ",hkl%nref
       write(unit=lun,fmt="(a,i9)") " => Total number of generated reflections is ",hkl%nref

       if (PPC%job == 1) then      !Neutrons
          call Init_Structure_Factors(hkl,A,Spg,mode="NUC",lun=lun)
       else if(PPC%job == 0) then !Xrays
          call Init_Structure_Factors(hkl,A,Spg,mode="XRA",lambda=PPC%lambda,lun=lun)
       end if

       call cpu_time(tini)
       write(*,*) " => Calculating structure factors ..."
       if (PPC%job == 1) then      !Neutrons
          call Structure_Factors(A,SpG,hkl,mode="NUC")
       else if(PPC%job == 0) then !X-rays
          call Structure_Factors(A,SpG,hkl,mode="XRA",lambda=PPC%lambda)
       end if
       call cpu_time(tfin)
       tim=tim+ tfin-tini
       write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for Structure_Factors: ",(tfin-tini)*ftim,units

       if (ERR_SFac) then
          write(*,*) " => Error in calculations of Structure Factors"
          write(*,*) " => "//trim(ERR_SFac_Mess)
          stop
      end if

      if (radiation(1:1) == "N") then
         call Write_Structure_Factors(lun,hkl,mode="NUC")
      else
         call Write_Structure_Factors(lun,hkl,mode="XRA")
      end if

      call cpu_time(tini)
      PPC%Scalef=cell%RCellVol
      call Calc_powder_pattern(PPC,hkl,Pat)
      call cpu_time(tfin)
      tim=tim+ tfin-tini
      write(*,*) " => CPU-time used for Calc_powder_pattern: ",(tfin-tini)*ftim,units
      write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time used for Calc_powder_pattern: ",(tfin-tini)*ftim,units
      write(*,*) " => CPU-time for all calculations: ",tim*ftim,units
      write(unit=lun,fmt="(a,f15.3,a)") " => CPU-time for all calculations: ",tim*ftim,units

      open(unit=lp,file=trim(powfile),status="replace",action="write")
      !  write(unit=lp,fmt="(a)") "!"//trim(Pat%title)
      !  write(unit=lp,fmt="(3f10.4)") Pat%xmin,Pat%step,Pat%xmax
      ! write(unit=lp,fmt="(8f16.4)") Pat%ycalc
      write(unit=lp,fmt="(a)") "#S 1 "//trim(Pat%title)
      write(unit=lp,fmt="(a)") "#N 2"
      write(unit=lp,fmt="(a)") "#L TwoTheta[Deg]  Intensity[a.u.]"
      do i=1,pat%npts
        write(unit=lp,fmt="(2f16.4)") Pat%xmin+(i-1)*Pat%step,Pat%ycalc(i)
      end do 
 
      
      
      close(unit=lp)
      
      
      
      open(unit=lp,file="xpowder_fml.ref",status="replace",action="write")
      write(unit=lp,fmt="(a)") "#S 1 "//trim(Pat%title)//" Reflections "
      write(unit=lp,fmt="(a)") "#N 10"
      write(unit=lp,fmt="(a)") "#L twoTheta  Intens  d-spacing  h  k  l  mult  F**2  fwhm  eta"
     
      Y= to_deg*PPC%Lambda/PPC%Ls
      
      do i=1,hkl%nref
          ss=PPC%Lambda*hkl%ref(i)%S
          cs=sqrt(abs(1.0-ss*ss))
          tt=ss/cs
          LorentzF=0.5/(ss*ss*cs)
          Bragg=2.0*asind(ss)
          HG=sqrt(tt*(PPC%U*tt+PPC%V)+PPC%W)
          HL=PPC%X*tt + Y/cs
          call TCH(hg,hl,fwhm,eta)
          Intens= LorentzF *hkl%ref(i)%mult * hkl%ref(i)%Fc**2 * PPC%Scalef
           write(unit=lp,fmt="(3f16.4,4I5,3f16.4)")  asin(ss)*2.0*to_deg, Intens, 0.5/hkl%ref(i)%s, &
                  hkl%ref(i)%h, hkl%ref(i)%mult, &
                  hkl%ref(i)%Fc**2, fwhm, eta
       end do
      
      
      close(unit=lp)
      
      close(unit=lun)
      
      
      stop

  End Program xpowder_fml

