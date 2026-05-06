MODULE Adjoint_Utils_Mod

        
! =====================================================================
! Utility to handle 
! 1. "Push/Pop" logic for Adjoint Operator Splitting
! 2.  Set the Initial conditions for adjoint state
! 3.  Integrate surface adjoint
! =====================================================================

#include "MAPL_Generic.h"
  USE MAPL_MOD
  USE ESMF
  USE Precision_Mod
  USE State_Chm_Mod,   ONLY: ChmState
  USE Input_Opt_Mod,   ONLY : OptInput
  USE State_Grid_Mod,  ONLY : GrdState
  USE State_Met_Mod,   ONLY : MetState
  USE netcdf

  IMPLICIT NONE
  PRIVATE

  ! Public Interface
  PUBLIC :: State_Snapshot
  PUBLIC :: Push_State
  PUBLIC :: Pop_State
  PUBLIC :: Setup_Adjoint_ForwardPert
  PUBLIC :: Integrate_Srf_Adjoint
  PUBLIC :: Load_OCO2_Adjoint_Forcing

  ! Container for the saved data
  ! We do NOT save the whole Type(ChmState), only the adjointed array.
  TYPE :: State_Snapshot
     REAL(fp), ALLOCATABLE :: Conc(:,:,:)   ! make sure precission is always correct
     ! Add other fields here if needed (e.g., Met_Pressure, Temperature)
  END TYPE State_Snapshot

CONTAINS

  ! -------------------------------------------------------------------
  ! PUSH: Saves a deep copy of the current state (Species(N)%Conc)
  ! -------------------------------------------------------------------
  SUBROUTINE Push_State( Current_State, Backup,NX,NY,NZ,NFD)
    TYPE(ChmState),     INTENT(IN)    :: Current_State
    TYPE(State_Snapshot), INTENT(OUT) :: Backup
    INTEGER, INTENT(IN)               :: NX,NY,NZ,NFD

     IF ( NFD < 1 .OR. NFD > Current_State%nSpecies ) THEN
       PRINT *, 'Error in Push_State: invalid NFD = ', NFD,                 &
             ' valid range is 1..', Current_State%nSpecies
       STOP
     ENDIF


    ! 1. Safety check: Ensure we aren't overwriting an existing backup
    IF ( ALLOCATED(Backup%Conc) ) THEN
       DEALLOCATE( Backup%Conc )
    END IF

    ! 2. Allocate Concentration
    ALLOCATE( Backup%Conc(NX,NY,NZ) )

    ! 3. Deep Copy the values
    Backup%Conc = Current_State%Species(NFD)%Conc

  END SUBROUTINE Push_State

  ! -------------------------------------------------------------------
  ! POP: Restores the state and frees memory
  ! -------------------------------------------------------------------
  SUBROUTINE Pop_State( Backup, Target_State,NFD )
    TYPE(State_Snapshot), INTENT(INOUT) :: Backup
    TYPE(ChmState),       INTENT(INOUT) :: Target_State
    INTEGER, INTENT(IN)               :: NFD

     IF ( NFD < 1 .OR. NFD > Target_State%nSpecies ) THEN
       PRINT *, 'Error in Pop_State: invalid NFD = ', NFD,                  &
             ' valid range is 1..', Target_State%nSpecies
       STOP
     ENDIF

    ! 1. Restore the values to the State object
    IF ( ALLOCATED(Backup%Conc) ) THEN
       Target_State%Species(NFD)%Conc = Backup%Conc
    ELSE
       ! Error handling if you pop an empty stack
       PRINT *, "Error: Attempted to Pop an unallocated State Snapshot!"
       STOP
    END IF

    ! 2. Clean up memory immediately (Crucial for Adjoints)
    DEALLOCATE( Backup%Conc )

  END SUBROUTINE Pop_State


   SUBROUTINE Setup_Adjoint_ForwardPert(State_Chm,State_Grid,Input_Opt)
    !=====================================================================
    ! PURPOSE: Setup initial adjoint state (forward perturbation or adjoint 
    !          adjoint sensitivity seed) based on FD_TYPE setting
    !
      ! COMMON INPUTS:
      !  - FD_SPEC: species name (e.g. CO2)
      !  (NFD is not a user input; it is obtained by GCHP as the species
      !    index corresponding to FD_SPEC)
      !
      !  - FD_STEP: only required for forward simulation
      !      = 0 : no perturbation
      !      = 1 : perturbation of +10%
      !      = 2 : perturbation of -10%
      !      = 3 : perturbation of +5%
      !      = 4 : perturbation of -5%
      !      = 6 : perturbation of +1%
      !      = 7 : perturbation of -1%
      !
    ! FD_TYPE OPTIONS:
    ! =====================================================================
    ! 1. GLOBAL  - Perturb/set adjoint for ALL grid cells
    !    Required inputs:  FD_SPEC (species name), FD_STEP (perturbation magnitude)           
    !                      NFD will be set from FD_SPEC
    !    Optional inputs:  None
   !    Forward mode:     Scale all concentrations in all cells as defined by FD_STEP
    !    Adjoint mode:     Set SpeciesAdj = 1 everywhere
    !
   ! 2. SPOT    - Perturb/set adjoint for a SINGLE grid cell
   !    Required inputs:  FD_SPEC, FD_STEP, (IFD, JFD, LFD)
    !                      NFD will be set from FD_SPEC
    !    Optional inputs:  None
   !    Forward mode:     Scale concentration at (IFD, JFD, LFD) as defined by FD_STEP
   !    Adjoint mode:     Set SpeciesAdj(IFD, JFD, LFD, NFD) = 1, rest = 0
    !
    ! 3. LAYER   - Perturb/set adjoint for all cells in ONE vertical layer
    !    Required inputs:  FD_SPEC, FD_STEP, LFD (layer index, 1-based)
    !                      NFD will be set from FD_SPEC
    !    Optional inputs:  None
   !    Forward mode:     Scale all concentrations at layer LFD as defined by FD_STEP
    !    Adjoint mode:     Set SpeciesAdj(:, :, LFD, NFD) = 1, rest = 0
    !
   ! 4. REGIONAL - Perturb/set adjoint for a RECTANGULAR geographic region
   !    Required inputs:  FD_SPEC, FD_STEP, FD_LAT_MIN, FD_LAT_MAX, FD_LON_MIN, FD_LON_MAX
    !                      NFD will be set from FD_SPEC
    !    Optional inputs:  LFD (if 0 or -999, uses all layers)
   !    Coordinate range: FD_LAT_MIN/MAX (latitude)  ∈ [-90, 90]
   !                      FD_LON_MIN/MAX (longitude) ∈ [-180, 180]
   !    Forward mode:     Scale concentrations where YMID ∈ [FD_LAT_MIN,FD_LAT_MAX]
   !                      AND XMID ∈ [FD_LON_MIN,FD_LON_MAX] as defined by FD_STEP
    !    Adjoint mode:     Set SpeciesAdj = 1 for cells within region, 0 outside
    ! =====================================================================
  
   TYPE(ChmState), INTENT(INOUT) ::  State_Chm
   TYPE(GrdState), INTENT(IN)    ::  State_Grid
   TYPE(OptInput), INTENT(INOUT) ::  Input_Opt
  
    ! local vars
   INTEGER :: IFD,JFD,LFD,NFD
   INTEGER :: I,J,L,L_START,L_END,N_FD_MODES
   REAL*8 :: CFN,Scale_Factor
   CHARACTER(len=ESMF_MAXSTR) :: FD_SPEC,Msg
   LOGICAL :: Is_Adj,Is_Root,Has_FD_Mode
  
    
   NFD     = Input_Opt%NFD
   IFD     = Input_Opt%IFD
   JFD     = Input_Opt%JFD
   LFD     = Input_Opt%LFD
   Is_Adj  = Input_Opt%IS_ADJOINT
   Is_Root  = Input_Opt%amIRoot
   ! IF not Had_FD_Mode, nothing needs to be changed
   Has_FD_Mode = Input_Opt%IS_FD_SPOT .OR. Input_Opt%IS_FD_GLOBAL .OR.   &
                           Input_Opt%IS_FD_LAYER .OR. Input_Opt%IS_FD_REGIONAL

   N_FD_MODES = 0
   IF ( Input_Opt%IS_FD_SPOT     ) N_FD_MODES = N_FD_MODES + 1
   IF ( Input_Opt%IS_FD_GLOBAL   ) N_FD_MODES = N_FD_MODES + 1
   IF ( Input_Opt%IS_FD_LAYER    ) N_FD_MODES = N_FD_MODES + 1
   IF ( Input_Opt%IS_FD_REGIONAL ) N_FD_MODES = N_FD_MODES + 1

   IF ( N_FD_MODES > 1 ) THEN
      WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: multiple FD modes are active.'
      WRITE(*,*) '   Exactly one of SPOT/GLOBAL/LAYER/REGIONAL must be TRUE.'
      WRITE(*,*) '   SPOT=', Input_Opt%IS_FD_SPOT,                           &
                 ' GLOBAL=', Input_Opt%IS_FD_GLOBAL,                         &
                 ' LAYER=', Input_Opt%IS_FD_LAYER,                           &
                 ' REGIONAL=', Input_Opt%IS_FD_REGIONAL
      STOP
   ENDIF
  
    Scale_Factor = 1.0d0
    Msg          = 'Not perturbing'

    ! ===== GENERAL VALIDATION =====
      IF ( Has_FD_Mode .AND. ( NFD < 1 .OR. NFD > State_Chm%nSpecies ) ) THEN
      WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: invalid NFD = ', NFD,      &
              ' valid range is 1..', State_Chm%nSpecies
       STOP
    ENDIF

     ! Validation for SPOT and LAYER: LFD must be valid
     IF ( ( Input_Opt%IS_FD_SPOT .OR. Input_Opt%IS_FD_LAYER ) .AND.         &
        ( LFD < 1 .OR. LFD > SIZE(State_Chm%SpeciesAdj,3) ) ) THEN
      WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: invalid LFD = ', LFD,      &
              ' valid range is 1..', SIZE(State_Chm%SpeciesAdj,3)
       WRITE(*,*) '   FD_TYPE=SPOT or LAYER requires LFD to be set to a valid layer index'
       STOP
     ENDIF

       ! Validation for SPOT: IFD and JFD must be valid on this PET
     IF ( Input_Opt%IS_FD_SPOT_THIS_PET .AND. Input_Opt%IS_FD_SPOT ) THEN
          IF ( IFD < 1 .OR. IFD > SIZE(State_Chm%SpeciesAdj,1) .OR.            &
               JFD < 1 .OR. JFD > SIZE(State_Chm%SpeciesAdj,2) ) THEN
             WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: invalid IFD/JFD on this PET:', &
                        IFD, JFD, ' valid I range 1..', SIZE(State_Chm%SpeciesAdj,1), &
                ' J range 1..', SIZE(State_Chm%SpeciesAdj,2)
             WRITE(*,*) '   FD_TYPE=SPOT requires IFD, JFD, LFD to be specified in GCHP.rc'
         STOP
       ENDIF
     ENDIF
    
   IF ( Has_FD_Mode .AND. .NOT. Is_Adj ) THEN
        SELECT CASE (Input_Opt%FD_STEP)
        CASE (0)
            ! No change
        CASE (1)
            Scale_Factor = 1.10d0
            Msg = 'Perturbing +0.1'
        CASE (2)
            Scale_Factor = 0.90d0
            Msg = 'Perturbing -0.1'
        CASE (3)
            Scale_Factor = 1.05d0
            Msg = 'Perturbing +0.05'
        CASE (4)
            Scale_Factor = 0.95d0
            Msg = 'Perturbing -0.05'
        CASE (6)
            Scale_Factor = 1.01d0
            Msg = 'Perturbing +0.01'
        CASE (7)
            Scale_Factor = 0.99d0
            Msg = 'Perturbing -0.01'
        CASE DEFAULT
            WRITE(*,*) ' FD_STEP = ', Input_Opt%FD_STEP, ' NOT SUPPORTED!'
        END SELECT
    END IF
  

   !  FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)


    ! ------------------------------------------------------------------
    ! LOGIC BLOCK: SPOT PERTURBATION
    ! ------------------------------------------------------------------
   ! Required: IFD, JFD, LFD (set via gchp_chunk_mod config reading)
    ! Check: Ensure inputs are valid before proceeding
    ! ------------------------------------------------------------------
  
    IF (Input_Opt%IS_FD_SPOT_THIS_PET .and.  Input_Opt%IS_FD_SPOT) THEN
       
       IF (Is_Root) THEN
          WRITE(*,*) '======== SETUP_ADJOINT_FORWARDPERT: FD_TYPE=SPOT ========'
          WRITE(*,*) 'Spot location: IFD=', IFD, ' JFD=', JFD, ' LFD=', LFD
       ENDIF
          
       IF (Is_Adj) THEN    
           State_Chm%SpeciesAdj(:,:,:,:) = 0.d0
            State_Chm%SpeciesAdj(IFD,JFD,LFD,NFD) = Input_Opt%Adjoint_Val
           IF (Is_Root) & 
            WRITE(*,*) ' Setting Single Forcing to 1 (ifd,jfd,lfd,nfd)',IFD,JFD,LFD,NFD
       ELSE
       
       ! Forward: Apply scale factor to a single point
          WRITE(*,*) TRIM(Msg)
                  State_Chm%Species(NFD)%Conc(IFD,JFD,LFD) = &
                        State_Chm%Species(NFD)%Conc(IFD,JFD,LFD) * Scale_Factor
       ENDIF
    ENDIF

  ! ------------------------------------------------------------------
  ! LOGIC BLOCK: GLOBAL PERTURBATION
  ! ------------------------------------------------------------------
  ! Required: FD_SPEC, FD_STEP
  ! Check: Ensure FD_SPEC is set (NFD > 0)
  ! ------------------------------------------------------------------
  
    IF (Input_Opt%IS_FD_GLOBAL) THEN

       IF (Is_Root) THEN
          WRITE(*,*) '======== SETUP_ADJOINT_FORWARDPERT: FD_TYPE=GLOBAL ========'
          IF (Is_Adj) THEN
             WRITE(*,*) 'Adjoint mode: setting all cells to 1'
          ELSE
             WRITE(*,*) 'Forward mode: perturbing all cells by ', TRIM(Msg)
          ENDIF
       ENDIF

       IF (Is_Adj) THEN    
           State_Chm%SpeciesAdj(:,:,:,:)   = 0.d0
           State_Chm%SpeciesAdj(:,:,:,NFD) = Input_Opt%Adjoint_Val
           IF (Is_Root) WRITE(*,*) ' Setting Global Adjoint Forcing to 1'
       ELSE
       
       ! Forward: Apply scale factor to all cells
          WRITE(*,*) TRIM(Msg)
            State_Chm%Species(NFD)%Conc(:,:,:) = &
                State_Chm%Species(NFD)%Conc(:,:,:) * Scale_Factor
       ENDIF
    ENDIF
      
  ! ------------------------------------------------------------------
  ! LOGIC BLOCK: LAYER PERTURBATION
  ! ------------------------------------------------------------------
  ! Required: FD_SPEC, FD_STEP, LFD (layer index, 1-based)
  ! Check: Ensure LFD is in valid range [1, nLevels]
  ! ------------------------------------------------------------------
    
    IF (Input_Opt%IS_FD_LAYER) THEN

       IF (LFD < 1) THEN
          WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: FD_TYPE=LAYER requires LFD > 0'
          WRITE(*,*) '   LFD found: ', LFD
          STOP
       ENDIF

       IF (Is_Root) THEN
          WRITE(*,*) '======== SETUP_ADJOINT_FORWARDPERT: FD_TYPE=LAYER ========'
          WRITE(*,*) 'Layer index: LFD=', LFD
          IF (Is_Adj) THEN
             WRITE(*,*) 'Adjoint mode: setting layer', LFD, 'to 1'
          ELSE
             WRITE(*,*) 'Forward mode: perturbing layer', LFD, 'by ', TRIM(Msg)
          ENDIF
       ENDIF

        IF (Is_Adj) THEN
            State_Chm%SpeciesAdj(:,:,:,:)     = 0.d0
            State_Chm%SpeciesAdj(:,:,LFD,NFD) = Input_Opt%Adjoint_Val
            IF (Is_Root) WRITE(*,*) ' Setting Layer ', LFD, ' Adjoint Forcing to 1'
        ELSE
            ! Forward: Apply Scale Factor to specific LAYER slice
            CALL WRITE_PARALLEL(TRIM(Msg))
            State_Chm%Species(NFD)%Conc(:,:,LFD) = &
                 State_Chm%Species(NFD)%Conc(:,:,LFD) * Scale_Factor
        END IF

    END IF

  ! ------------------------------------------------------------------
  ! LOGIC BLOCK: REGIONAL PERTURBATION
  ! ------------------------------------------------------------------
   ! Required: FD_SPEC, FD_STEP, FD_LAT_MIN, FD_LAT_MAX, FD_LON_MIN, FD_LON_MAX
  ! Optional: LFD (if 0 or -999, perturbation applies to all layers)
  ! Check:    All four bounds must be valid and within [-180,180] and [-90,90]
  ! ------------------------------------------------------------------

    IF (Input_Opt%IS_FD_REGIONAL) THEN

       IF (Input_Opt%FD_LAT_MIN == -999.0_fp .OR. Input_Opt%FD_LAT_MAX == -999.0_fp .OR. &
           Input_Opt%FD_LON_MIN == -999.0_fp .OR. Input_Opt%FD_LON_MAX == -999.0_fp) THEN
          WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: FD_TYPE=REGIONAL requires all bounds'
          WRITE(*,*) '   FD_LAT_MIN=', Input_Opt%FD_LAT_MIN, ' FD_LAT_MAX=', Input_Opt%FD_LAT_MAX
          WRITE(*,*) '   FD_LON_MIN=', Input_Opt%FD_LON_MIN, ' FD_LON_MAX=', Input_Opt%FD_LON_MAX
          STOP
       ENDIF

       ! Validate latitude bounds
       IF (Input_Opt%FD_LAT_MIN < -90.0_fp .OR. Input_Opt%FD_LAT_MIN > 90.0_fp .OR. &
           Input_Opt%FD_LAT_MAX < -90.0_fp .OR. Input_Opt%FD_LAT_MAX > 90.0_fp) THEN
          WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: FD_LAT bounds must be in [-90, 90]'
          WRITE(*,*) '   FD_LAT_MIN=', Input_Opt%FD_LAT_MIN, ' FD_LAT_MAX=', Input_Opt%FD_LAT_MAX
          STOP
       ENDIF

       ! Validate longitude bounds
       IF (Input_Opt%FD_LON_MIN < -180.0_fp .OR. Input_Opt%FD_LON_MIN > 180.0_fp .OR. &
           Input_Opt%FD_LON_MAX < -180.0_fp .OR. Input_Opt%FD_LON_MAX > 180.0_fp) THEN
          WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: FD_LON bounds must be in [-180, 180]'
          WRITE(*,*) '   FD_LON_MIN=', Input_Opt%FD_LON_MIN, ' FD_LON_MAX=', Input_Opt%FD_LON_MAX
          STOP
       ENDIF

       IF (Input_Opt%FD_LAT_MIN > Input_Opt%FD_LAT_MAX) THEN
          WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: FD_LAT_MIN must be <= FD_LAT_MAX'
          STOP
       ENDIF

       IF (Input_Opt%FD_LON_MIN > Input_Opt%FD_LON_MAX) THEN
          WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: FD_LON_MIN must be <= FD_LON_MAX'
          STOP
       ENDIF

       ! Vertical selection for REGIONAL:
       !   LFD = 0 or -999 => all levels
       !   otherwise        => only level LFD
       IF ( LFD == 0 .OR. LFD == -999 ) THEN
          L_START = 1
          L_END   = SIZE(State_Chm%SpeciesAdj,3)
       ELSEIF ( LFD >= 1 .AND. LFD <= SIZE(State_Chm%SpeciesAdj,3) ) THEN
          L_START = LFD
          L_END   = LFD
       ELSE
          WRITE(*,*) 'ERROR in Setup_Adjoint_ForwardPert: invalid LFD for FD_TYPE=REGIONAL:', LFD
          WRITE(*,*) '   Use LFD=0 or -999 for all levels, or a valid level in range 1..', &
                     SIZE(State_Chm%SpeciesAdj,3)
          STOP
       ENDIF

       IF (Is_Root) THEN
          WRITE(*,*) '======== SETUP_ADJOINT_FORWARDPERT: FD_TYPE=REGIONAL ========'
          WRITE(*,*) 'Region bounds: LAT [', Input_Opt%FD_LAT_MIN, ',', Input_Opt%FD_LAT_MAX, ']'
          WRITE(*,*) '               LON [', Input_Opt%FD_LON_MIN, ',', Input_Opt%FD_LON_MAX, ']'
          IF ( LFD == 0 .OR. LFD == -999 ) THEN
             WRITE(*,*) 'Levels: all'
          ELSE
             WRITE(*,*) 'Level: LFD=', LFD
          ENDIF
          IF (Is_Adj) THEN
             WRITE(*,*) 'Adjoint mode: setting region to 1'
          ELSE
             WRITE(*,*) 'Forward mode: perturbing region by ', TRIM(Msg)
          ENDIF
       ENDIF

        IF (Is_Adj) THEN
            ! Adjoint: Set SpeciesAdj to 1 for all cells within the specified region
            State_Chm%SpeciesAdj(:,:,:,:) = 0.d0
            
            ! Loop over all grid cells and set adjoint to 1 if within region bounds
            DO L = L_START, L_END
               DO J = 1, SIZE(State_Chm%SpeciesAdj,2)
                  DO I = 1, SIZE(State_Chm%SpeciesAdj,1)
                     ! Check if grid cell (I,J) is within regional bounds
                     ! NOTE: Exact bounds checking depends on your coordinate system
                     ! This assumes a simple rectangular region in lat/lon space
                      IF (State_Grid%YMid(I,J) >= Input_Opt%FD_LAT_MIN .AND. &
                         State_Grid%YMid(I,J) <= Input_Opt%FD_LAT_MAX .AND. &
                         State_Grid%XMid(I,J) >= Input_Opt%FD_LON_MIN .AND. &
                         State_Grid%XMid(I,J) <= Input_Opt%FD_LON_MAX) THEN
                        State_Chm%SpeciesAdj(I,J,L,NFD) = Input_Opt%Adjoint_Val
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            
            IF (Is_Root) THEN
               WRITE(*,*) ' Setting REGIONAL Adjoint Forcing to 1'
               WRITE(*,*) '   Region: LAT=[', Input_Opt%FD_LAT_MIN, ',', Input_Opt%FD_LAT_MAX, ']'
               WRITE(*,*) '           LON=[', Input_Opt%FD_LON_MIN, ',', Input_Opt%FD_LON_MAX, ']'
               WRITE(*,*) '         LEVELS=[', L_START, ',', L_END, ']'
               WRITE(*,*) '         MAX(SpeciesAdj(:,:,:,NFD))=', MAXVAL( State_Chm%SpeciesAdj(:,:,:,NFD) )
            ENDIF
        ELSE
            ! Forward: Apply Scale Factor to all grid cells within the region
            CALL WRITE_PARALLEL(TRIM(Msg)//' REGIONAL')
            DO L = L_START, L_END
               DO J = 1, SIZE(State_Chm%Species(NFD)%Conc,2)
                  DO I = 1, SIZE(State_Chm%Species(NFD)%Conc,1)
                     ! Check if grid cell (I,J) is within regional bounds
                      IF (State_Grid%YMid(I,J) >= Input_Opt%FD_LAT_MIN .AND. &
                         State_Grid%YMid(I,J) <= Input_Opt%FD_LAT_MAX .AND. &
                         State_Grid%XMid(I,J) >= Input_Opt%FD_LON_MIN .AND. &
                         State_Grid%XMid(I,J) <= Input_Opt%FD_LON_MAX) THEN
                        State_Chm%Species(NFD)%Conc(I,J,L) = &
                             State_Chm%Species(NFD)%Conc(I,J,L) * Scale_Factor
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
        ENDIF

    END IF


!
! Set all adjoint to 0
!
  IF (Is_Adj) THEN
        State_Chm%SurfaceFluxAdj(:,:,:)=0.d0
  END IF        

END SUBROUTINE Setup_Adjoint_ForwardPert
  
  
  
SUBROUTINE  Integrate_Srf_Adjoint(Input_Opt,State_Chm,State_Grid,State_Met, &
                          SurfaceFluxForward, dt) 
   USE UnitConv_Mod
   
    TYPE(OptInput),      INTENT(IN) :: Input_Opt      ! Input Options object
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm      ! Chem State object
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid     ! Grid State object
    TYPE(MetState),      INTENT(INOUT) :: State_Met      ! Met State object
   REAL(fp),            INTENT(IN) :: SurfaceFluxForward(:,:)
   REAL(fp),            INTENT(IN) :: dt                 ! Emission timestep [s]


! INTERNALS
    INTEGER                        :: previous_units
   INTEGER :: NFD, N
    INTEGER    :: RC 
    
    REAL(fp), POINTER :: rho_dry(:,:)   => NULL()
    REAL(fp), POINTER :: dz(:,:)        => NULL()
   REAL(fp) :: max_speciesadj_n, max_surfacefluxadj_n
    
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = KG_SPECIES_PER_KG_DRY_AIR,                         &
         previous_units = previous_units,                                    &
         RC             = RC                                                )
 !   _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')


       NFD=Input_Opt%NFD

       IF ( NFD < 1 .OR. NFD > State_Chm%nSpecies ) THEN
          WRITE(*,*) 'ERROR in Integrate_Srf_Adjoint: invalid NFD=', NFD,      &
                     ' valid range is 1..', State_Chm%nSpecies
          STOP
       ENDIF

       IF ( .NOT. ASSOCIATED(State_Chm%SpeciesAdj) ) THEN
          WRITE(*,*) 'ERROR in Integrate_Srf_Adjoint: SpeciesAdj is not allocated'
          STOP
       ENDIF

       IF ( .NOT. ASSOCIATED(State_Chm%SurfaceFluxAdj) ) THEN
          WRITE(*,*) 'ERROR in Integrate_Srf_Adjoint: SurfaceFluxAdj is not allocated'
          STOP
       ENDIF

       IF ( .NOT. ASSOCIATED(State_Met%AIRDEN) .OR. .NOT. ASSOCIATED(State_Met%BXHEIGHT) ) THEN
          WRITE(*,*) 'ERROR in Integrate_Srf_Adjoint: required met fields AIRDEN/BXHEIGHT are not allocated'
          STOP
       ENDIF

       IF ( NFD > SIZE(State_Chm%SurfaceFluxAdj,3) .OR. NFD > SIZE(State_Chm%SpeciesAdj,4) ) THEN
          WRITE(*,*) 'ERROR in Integrate_Srf_Adjoint: NFD exceeds array bounds. NFD=', NFD
          WRITE(*,*) '       SIZE(SurfaceFluxAdj,3)=', SIZE(State_Chm%SurfaceFluxAdj,3), &
                     ' SIZE(SpeciesAdj,4)=', SIZE(State_Chm%SpeciesAdj,4)
          STOP
       ENDIF

       ! Ensure only target species NFD carries surface adjoint values in this routine.
       DO N = 1, State_Chm%nSpecies
          IF ( N == NFD ) CYCLE

          max_speciesadj_n     = MAXVAL( ABS( State_Chm%SpeciesAdj(:,:,1,N) ) )
          max_surfacefluxadj_n = MAXVAL( ABS( State_Chm%SurfaceFluxAdj(:,:,N) ) )

          IF ( max_speciesadj_n /= 0.0_fp .OR. max_surfacefluxadj_n /= 0.0_fp ) THEN
             WRITE(*,*) 'ERROR in Integrate_Srf_Adjoint: non-target species has non-zero adjoint values'
             WRITE(*,*) '       target NFD=', NFD, ' offending N=', N
             WRITE(*,*) '       max|SpeciesAdj(:,:,1,N)|=', max_speciesadj_n
             WRITE(*,*) '       max|SurfaceFluxAdj(:,:,N)|=', max_surfacefluxadj_n
             STOP
          ENDIF
       ENDDO

    rho_dry=>State_Met%AIRDEN(:,:,1)
    dz=>State_Met%BXHEIGHT(:,:,1)

    ! Instantaneous surface-flux scaling-factor adjoint at current reverse step:
    !
    ! Forward:  dConc = (SrfFlux_base * scale) * dt / (rho_dry * dz)
    ! Adjoint:  dJ/d(scale) = SpeciesAdj(:,:,1) * SrfFlux_base * dt / (rho_dry * dz)
    !
    ! Units:  SpeciesAdj [J / (kg_spc/kg_dry)] * surf_flux [kg_spc/m2/s] * dt [s]
    !          / (rho_dry [kg_dry/m3] * dz [m])
    !         = SpeciesAdj * [kg_spc/kg_dry]  =>  dimensionless (sensitivity to scale factor)
    
        State_Chm%SurfaceFluxAdj(:,:,NFD) =      State_Chm%SurfaceFluxAdj(:,:,NFD) +                                 &
               State_Chm%SpeciesAdj(:,:,1,NFD) *                                  &
               ( SurfaceFluxForward * dt / ( rho_dry * dz ) )

        IF ( Input_Opt%amIRoot ) THEN
                WRITE(*,*) 'Integrate_Srf_Adjoint: NFD=', NFD,                    &
                   ' dt=', dt,                                                     &
                   ' max|surf_flux|=', MAXVAL( ABS( SurfaceFluxForward ) ),        &
                   ' max|rho_dry|=', MAXVAL( ABS( rho_dry ) ),                &
                   ' max|dz|=', MAXVAL( ABS( dz ) ),                          &
             ' max|SpeciesAdj(sfc)|=',                                   &
             MAXVAL( ABS( State_Chm%SpeciesAdj(:,:,1,NFD) ) ),          &
                   ' max|SpeciesAdj(global)|=',                                &
                   MAXVAL( ABS( State_Chm%SpeciesAdj(:,:,:,NFD) ) ),          &
             ' max|SurfaceFluxAdj|=',                                    &
               MAXVAL( ABS( State_Chm%SurfaceFluxAdj(:,:,NFD) ) ),       &
                   ' max|SurfaceFluxAdj(global)|=',                            &
                                 MAXVAL( ABS( State_Chm%SurfaceFluxAdj(:,:,:) ) )
        ENDIF

      CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = previous_units,                         &
         RC             = RC                                                )
 !   _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')

      
END SUBROUTINE Integrate_Srf_Adjoint


  ! -----------------------------------------------------------------------
  ! Load OCO-2 adjoint forcing from a per-obs sparse netCDF file (written
  ! by co2_adjoint_forcing.py) and ADD it to State_Chm%SpeciesAdj.
  !
  ! Rank 0 reads the file; n_obs, lat/lon arrays, and the forcing array
  ! are broadcast to all ranks.  Each rank finds the nearest local cell
  ! (minimum angular distance) for each obs and accumulates directly into
  ! SpeciesAdj, without interpolation.  Since cubed-sphere tiles are
  ! non-overlapping, each obs is accumulated by exactly one rank.
  !
  ! Forcing file: CO2_adjoint_forcing_YYYYMMDD_HHMMz.nc4
  !   Dimensions: obs (n_obs), lev (nlev)
  !   lat_obs  : (obs)      observation latitudes  [degrees_north]
  !   lon_obs  : (obs)      observation longitudes [-180,180) [degrees_east]
  !   forcing  : (obs, lev) dJ/d(CO2_mmr) [1/(kg_CO2/kg_dry_air)]
  !              (Fortran reads as forcing_f(nlev, n_obs) — column-major)
  ! -----------------------------------------------------------------------
  SUBROUTINE Load_OCO2_Adjoint_Forcing( State_Chm, State_Grid, Input_Opt, &
                                         year, month, day, hour, minute, RC )

    TYPE(ChmState),  INTENT(INOUT) :: State_Chm
    TYPE(GrdState),  INTENT(IN)    :: State_Grid
    TYPE(OptInput),  INTENT(IN)    :: Input_Opt
    INTEGER,         INTENT(IN)    :: year, month, day, hour, minute
    INTEGER,         INTENT(OUT)   :: RC

    ! Local variables
    CHARACTER(LEN=512)              :: fname
    TYPE(ESMF_VM)                   :: vm
    INTEGER                         :: STATUS
    INTEGER                         :: ncid, varid, dimid, nf_rc
    INTEGER                         :: n_obs, nlev
    INTEGER                         :: start1(1), count1(1)
    INTEGER                         :: start2(2), count2(2)
    REAL(f4), ALLOCATABLE           :: lat_obs(:), lon_obs(:)
    ! forcing_1d holds forcing_obs(nlev, n_obs) in column-major as 1D for broadcast
    REAL(f4), ALLOCATABLE, TARGET   :: forcing_1d(:)
    REAL(f4), POINTER               :: forcing_obs(:,:)
    INTEGER                         :: file_flag
    LOGICAL                         :: file_exists
    INTEGER                         :: I, J, L, K, N
    REAL(fp)                        :: obs_lat, obs_lon
    REAL(fp)                        :: obs_offset, lon_width
    ! Tile bounding box built from cell edges; obs outside it belong to another rank.
    REAL(fp)                        :: tile_lat_min, tile_lat_max

    RC = 0

    WRITE(fname,'(A,"/CO2_adjoint_forcing_",I4.4,I2.2,I2.2,"_",I2.2,I2.2,"z.nc4")') &
         TRIM(Input_Opt%OCO2_FORCING_DIR), year, month, day, hour, minute

    CALL ESMF_VMGetCurrent(vm, RC=STATUS)
    IF (STATUS /= ESMF_SUCCESS) THEN; RC = STATUS; RETURN; ENDIF

    file_flag = 0
    n_obs = 0; nlev = 0

    IF (MAPL_AM_I_ROOT(vm)) THEN
       INQUIRE(FILE=TRIM(fname), EXIST=file_exists)
       IF (file_exists) THEN
          nf_rc = nf90_open(TRIM(fname), NF90_NOWRITE, ncid)
          IF (nf_rc == NF90_NOERR) THEN
             file_flag = 1
             nf_rc = nf90_inq_dimid(ncid, 'obs', dimid)
             nf_rc = nf90_inquire_dimension(ncid, dimid, len=n_obs)
             nf_rc = nf90_inq_dimid(ncid, 'lev', dimid)
             nf_rc = nf90_inquire_dimension(ncid, dimid, len=nlev)
          ELSE
             WRITE(*,'(A,I0,1X,A)') 'Load_OCO2: nf90_open error ', nf_rc, TRIM(fname)
          ENDIF
       ELSE
          WRITE(*,*) 'Load_OCO2: no forcing file at this checkpoint: ', TRIM(fname)
       ENDIF
    ENDIF

    CALL MAPL_CommsBcast(vm, file_flag, 1, 0, RC=STATUS)
    IF (STATUS /= ESMF_SUCCESS) THEN; RC = STATUS; RETURN; ENDIF
    IF (file_flag == 0) RETURN

    CALL MAPL_CommsBcast(vm, n_obs, 1, 0, RC=STATUS)
    IF (STATUS /= ESMF_SUCCESS) THEN; RC = STATUS; RETURN; ENDIF
    CALL MAPL_CommsBcast(vm, nlev, 1, 0, RC=STATUS)
    IF (STATUS /= ESMF_SUCCESS) THEN; RC = STATUS; RETURN; ENDIF

    IF (n_obs == 0) THEN
       IF (MAPL_AM_I_ROOT(vm)) nf_rc = nf90_close(ncid)
       RETURN
    ENDIF

    ALLOCATE(lat_obs(n_obs), lon_obs(n_obs))
    ALLOCATE(forcing_1d(nlev*n_obs))
    lat_obs    = 0.0_f4
    lon_obs    = 0.0_f4
    forcing_1d = 0.0_f4
    forcing_obs => NULL()

    IF (MAPL_AM_I_ROOT(vm)) THEN
       start1 = (/1/); count1 = (/n_obs/)
       nf_rc = nf90_inq_varid(ncid, 'lat_obs', varid)
       nf_rc = nf90_get_var(ncid, varid, lat_obs, start=start1, count=count1)
       nf_rc = nf90_inq_varid(ncid, 'lon_obs', varid)
       nf_rc = nf90_get_var(ncid, varid, lon_obs, start=start1, count=count1)
       ! Python dim order ('obs','lev') → Fortran reads as forcing_f(nlev, n_obs)
       forcing_obs(1:nlev, 1:n_obs) => forcing_1d
       start2 = (/1, 1/); count2 = (/nlev, n_obs/)
       nf_rc = nf90_inq_varid(ncid, 'forcing', varid)
       nf_rc = nf90_get_var(ncid, varid, forcing_obs, start=start2, count=count2)
       IF (nf_rc /= NF90_NOERR) &
            WRITE(*,'(A,I0)') 'Load_OCO2: nf90_get_var forcing error ', nf_rc
       NULLIFY(forcing_obs)
       nf_rc = nf90_close(ncid)
    ENDIF

    CALL MAPL_CommsBcast(vm, lat_obs, n_obs, 0, RC=STATUS)
    IF (STATUS /= ESMF_SUCCESS) THEN
       DEALLOCATE(lat_obs, lon_obs, forcing_1d); RC = STATUS; RETURN
    ENDIF
    CALL MAPL_CommsBcast(vm, lon_obs, n_obs, 0, RC=STATUS)
    IF (STATUS /= ESMF_SUCCESS) THEN
       DEALLOCATE(lat_obs, lon_obs, forcing_1d); RC = STATUS; RETURN
    ENDIF
    CALL MAPL_CommsBcast(vm, forcing_1d, nlev*n_obs, 0, RC=STATUS)
    IF (STATUS /= ESMF_SUCCESS) THEN
       DEALLOCATE(lat_obs, lon_obs, forcing_1d); RC = STATUS; RETURN
    ENDIF

    ! Re-associate pointer after broadcast so loop can use (L,K) indexing
    forcing_obs(1:nlev, 1:n_obs) => forcing_1d
    N = Input_Opt%NFD

    ! State_Grid%YMin / YMax are scalars giving this rank's tile lat envelope.
    ! YEdge(NX, NY+1): cell (I,J) south edge = YEdge(I,J), north = YEdge(I,J+1)
    ! XEdge(NX+1, NY): cell (I,J) west  edge = XEdge(I,J), east  = XEdge(I+1,J)
    tile_lat_min = State_Grid%YMin
    tile_lat_max = State_Grid%YMax

    ! Point-in-cell assignment: for each obs, find the unique local cell whose
    ! lat/lon bounds contain it and accumulate directly.
    ! Because cubed-sphere cells tile the sphere without overlap, each obs
    ! belongs to exactly one rank — no double-counting.
    DO K = 1, n_obs
       obs_lat = REAL(lat_obs(K), fp)
       obs_lon = REAL(lon_obs(K), fp)

       ! Skip obs whose lat is outside this rank's tile entirely
       IF (obs_lat < tile_lat_min .OR. obs_lat >= tile_lat_max) CYCLE

       cell_search: DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          ! Lat containment: [YEdge(I,J), YEdge(I,J+1))
          IF (obs_lat < State_Grid%YEdge(I,J) .OR. &
              obs_lat >= State_Grid%YEdge(I,J+1)) CYCLE

          ! Lon containment: normalize obs offset from west edge to [0,360)
          obs_offset = obs_lon - State_Grid%XEdge(I,J)
          IF (obs_offset <    0.0_fp) obs_offset = obs_offset + 360.0_fp
          IF (obs_offset >= 360.0_fp) obs_offset = obs_offset - 360.0_fp
          lon_width  = State_Grid%XEdge(I+1,J) - State_Grid%XEdge(I,J)
          IF (lon_width <= 0.0_fp) lon_width = lon_width + 360.0_fp
          IF (obs_offset > lon_width) CYCLE

          ! Obs is inside cell (I,J): accumulate and move to next obs
          DO L = 1, MIN(State_Grid%NZ, nlev)
             State_Chm%SpeciesAdj(I,J,L,N) = State_Chm%SpeciesAdj(I,J,L,N) + &
                  REAL(forcing_obs(L,K), fp)
          ENDDO
          EXIT cell_search
       ENDDO
       ENDDO cell_search
    ENDDO

    NULLIFY(forcing_obs)
    DEALLOCATE(lat_obs, lon_obs, forcing_1d)

  END SUBROUTINE Load_OCO2_Adjoint_Forcing


END MODULE Adjoint_Utils_Mod
