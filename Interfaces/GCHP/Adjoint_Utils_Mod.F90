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

  IMPLICIT NONE
  PRIVATE

  ! Public Interface
  PUBLIC :: State_Snapshot
  PUBLIC :: Push_State
  PUBLIC :: Pop_State
  PUBLIC :: Setup_Adjoint_State
  PUBLIC :: Integrate_Srf_Adjoint

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


  SUBROUTINE Setup_Adjoint_State(State_Chm,Input_Opt)
  
    TYPE(ChmState), INTENT(INOUT) ::  State_Chm
    TYPE(OptInput), INTENT(INOUT) ::  Input_Opt
  
   ! local vars
    INTEGER :: IFD,JFD,LFD,NFD
    INTEGER :: I,J,L
    REAL*8 :: CFN,Scale_Factor
    CHARACTER(len=ESMF_MAXSTR) :: FD_SPEC,Msg
    LOGICAL :: Is_Adj,Is_Root
  
    
    NFD     = Input_Opt%NFD
    IFD     = Input_Opt%IFD
    JFD     = Input_Opt%JFD
    LFD     = Input_Opt%LFD
    Is_Adj  = Input_Opt%IS_ADJOINT
    Is_Root  = Input_Opt%amIRoot
  
    Scale_Factor = 1.0d0
    Msg          = 'Not perturbing'
    
    IF (.NOT. Is_Adj) THEN
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
        CASE DEFAULT
            WRITE(*,*) ' FD_STEP = ', Input_Opt%FD_STEP, ' NOT SUPPORTED!'
        END SELECT
    END IF
  

   !  FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)


    ! ------------------------------------------------------------------
    ! LOGIC BLOCK: SPOT PERTURBATION
    ! ------------------------------------------------------------------
  
    IF (Input_Opt%IS_FD_SPOT_THIS_PET .and.  Input_Opt%IS_FD_SPOT) THEN
          
          
       IF (Is_Adj) THEN    
           State_Chm%SpeciesAdj(:,:,:,:) = 0.d0
           State_Chm%SpeciesAdj(IFD,JFD,LFD,NFD)=1.0d0
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
  
    IF (Input_Opt%IS_FD_GLOBAL) THEN
 
    IF (Is_Adj) THEN    
           State_Chm%SpeciesAdj(:,:,:,:)   = 0.d0
           State_Chm%SpeciesAdj(:,:,:,NFD) = 1.0d0
           IF (Is_Root) WRITE(*,*) ' Setting Global Adjoint Forcing to 1'
       ELSE
       
       ! Forward: Apply scale factor to a single point
          WRITE(*,*) TRIM(Msg)
            State_Chm%Species(NFD)%Conc(:,:,:) = &
                State_Chm%Species(NFD)%Conc(:,:,:) * Scale_Factor
       ENDIF
    ENDIF
      
  ! ------------------------------------------------------------------
  ! LOGIC BLOCK: LAYER PERTURBATION
  ! ------------------------------------------------------------------
    
    IF (Input_Opt%IS_FD_LAYER) THEN

        IF (Is_Adj) THEN
            State_Chm%SpeciesAdj(:,:,:,:)     = 0.d0
            State_Chm%SpeciesAdj(:,:,LFD,NFD) = 1.0d0
            IF (Is_Root) WRITE(*,*) ' Setting Layer ', LFD, ' Adjoint Forcing to 1'
        ELSE
            ! Forward: Apply Scale Factor to specific LAYER slice
            CALL WRITE_PARALLEL(TRIM(Msg))
            State_Chm%Species(NFD)%Conc(:,:,LFD) = &
                 State_Chm%Species(NFD)%Conc(:,:,LFD) * Scale_Factor
        END IF

    END IF


!
! Set all adjoint to 0
!
  IF (Is_Adj) THEN
        State_Chm%SurfaceFluxAdj(:,:,:)=0.d0
  END IF        

END SUBROUTINE Setup_Adjoint_State
  
  
  
SUBROUTINE  Integrate_Srf_Adjoint(Input_Opt,State_Chm,State_Grid,State_Met,DT) 
   USE UnitConv_Mod
   
    TYPE(OptInput),      INTENT(IN) :: Input_Opt      ! Input Options object
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm      ! Chem State object
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid     ! Grid State object
    TYPE(MetState),      INTENT(INOUT) :: State_Met      ! Met State object
    REAL*8       :: DT

! INTERNALS
    INTEGER                        :: previous_units
    INTEGER :: NFD
    INTEGER    :: RC 
    
    REAL(fp), POINTER :: surf_flux(:,:) => NULL()
    REAL(fp), POINTER :: rho_dry(:,:)   => NULL()
    REAL(fp), POINTER :: dz(:,:)        => NULL()
    
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
    
    ! Accumulate gradients for surface flux scaling factors:
    ! The sensitivity at the previous timestep (n-1) is the sum of:
    !   1. The forcing from current concentrations (mapped to flux space via F)
    !   2. The accumulated sensitivity from all future timesteps (n..N)
    !
    ! Eq: Adj_SrfFlux = Adj_SrfFlux + (Adj_Conc * Sensitivity_Factor)
    
    surf_flux=>State_Chm%SurfaceFlux(:,:,NFD)
    rho_dry=>State_Met%AIRDEN(:,:,1)
    dz=>State_Met%BXHEIGHT(:,:,1)

    ! Accumulate adjoint of surface flux scaling factor:
    !
    ! Forward:  dConc = (SrfFlux_base * scale) * DT / (rho_dry * dz)
    ! Adjoint:  dJ/d(scale) += SpeciesAdj(:,:,1) * SrfFlux_base * DT / (rho_dry * dz)
    !
    ! Units:  SpeciesAdj [J / (kg_spc/kg_dry)] * surf_flux [kg_spc/m2/s]
    !         * DT [s] / (rho_dry [kg_dry/m3] * dz [m])
    !         = SpeciesAdj * [kg_spc/kg_dry]  =>  dimensionless (sensitivity to scale factor)
    State_Chm%SurfaceFluxAdj(:,:,NFD) = State_Chm%SurfaceFluxAdj(:,:,NFD) + &
         State_Chm%SpeciesAdj(:,:,1,NFD) * ( surf_flux * DT / ( rho_dry * dz ) )

      CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = previous_units,                         &
         RC             = RC                                                )
 !   _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')

      
END SUBROUTINE Integrate_Srf_Adjoint


END MODULE Adjoint_Utils_Mod
