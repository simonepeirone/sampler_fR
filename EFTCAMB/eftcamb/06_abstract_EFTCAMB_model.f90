!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2017 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 06_abstract_EFTCAMB_model.f90
!! This file contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB. All EFTCAMB models should inherit from this class or the two derived
!! classes contained in 06p2_abstract_EFTCAMB_full.f90 or 06p3_abstract_EFTCAMB_designer.f90


!----------------------------------------------------------------------------------------
!> This module contains the abstract definition of all the places where EFTCAMB interacts
!! with CAMB. All EFTCAMB models should inherit from this class or the two derived
!! classes contained in 06p2_abstract_EFTCAMB_full.f90 or 06p3_abstract_EFTCAMB_designer.f90

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_abstract_model

    use precision
    use IniFile
    use EFTCAMB_cache

    implicit none

    private

    public EFTCAMB_model

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for EFTCAMB models. As a rule, when there is a
    !! new model it should be declared as a class inheriting from EFTCAMB_model.
    !! This guarantees maximum performances as well as maximum flexibility.
    type, abstract :: EFTCAMB_model

        integer                       :: parameter_number !< number of parameters of the model.
        character(len=:), allocatable :: name             !< name of the model.
        character(len=:), allocatable :: name_latex       !< latex name of the model.

    contains

        ! initialization of the model:
        procedure :: init                        => EFTCAMBModelInitialize                                !< subroutine that initializes the name and latex name of the model.
        procedure(EFTCAMBModelReadModelSelectionFromFile  ), deferred :: read_model_selection             !< subroutine that reads the parameters of the model from file.
        procedure(EFTCAMBModelAllocateModelSelection      ), deferred :: allocate_model_selection         !< subroutine that allocates the model selection.
        procedure(EFTCAMBModelInitModelParameters         ), deferred :: init_model_parameters            !< subroutine taht initializes the model parameters based on the values found in an input array.
        procedure(EFTCAMBModelInitModelParametersFromFile ), deferred :: init_model_parameters_from_file  !< subroutine that reads the parameters of the model from file.

        ! utility functions:
        procedure(EFTCAMBModelComputeParametersNumber     ), deferred :: compute_param_number             !< subroutine that computes the number of parameters of the model.
        procedure(EFTCAMBModelFeedback                    ), deferred :: feedback                         !< subroutine that prints on the screen feedback information about the model.
        procedure(EFTCAMBModelParameterNames              ), deferred :: parameter_names                  !< subroutine that returns the i-th parameter name of the model.
        procedure(EFTCAMBModelParameterNamesLatex         ), deferred :: parameter_names_latex            !< subroutine that returns the i-th parameter name of the model.
        procedure(EFTCAMBModelParameterValues             ), deferred :: parameter_values                 !< subroutine that returns the i-th parameter value.

        ! background initialization functions:
        procedure :: initialize_background       => EFTCAMBModelInitBackground                            !< subroutine that initializes the background of the model, if needed.

        ! CAMB related procedures:
        procedure(EFTCAMBModelBackgroundEFTFunctions ), deferred :: compute_background_EFT_functions      !< subroutine that computes the value of the background EFT functions at a given time.
        procedure(EFTCAMBModelSecondOrderEFTFunctions), deferred :: compute_secondorder_EFT_functions     !< subroutine that computes the value of the second order EFT functions at a given time.
        procedure(EFTCAMBModelComputeDtauda          ), deferred :: compute_dtauda                        !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure(EFTCAMBModelComputeAdotoa          ), deferred :: compute_adotoa                        !< subroutine that computes adotoa = H.
        procedure(EFTCAMBModelComputeHubbleDer       ), deferred :: compute_H_derivs                      !< subroutine that computes the two derivatives wrt conformal time of H.

        procedure :: compute_rhoQPQ              => EFTCAMBModelComputeRhoQPQ                             !< subroutine that computes \rho_Q and P_Q. For details refer to the numerical notes.
        procedure :: compute_Einstein_factors    => EFTCAMBModelComputeEinsteinFactors                    !< subroutine that computes the Einstein equations factors. For details refer to the numerical notes.
        procedure :: compute_pi_factors          => EFTCAMBModelComputePiFactors                          !< subroutine that computes the pi field equations factors. For details refer to the numerical notes.
        procedure :: compute_tensor_factors      => EFTCAMBModelComputeTensorFactors                      !< subroutine that computes the factors for the tensor propagation equation. For details refer to the numerical notes.
        procedure :: compute_stability_factors   => EFTCAMBModelComputeStabilityFactors                   !< subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.

        ! stability procedures:
        procedure :: additional_model_stability  => EFTCAMBModelAdditionalModelStability                  !< function that computes model specific stability requirements.

    end type EFTCAMB_model

    ! ---------------------------------------------------------------------------------------------
    ! EFTCAMB abstract interfaces: these are all the model procedures that the user HAS to override
    ! when writing its own model.
    ! ---------------------------------------------------------------------------------------------

    abstract interface

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes the number of parameters of the model.
        subroutine EFTCAMBModelComputeParametersNumber( self )
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)       :: self   !< the base class
        end subroutine EFTCAMBModelComputeParametersNumber

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine EFTCAMBModelReadModelSelectionFromFile( self, Ini )
            use IniFile
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)  :: self   !< the base class
            type(TIniFile)        :: Ini    !< Input ini file
        end subroutine EFTCAMBModelReadModelSelectionFromFile

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine EFTCAMBModelAllocateModelSelection( self )
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)  :: self   !< the base class
        end subroutine EFTCAMBModelAllocateModelSelection

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine EFTCAMBModelInitModelParameters( self, array )
            use precision
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)                                   :: self   !< the base class
            real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.
        end subroutine EFTCAMBModelInitModelParameters

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine EFTCAMBModelInitModelParametersFromFile( self, Ini )
            use IniFile
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)  :: self   !< the base class
            type(TIniFile)        :: Ini    !< Input ini file
        end subroutine EFTCAMBModelInitModelParametersFromFile

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that prints on the screen feedback information about the model.
        subroutine EFTCAMBModelFeedback( self, print_params )
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)  :: self         !< the base class
            logical, optional     :: print_params !< optional flag that decised whether to print numerical values
                                                  !! of the parameters.
        end subroutine EFTCAMBModelFeedback

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that returns the i-th parameter name of the model.
        subroutine EFTCAMBModelParameterNames( self, i, name )
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)      :: self   !< the base class
            integer     , intent(in)  :: i      !< the index of the parameter
            character(*), intent(out) :: name   !< the output name of the i-th parameter
        end subroutine EFTCAMBModelParameterNames

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that returns the i-th parameter name of the model.
        subroutine EFTCAMBModelParameterNamesLatex( self, i, latexname )
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)      :: self       !< the base class
            integer     , intent(in)  :: i          !< The index of the parameter
            character(*), intent(out) :: latexname  !< the output latex name of the i-th parameter
        end subroutine EFTCAMBModelParameterNamesLatex

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that returns the i-th parameter name of the model.
        subroutine EFTCAMBModelParameterValues( self, i, value )
            use precision
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)  :: self   !< the base class
            integer , intent(in)  :: i      !< The index of the parameter
            real(dl), intent(out) :: value  !< the output value of the i-th parameter
        end subroutine EFTCAMBModelParameterValues

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes the value of the background EFT functions at a given time.
        subroutine EFTCAMBModelBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )
            use precision
            use EFTCAMB_cache
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                         :: a             !< the input scale factor.
            type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
            type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
        end subroutine EFTCAMBModelBackgroundEFTFunctions

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes the value of the second order EFT functions at a given time.
        subroutine EFTCAMBModelSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )
            use precision
            use EFTCAMB_cache
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                         :: a             !< the input scale factor.
            type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
            type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
        end subroutine EFTCAMBModelSecondOrderEFTFunctions

        ! ---------------------------------------------------------------------------------------------
        !> Function that computes dtauda = 1/sqrt(a^2H^2).
        function EFTCAMBModelComputeDtauda( self, a, eft_par_cache, eft_cache )
            use precision
            use EFTCAMB_cache
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                         :: a             !< the input scale factor.
            type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
            type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
            real(dl)                                     :: EFTCAMBModelComputeDtauda !< the output dtauda
        end function EFTCAMBModelComputeDtauda

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes adotoa = H and its two derivatives wrt conformal time.
        subroutine EFTCAMBModelComputeAdotoa( self, a, eft_par_cache, eft_cache )
            use precision
            use EFTCAMB_cache
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                         :: a             !< the input scale factor.
            type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
            type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
        end subroutine EFTCAMBModelComputeAdotoa

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes the two derivatives wrt conformal time of H.
        subroutine EFTCAMBModelComputeHubbleDer( self, a, eft_par_cache, eft_cache )
            use precision
            use EFTCAMB_cache
            import EFTCAMB_model
            implicit none
            class(EFTCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                         :: a             !< the input scale factor.
            type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
            type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
        end subroutine EFTCAMBModelComputeHubbleDer

    ! ---------------------------------------------------------------------------------------------

    end interface

contains

    ! ---------------------------------------------------------------------------------------------
    ! EFTCAMB abstract model implementation: the following are all the procedures that can be
    ! be safely implemented for the abstract class and are not harmful if not overritten.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the name and latex name of the model.
    subroutine EFTCAMBModelInitialize( self, name, latexname )

        implicit none

        class(EFTCAMB_model)     :: self      !< the base class
        character(*), intent(in) :: name      !< the name of the function
        character(*), intent(in) :: latexname !< the latex name of the function

        self%name       = TRIM(name)
        self%name_latex = TRIM(latexname)

    end subroutine EFTCAMBModelInitialize

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of the model, if needed.
    subroutine EFTCAMBModelInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(EFTCAMB_model)                         :: self           !< the base class
        type(EFTCAMB_parameter_cache), intent(in)    :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        success = .True.

    end subroutine EFTCAMBModelInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes \rho_Q and P_Q. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeRhoQPQ( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_model)                         :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: a2, adotoa2, aomegaP

        ! precompute some parts:
        a2      = a*a
        adotoa2 = eft_cache%adotoa**2
        aomegaP = a*eft_cache%EFTOmegaP

        ! do the computations:
        eft_cache%grhoq     = 2._dl*eft_cache%EFTc -eft_cache%EFTLambda -3._dl*adotoa2*aomegaP
        eft_cache%gpresq    = eft_cache%EFTLambda + a2*adotoa2*eft_cache%EFTOmegaPP +aomegaP*(eft_cache%Hdot+2._dl*adotoa2)
        eft_cache%grhodotq  = 3._dl*eft_cache%adotoa*(-eft_cache%grhoq-eft_cache%gpresq+adotoa2*aomegaP )
        eft_cache%gpresdotq = eft_cache%EFTLambdadot &
            & +adotoa2*eft_cache%adotoa*(a*a2*eft_cache%EFTOmegaPPP-2._dl*aomegaP+2._dl*a2*eft_cache%EFTOmegaPP) &
            & +aomegaP*eft_cache%Hdotdot &
            & +3._dl*eft_cache%adotoa*eft_cache%Hdot*( aomegaP+a2*eft_cache%EFTOmegaPP )

    end subroutine EFTCAMBModelComputeRhoQPQ

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the Einstein equations factors. For details refer to the numerical notes.
    !! The implementation might look a bit intricated because this is a crucial part performance wise.
    subroutine EFTCAMBModelComputeEinsteinFactors( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_model)                         :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: one_plus_omega, adotoa2, aomegaP, k2, a2, hdot_m_adotoa2, pidot_p_H_pi

        ! precompute some common parts:
        one_plus_omega = 1._dl+eft_cache%EFTOmegaV
        adotoa2        = eft_cache%adotoa**2
        hdot_m_adotoa2 = eft_cache%Hdot-adotoa2
        aomegaP        = a*eft_cache%EFTOmegaP
        k2             = eft_cache%k**2
        a2             = a**2
        pidot_p_H_pi   = eft_cache%pidot+eft_cache%adotoa*eft_cache%pi

        ! compute the coefficients:
        eft_cache%EFTeomF     = 1.5_dl/(eft_cache%k*one_plus_omega)*( (eft_cache%grhoq+eft_cache%gpresq)*eft_cache%pi &
            & + (aomegaP*eft_cache%adotoa+a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)*pidot_p_H_pi &
            & + eft_cache%pi*( k2*(eft_cache%EFTGamma3V+eft_cache%EFTGamma4V) -(3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V)*hdot_m_adotoa2  ) )
        eft_cache%EFTeomG     = +1._dl +0.5/one_plus_omega*( aomegaP &
            & +a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V/eft_cache%adotoa +3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V )
        eft_cache%EFTeomL     = +0.5_dl/one_plus_omega*( (2._dl*eft_cache%EFTc*pidot_p_H_pi +eft_cache%grhodotq*eft_cache%pi)/eft_cache%adotoa &
            & -3._dl*aomegaP*( (3._dl*adotoa2-eft_cache%Hdot+k2/3._dl)*eft_cache%pi +eft_cache%adotoa*eft_cache%pidot) &
            & +4._dl*a2*eft_par_cache%h0_mpc**2*eft_cache%EFTGamma1V/eft_cache%adotoa*pidot_p_H_pi &
            & +3._dl*a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V*( (eft_cache%Hdot-2._dl*adotoa2-k2/3._dl)*eft_cache%pi/eft_cache%adotoa -eft_cache%pidot ) &                                                                                  ! Gamma2
            & -eft_cache%pi*( k2-3._dl*hdot_m_adotoa2)*(3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V) +8._dl*eft_cache%EFTGamma6V*k2*pidot_p_H_pi/eft_cache%adotoa )&
            &+2._dl*eft_cache%EFTGamma5V/one_plus_omega*eft_cache%k**2*eft_cache%pi
        eft_cache%EFTeomM     = eft_cache%gpresdotq*eft_cache%pi +(eft_cache%grhoq+eft_cache%gpresq+a2*adotoa2*eft_cache%EFTOmegaPP)*pidot_p_H_pi &
            & +aomegaP*eft_cache%adotoa*( eft_cache%pidotdot +(eft_cache%Hdot+4._dl*adotoa2)*eft_cache%pidot/eft_cache%adotoa +2._dl*(eft_cache%Hdot+3._dl*adotoa2+k2/3._dl)*eft_cache%pi ) &
            & +a*eft_par_cache%h0_mpc*( eft_cache%EFTGamma2V*eft_cache%pidotdot &
            & +(4._dl*eft_cache%EFTGamma2V +a*eft_cache%EFTGamma2P)*eft_cache%adotoa*eft_cache%pidot +(3._dl*adotoa2*eft_cache%EFTGamma2V &
            & +eft_cache%Hdot*eft_cache%EFTGamma2V +a*adotoa2*eft_cache%EFTGamma2P)*eft_cache%pi) &
            & -(hdot_m_adotoa2-k2/3._dl)*( (3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V)*eft_cache%pidot &
            & +2._dl*eft_cache%adotoa*(+3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V+1.5_dl*a*eft_cache%EFTGamma3P+0.5_dl*a*eft_cache%EFTGamma4P)*eft_cache%pi ) &
            & -(3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V)*(eft_cache%Hdotdot-2._dl*eft_cache%adotoa*eft_cache%Hdot)*eft_cache%pi&
            & -4._dl*eft_cache%EFTGamma5V*k2*pidot_p_H_pi/3._dl
        eft_cache%EFTeomN     = eft_cache%k/one_plus_omega*( eft_cache%adotoa*eft_cache%pi*(-aomegaP+2._dl*eft_cache%EFTGamma4V+a*eft_cache%EFTGamma4P) &
            & +eft_cache%EFTGamma4V*eft_cache%pidot +2._dl*eft_cache%EFTGamma5V*pidot_p_H_pi )
        eft_cache%EFTeomNdot  =  eft_cache%k/one_plus_omega*( -eft_cache%Hdot*aomegaP*eft_cache%pi &
            & -eft_cache%adotoa*aomegaP*eft_cache%pidot &
            & -adotoa2*(aomegaP+a2*eft_cache%EFTOmegaPP-aomegaP**2/one_plus_omega)*eft_cache%pi &
            & +eft_cache%EFTGamma4V*eft_cache%pidotdot +a*eft_cache%adotoa*eft_cache%pidot*&
            &( +eft_cache%EFTGamma4P -eft_cache%EFTGamma4V*eft_cache%EFTOmegaP/one_plus_omega)&
            & +2._dl*(eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTGamma4P)*( eft_cache%Hdot*eft_cache%pi +eft_cache%adotoa*eft_cache%pidot)&
            & +2._dl*a*adotoa2*eft_cache%pi*(+0.5_dl*a*eft_cache%EFTGamma4PP +1.5_dl*eft_cache%EFTGamma4P&
            & -eft_cache%EFTOmegaP/one_plus_omega*(eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTGamma4P))&
            & +2._dl*eft_cache%EFTGamma5V*( eft_cache%pidotdot+eft_cache%adotoa*eft_cache%pidot+eft_cache%Hdot*eft_cache%pi)&
            & +2._dl*eft_cache%adotoa*pidot_p_H_pi*( +a*eft_cache%EFTGamma5P-eft_cache%EFTGamma5V*aomegaP/one_plus_omega) )
        eft_cache%EFTeomU     = 1._dl +(+1.5_dl*eft_cache%EFTGamma3V+0.5_dl*eft_cache%EFTGamma4V)/one_plus_omega
        eft_cache%EFTeomV     = +0.5_dl/one_plus_omega*( aomegaP -2._dl*eft_cache%EFTGamma4V -a*eft_cache%EFTGamma4P )
        eft_cache%EFTeomVdot  = 0.5_dl*eft_cache%adotoa/one_plus_omega*( aomegaP-3._dl*a*eft_cache%EFTGamma4P &
            & +a2*(eft_cache%EFTOmegaPP-eft_cache%EFTGamma4PP) +aomegaP/one_plus_omega*(-aomegaP+2._dl*eft_cache%EFTGamma4V+a*eft_cache%EFTGamma4P))
        eft_cache%EFTeomX     = 1._dl -eft_cache%EFTGamma4V/one_plus_omega
        eft_cache%EFTeomXdot  = -a*eft_cache%adotoa/one_plus_omega*( +eft_cache%EFTGamma4P &
            & -eft_cache%EFTGamma4V*eft_cache%EFTOmegaP/one_plus_omega)
        eft_cache%EFTeomY     = +0.5_dl/one_plus_omega*( aomegaP &
            & +3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V   &
            & +0.5_dl*a*(3._dl*eft_cache%EFTGamma3P+eft_cache%EFTGamma4P) )
        eft_cache%EFTeomQ =1._dl+2._dl*eft_cache%EFTGamma5V/one_plus_omega

    end subroutine EFTCAMBModelComputeEinsteinFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the pi field equations factors. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputePiFactors( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_model)                         :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: one_plus_omega, adotoa2, aomegaP, k2, a2, hdot_m_adotoa2, adotoa02

        ! precompute some common parts:
        k2             = eft_cache%k**2
        a2             = a**2
        one_plus_omega = 1._dl+eft_cache%EFTOmegaV
        adotoa2        = eft_cache%adotoa**2
        adotoa02       = eft_par_cache%h0_mpc**2
        hdot_m_adotoa2 = eft_cache%Hdot-adotoa2
        aomegaP        = a*eft_cache%EFTOmegaP

        ! compute the coefficients:
        eft_cache%EFTpiA1  = eft_cache%EFTc +2._dl*a2*adotoa02*eft_cache%EFTGamma1V +1.5_dl*a2*( eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V )**2&
            &/(2._dl*one_plus_omega +3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V)
            !
        eft_cache%EFTpiA2  = +4._dl*eft_cache%EFTGamma6V
            !
        eft_cache%EFTpiB1  = eft_cache%EFTcdot +4._dl*eft_cache%adotoa*eft_cache%EFTc +8._dl*a2*eft_cache%adotoa*adotoa02*(eft_cache%EFTGamma1V +0.25_dl*a*eft_cache%EFTGamma1P)&
            & -a*(eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/( 4._dl*one_plus_omega +6._dl*eft_cache%EFTGamma3V +2._dl*eft_cache%EFTGamma4V)*&
            &(-3._dl*( eft_cache%grhoq +eft_cache%gpresq ) -3._dl*adotoa2*aomegaP*(4._dl+eft_cache%Hdot/adotoa2) -3._dl*a2*adotoa2*eft_cache%EFTOmegaPP&
            & -3._dl*a*eft_cache%adotoa*eft_par_cache%h0_mpc*(4._dl*eft_cache%EFTGamma2V +a*eft_cache%EFTGamma2P) +( 9._dl*eft_cache%EFTGamma3V +3._dl*eft_cache%EFTGamma4V)*hdot_m_adotoa2 )&
            & +1._dl/(one_plus_omega+2._dl*eft_cache%EFTGamma5V)*( eft_cache%adotoa*aomegaP +2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P)&
            & -one_plus_omega*(eft_cache%adotoa*aomegaP +a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/( 2._dl*one_plus_omega +3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V))*&
            &(-eft_cache%EFTc +1.5_dl*a*adotoa2*eft_cache%EFTOmegaP -2._dl*a2*eft_par_cache%h0_mpc**2*eft_cache%EFTGamma1V +1.5_dl*a*eft_cache%adotoa*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)
            !
        eft_cache%EFTpiB2  = +4._dl*eft_cache%adotoa*(2._dl*eft_cache%EFTGamma6V +a*eft_cache%EFTGamma6P) +a*(eft_cache%EFTGamma4V&
            & +2._dl*eft_cache%EFTGamma5V)/(2._dl*(1._dl+eft_cache%EFTOmegaV)-2._dl*eft_cache%EFTGamma4V)*(eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V) &
            & -a*(eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/( 4._dl*one_plus_omega +6._dl*eft_cache%EFTGamma3V +2._dl*eft_cache%EFTGamma4V)*&
            &( +(-3._dl*eft_cache%EFTGamma3V -eft_cache%EFTGamma4V +4._dl*eft_cache%EFTGamma5V ))&
            & +1._dl/(1._dl+eft_cache%EFTOmegaV +2._dl*eft_cache%EFTGamma5V)*( a*eft_cache%adotoa*eft_cache%EFTOmegaP +2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V + a*eft_cache%EFTGamma5P)&
            & -(1._dl+eft_cache%EFTOmegaV)*(a*eft_cache%adotoa*eft_cache%EFTOmegaP +a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/( 2._dl*one_plus_omega +3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V))*&
            &( -4._dl*eft_cache%EFTGamma6V )
            !
        eft_cache%EFTpiC = +eft_cache%adotoa*eft_cache%EFTcdot + ( 6._dl*adotoa2 -2._dl*eft_cache%Hdot)*eft_cache%EFTc +1.5_dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP*( eft_cache%Hdotdot -2._dl*eft_cache%adotoa**3) &
            & +6._dl*(a*eft_cache%adotoa*eft_par_cache%h0_mpc)**2*eft_cache%EFTGamma1V +2._dl*a2*eft_cache%Hdot*adotoa02*eft_cache%EFTGamma1V &
            & +2._dl*a**3*adotoa2*adotoa02*eft_cache%EFTGamma1P +1.5_dl*hdot_m_adotoa2**2*(eft_cache%EFTGamma4V +3._dl*eft_cache%EFTGamma3V )&
            & +4.5_dl*eft_cache%adotoa*eft_par_cache%h0_mpc*a*hdot_m_adotoa2*( eft_cache%EFTGamma2V +a*eft_cache%EFTGamma2P/3._dl )&
            & +0.5_dl*a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V*( 3._dl*eft_cache%Hdotdot -12._dl*eft_cache%Hdot*eft_cache%adotoa +6._dl*eft_cache%adotoa**3) &
            & -a*( eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/(4._dl*(1._dl+eft_cache%EFTOmegaV)+6._dl*eft_cache%EFTGamma3V +2._dl*eft_cache%EFTGamma4V)*&
            &(-3._dl*eft_cache%gpresdotq -3._dl*eft_cache%adotoa*( eft_cache%grhoq +eft_cache%gpresq) -3._dl*a*eft_cache%adotoa**3*( a*eft_cache%EFTOmegaPP +6._dl*eft_cache%EFTOmegaP) &
            & -6._dl*a*eft_cache%adotoa*eft_cache%Hdot*eft_cache%EFTOmegaP +3._dl*(eft_cache%Hdotdot -2._dl*eft_cache%adotoa*eft_cache%Hdot)*(eft_cache%EFTGamma4V +3._dl*eft_cache%EFTGamma3V)&
            & +6._dl*eft_cache%adotoa*hdot_m_adotoa2*( 3._dl*eft_cache%EFTGamma3V +1.5_dl*a*eft_cache%EFTGamma3P +eft_cache%EFTGamma4V + 0.5_dl*a*eft_cache%EFTGamma4P)&
            & -3._dl*a*eft_par_cache%h0_mpc*(3._dl*adotoa2*eft_cache%EFTGamma2V +eft_cache%Hdot*eft_cache%EFTGamma2V +a*adotoa2*eft_cache%EFTGamma2P))&
            & +1._dl/(1._dl +eft_cache%EFTOmegaV +2._dl*eft_cache%EFTGamma5V)*( a*eft_cache%adotoa*eft_cache%EFTOmegaP +2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V +a*eft_cache%EFTGamma5P)&
            & -(1._dl+eft_cache%EFTOmegaV)*( a*eft_cache%adotoa*eft_cache%EFTOmegaP+a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/(2._dl*(1._dl+eft_cache%EFTOmegaV)+3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V))*&
            &(-0.5*eft_cache%grhodotq -eft_cache%adotoa*eft_cache%EFTc +1.5_dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP*(3._dl*adotoa2 -eft_cache%Hdot) -2._dl*a2*eft_cache%adotoa*adotoa02*eft_cache%EFTGamma1V&
            & -1.5_dl*a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V*(eft_cache%Hdot-2._dl*adotoa2) -3._dl*eft_cache%adotoa*hdot_m_adotoa2*(1.5_dl*eft_cache%EFTGamma3V +0.5_dl*eft_cache%EFTGamma4V))
            !
        eft_cache%EFTpiD1 = eft_cache%EFTc -0.5_dl*a*eft_cache%adotoa*eft_par_cache%h0_mpc*(eft_cache%EFTGamma2V +a*eft_cache%EFTGamma2P) -hdot_m_adotoa2*(3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V)&
            & +4._dl*( eft_cache%Hdot*eft_cache%EFTGamma6V + adotoa2*eft_cache%EFTGamma6V + a*adotoa2*eft_cache%EFTGamma6P)&
            & +2._dl*( eft_cache%Hdot*eft_cache%EFTGamma5V +a*adotoa2*eft_cache%EFTGamma5P)&
            & -a*(eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/(4._dl*(1._dl+eft_cache%EFTOmegaV)+6._dl*eft_cache%EFTGamma3V +2._dl*eft_cache%EFTGamma4V)*&
            &(-2._dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP +4._dl*eft_cache%adotoa*eft_cache%EFTGamma5V -2._dl*eft_cache%adotoa*(3._dl*eft_cache%EFTGamma3V +1.5_dl*a*eft_cache%EFTGamma3P &
            & +eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTGamma4P))&
            & +1._dl/(1._dl+eft_cache%EFTOmegaV+2._dl*eft_cache%EFTGamma5V)*(a*eft_cache%adotoa*eft_cache%EFTOmegaP+2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V +a*eft_cache%EFTGamma5P)&
            & -(1._dl+eft_cache%EFTOmegaV)*(a*eft_cache%adotoa*eft_cache%EFTOmegaP +a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/(2._dl*(1._dl+eft_cache%EFTOmegaV)+3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V))*&
            &(+0.5_dl*a*eft_cache%adotoa*eft_cache%EFTOmegaP -2._dl*eft_cache%adotoa*eft_cache%EFTGamma5V +0.5_dl*a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V +1.5_dl*eft_cache%adotoa*eft_cache%EFTGamma3V&
            & +0.5_dl*eft_cache%adotoa*eft_cache%EFTGamma4V -4._dl*eft_cache%adotoa*eft_cache%EFTGamma6V)&
            & +(eft_cache%EFTGamma4V +2._dl*eft_cache%EFTGamma5V)/(2._dl*(1._dl+eft_cache%EFTOmegaV) -2._dl*eft_cache%EFTGamma4V)*(eft_cache%grhoq +eft_cache%gpresq +a*adotoa2*eft_cache%EFTOmegaP&
            & -eft_cache%EFTGamma4V*hdot_m_adotoa2 +a*eft_cache%adotoa*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V -3._dl*eft_cache%EFTGamma3V*hdot_m_adotoa2)
            !
        eft_cache%EFTpiD2 = +(+0.5_dl*eft_cache%EFTGamma3V +0.5_dl*eft_cache%EFTGamma4V &
            & +(eft_cache%EFTGamma4V +2._dl*eft_cache%EFTGamma5V)/(2._dl*(1._dl+eft_cache%EFTOmegaV) -2._dl*eft_cache%EFTGamma4V)*(eft_cache%EFTGamma3V +eft_cache%EFTGamma4V))
            !
        eft_cache%EFTpiE = (eft_cache%EFTc -1.5_dl*a*adotoa2*eft_cache%EFTOmegaP -0.5_dl*a*eft_cache%adotoa*eft_par_cache%h0_mpc*(2._dl*eft_cache%EFTGamma2V +a*eft_cache%EFTGamma2P)&
            & +0.5_dl*eft_cache%EFTGamma3V*(k2 -3._dl*eft_cache%Hdot +3._dl*adotoa2) +0.5_dl*eft_cache%EFTGamma4V*(k2 -eft_cache%Hdot +adotoa2)&
            & -a*(eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/( 4._dl*(1._dl+eft_cache%EFTOmegaV) +6._dl*eft_cache%EFTGamma3V +2._dl*eft_cache%EFTGamma4V)*&
            &(-2._dl*eft_cache%adotoa*(a*eft_cache%EFTOmegaP +2._dl*(1._dl+eft_cache%EFTOmegaV)) -2._dl*eft_cache%adotoa*(3._dl*eft_cache%EFTGamma3V +1.5_dl*a*eft_cache%EFTGamma3P&
            & +eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTGamma4P))&
            & +1._dl/(1._dl+eft_cache%EFTOmegaV+2._dl*eft_cache%EFTGamma5V)*(a*eft_cache%adotoa*eft_cache%EFTOmegaP+2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V +a*eft_cache%EFTGamma5P)&
            & -(1._dl+eft_cache%EFTOmegaV)*(a*eft_cache%adotoa*eft_cache%EFTOmegaP+a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/(2._dl*(1._dl+eft_cache%EFTOmegaV)+3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V))*&
            &( +eft_cache%adotoa*(1._dl +eft_cache%EFTOmegaV +0.5_dl*a*eft_cache%EFTOmegaP) +0.5_dl*a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V +1.5_dl*eft_cache%adotoa*eft_cache%EFTGamma3V +0.5_dl*eft_cache%adotoa*eft_cache%EFTGamma4V)&
            & +(eft_cache%EFTGamma4V +2._dl*eft_cache%EFTGamma5V)/(2._dl*(1._dl +eft_cache%EFTOmegaV) -2._dl*eft_cache%EFTGamma4V)*k2*(eft_cache%EFTGamma4V +eft_cache%EFTGamma3V))*eft_cache%k*eft_cache%z&
            & +1._dl*a*(eft_cache%adotoa*eft_cache%EFTOmegaP +eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/(4._dl*(1._dl+eft_cache%EFTOmegaV)+6._dl*eft_cache%EFTGamma3V +2._dl*eft_cache%EFTGamma4V)*&
            &(eft_cache%grhog_t*eft_cache%clxg +eft_cache%grhor_t*eft_cache%clxr +3._dl*eft_cache%dgpnu ) +(eft_cache%EFTGamma4V +2._dl*eft_cache%EFTGamma5V)/(2._dl*(1._dl+eft_cache%EFTOmegaV) -2._dl*eft_cache%EFTGamma4V)*eft_cache%k*eft_cache%dgq&
            & -0.5_dl/(1._dl+eft_cache%EFTOmegaV +2._dl*eft_cache%EFTGamma5V)*(a*eft_cache%adotoa*eft_cache%EFTOmegaP+2._dl*eft_cache%adotoa*(eft_cache%EFTGamma5V +a*eft_cache%EFTGamma5P)&
            & -(1._dl+eft_cache%EFTOmegaV)*(a*eft_cache%adotoa*eft_cache%EFTOmegaP+a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)/(2._dl*(1._dl+eft_cache%EFTOmegaV)+3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V))*eft_cache%dgrho

    end subroutine EFTCAMBModelComputePiFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the factors for the tensor propagation equation. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeTensorFactors( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_model)                         :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl) :: one_plus_omega

        ! precompute some common parts:
        one_plus_omega = 1._dl+eft_cache%EFTOmegaV

        ! compute the coefficients:
        eft_cache%EFTAT = one_plus_omega -eft_cache%EFTGamma4V
        eft_cache%EFTBT = 2._dl*eft_cache%adotoa*( one_plus_omega -eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTOmegaP -0.5_dl*a*eft_cache%EFTGamma4P )
        eft_cache%EFTDT = one_plus_omega

    end subroutine EFTCAMBModelComputeTensorFactors

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeStabilityFactors( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_model)                         :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        real(dl):: F1, F2, F3, F1dot, F2dot, rho_d, rho_d_dot, F1dotdot, F2dotdot, F3dot, F3dotdot, value1, value2, EFTGamma3PP, EFTcP, EFTcPP
        real(dl):: gpinudotdot_tot, w_m, EFTOmegaPPPP,Hddd
        !(now just for GBD!!!)
        rho_d = (eft_cache%grhob_t+ eft_cache%grhoc_t)/a**2
        rho_d_dot = -3._dl*eft_cache%adotoa/a*rho_d

        !< Define some EFT functions needed by stability
        EFTGamma3PP = -eft_cache%EFTGamma4PP
        EFTcP = eft_cache%EFTcdot/(a*eft_cache%adotoa) + 2._dl*eft_cache%EFTc/a

        gpinudotdot_tot= 0._dl !SP not massive nu implemented yet!

        w_m = 1._dl/3._dl

        ! EFTOmegaPPPP = self%PureEFTOmega%fourth_derivative(a)

        ! Hddd = ((-3*eft_cache%gpinu_tot*eft_cache%Hdot + eft_cache%grhov_t*eft_cache%Hdot +6*self%PureEFTwDE%value(a)*eft_cache%grhov_t*eft_cache%Hdot &
        !     &+ 9*self%PureEFTwDE%value(a)**2*eft_cache%grhov_t*eft_cache%Hdot +( eft_cache%grhob_t +eft_cache%grhoc_t)*eft_cache%Hdot &
        !     &+ 6*w_m*( eft_cache%grhob_t +eft_cache%grhoc_t)*eft_cache%Hdot +9*w_m**2*( eft_cache%grhob_t +eft_cache%grhoc_t)*eft_cache%Hdot &
        !     &+ eft_cache%grhonu_tot*eft_cache%Hdot -3*a*eft_cache%grhov_t*eft_cache%Hdot*self%PureEFTwDE%first_derivative(a)) +eft_cache%adotoa*(&
        !     &-9*eft_cache%gpinudot_tot -3*eft_cache%adotoa*(eft_cache%gpinu_tot + eft_cache%grhonu_tot)) - 3*gpinudotdot_tot &
        !     &-eft_cache%adotoa**2*(6*eft_cache%gpinu_tot + ( eft_cache%grhob_t +eft_cache%grhoc_t) + 9*w_m*( eft_cache%grhob_t +eft_cache%grhoc_t)&
        !     & + 27*w_m**2*( eft_cache%grhob_t +eft_cache%grhoc_t) +27*w_m**3*( eft_cache%grhob_t +eft_cache%grhoc_t) - 2*eft_cache%grhonu_tot &
        !     &+eft_cache%grhov_t*(1 + 27*self%PureEFTwDE%value(a)**2 + 27*self%PureEFTwDE%value(a)**3 - 6*a*self%PureEFTwDE%first_derivative(a)) &
        !     &+self%PureEFTwDE%value(a)*(9 - 27*a*self%PureEFTwDE%first_derivative(a))) + 3*a**2*self%PureEFTwDE%second_derivative(a))/(6._dl)
        !
        !
        ! EFTcdotdot = (a**2*eft_cache%grhov_t*eft_cache%Hdot*(-3*(1 + self%PureEFTwDE%value(a))**2 + a*self%PureEFTwDE%first_derivative(a)) &
        !     &+  (eft_cache%EFTOmegaV*(8*eft_cache%Hdot**2 - 2*Hddd) - a*(-(eft_cache%Hdot**2*(eft_cache%EFTOmegaP &
        !     &- 3*a*eft_cache%EFTOmegaPP)) + eft_cache%adotoa*eft_cache%Hdotdot*(eft_cache%EFTOmegaP + a*eft_cache%EFTOmegaPP) + eft_cache%EFTOmegaP*Hddd)) &
        !     &+ a*eft_cache%adotoa*eft_cache%adotoa**3* (eft_cache%EFTOmegaP + a*(3*eft_cache%EFTOmegaPP + a*eft_cache%EFTOmegaPPP)) - eft_cache%adotoa* (4*eft_cache%EFTcdot &
        !     &+ eft_cache%Hdotdot*(-8*eft_cache%EFTOmegaV + a*(eft_cache%EFTOmegaP + 3*a*eft_cache%EFTOmegaPP)) + a*eft_cache%adotoa*eft_cache%Hdot*(-eft_cache%EFTOmegaP &
        !     &+ a*(5*eft_cache%EFTOmegaPP + 3*a*eft_cache%EFTOmegaPPP))) + eft_cache%adotoa**2*(a**2*eft_cache%grhov_t*(3 - 6*a*self%PureEFTwDE%first_derivative(a) &
        !     &+ 3*self%PureEFTwDE%value(a)*(5 + self%PureEFTwDE%value(a)*(7 + 3*self%PureEFTwDE%value(a)) - 3*a*self%PureEFTwDE%first_derivative(a)) + a**2*self%PureEFTwDE%second_derivative(a)) &
        !     &+  eft_cache%Hdot* (-12*eft_cache%EFTOmegaV + a*(11*eft_cache%EFTOmegaP + 3*a*(eft_cache%EFTOmegaPP - a*eft_cache%EFTOmegaPPP)))) &
        !     &-  a*eft_cache%adotoa**4* (4*eft_cache%EFTOmegaP + a**2*(3*eft_cache%EFTOmegaPPP + a*EFTOmegaPPPP)))/(2.)

        EFTcPP = ( eft_cache%EFTcdotdot/eft_cache%adotoa +(3._dl-eft_cache%Hdot/eft_cache%adotoa)*eft_cache%EFTcdot +2._dl*eft_cache%adotoa*eft_cache%EFTc )/(a**2*eft_cache%adotoa**2)

        ! expressions for GBD
        ! F1 = 2._dl*(1._dl+eft_cache%EFTOmegaV)
        ! F2 = eft_cache%adotoa/a*F1 + eft_cache%adotoa*eft_cache%EFTOmegaP
        ! F3 = (-6* (1 + eft_cache%EFTOmegaV)*eft_cache%adotoa**2)/a**2 - (6* eft_cache%adotoa**2*eft_cache%EFTOmegaP)/a
        ! F1dot = eft_cache%adotoa*(2._dl*eft_cache%EFTOmegaP)
        ! F2dot = ( ((2 + 2*eft_cache%EFTOmegaV + a*eft_cache%EFTOmegaP)*eft_cache%Hdot +eft_cache%adotoa**2*(-2 &
        !       &- 2*eft_cache%EFTOmegaV + 2*a*eft_cache%EFTOmegaP + a**2*eft_cache%EFTOmegaPP)))/a**2
        ! F3dot =  -(( (6*eft_cache%adotoa*(2 + 2*eft_cache%EFTOmegaV + 2*a*eft_cache%EFTOmegaP)*eft_cache%Hdot &
        !     &-3*eft_cache%adotoa**3*(4 + 4*eft_cache%EFTOmegaV - 2*a**2*eft_cache%EFTOmegaPP)))/a**3)
        ! F1dotdot = (2*eft_cache%EFTOmegaP*eft_cache%Hdot  +a*eft_cache%adotoa**2*2*eft_cache%EFTOmegaPP )/a
        ! F2dotdot =( (eft_cache%adotoa*eft_cache%Hdot*(-8 - 8*eft_cache%EFTOmegaV + 5*a*eft_cache%EFTOmegaP + 3*a**2*eft_cache%EFTOmegaPP) &
        !       &+(2 + 2*eft_cache%EFTOmegaV + a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot +eft_cache%adotoa**3*(4 + 4*eft_cache%EFTOmegaV &
        !       &- 4*a*eft_cache%EFTOmegaP + 2*a**2*eft_cache%EFTOmegaPP + a**3*eft_cache%EFTOmegaPPP)))/a**3
        ! F3dotdot = ( (-6*(2 + 2*eft_cache%EFTOmegaV + 2*a*eft_cache%EFTOmegaP)*eft_cache%Hdot**2 +3*eft_cache%adotoa**2*eft_cache%Hdot*(24 &
        !       &+ 24*eft_cache%EFTOmegaV + 4*a*eft_cache%EFTOmegaP - 10*a**2*eft_cache%EFTOmegaPP) -6*eft_cache%adotoa*(2 + 2*eft_cache%EFTOmegaV &
        !       &+ 2*a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot -3*eft_cache%adotoa**4*(12 + 12*eft_cache%EFTOmegaV - 4*a*eft_cache%EFTOmegaP &
        !       &- 2*a**2*eft_cache%EFTOmegaPP + 2*a**3*eft_cache%EFTOmegaPPP)))/a**4

        !< general expressions for Fs
        F1 = 2*(1 + eft_cache%EFTOmegaV) + 3*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V

        F2 = eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V + (eft_cache%adotoa*(2*(1 + eft_cache%EFTOmegaV) + 3*eft_cache%EFTGamma3V &
            &+ eft_cache%EFTGamma4V))/a + eft_cache%adotoa*eft_cache%EFTOmegaP

        F3 = (2*eft_cache%EFTc)/a**2 + 4*eft_par_cache%h0_Mpc**2*eft_cache%EFTGamma1V - (6*eft_par_cache%h0_Mpc*eft_cache%adotoa*eft_cache%EFTGamma2V)/a &
            &- (3*eft_cache%adotoa**2*(2*(1 + eft_cache%EFTOmegaV) + 3*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V))/a**2 &
            &- (6*eft_cache%adotoa**2*eft_cache%EFTOmegaP)/a

        F1dot = eft_cache%adotoa*(2*eft_cache%EFTOmegaP + 3*eft_cache%EFTGamma3P + eft_cache%EFTGamma4P)

        F2dot = ((2 + 2*eft_cache%EFTOmegaV + 3*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + a*eft_cache%EFTOmegaP)*eft_cache%Hdot &
            &+  a**2*eft_par_cache%h0_Mpc**2*eft_cache%adotoa*eft_cache%EFTGamma2P + eft_cache%adotoa**2*(-2 - 2*eft_cache%EFTOmegaV &
            &- 3*eft_cache%EFTGamma3V - eft_cache%EFTGamma4V + 2*a*eft_cache%EFTOmegaP + 3*a*eft_cache%EFTGamma3P + a*eft_cache%EFTGamma4P &
            &+ a**2*eft_cache%EFTOmegaPP))/a**2

        F3dot = -((4*eft_cache%EFTc*eft_cache%adotoa + 6*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*eft_cache%Hdot - 2*eft_cache%adotoa*(a*EFTcP &
            &- 3*(2 + 2*eft_cache%EFTOmegaV + 3*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + 2*a*eft_cache%EFTOmegaP)*eft_cache%Hdot &
            &+  2*a**3*eft_par_cache%h0_Mpc**2*eft_cache%EFTGamma1P) +  6*a*eft_par_cache%h0_Mpc*eft_cache%adotoa**2*(-eft_cache%EFTGamma2V &
            &+ a*eft_cache%EFTGamma2P) -  3*eft_cache%adotoa**3*(4 + 4*eft_cache%EFTOmegaV + 6*eft_cache%EFTGamma3V + 2*eft_cache%EFTGamma4V &
            &- 3*a*eft_cache%EFTGamma3P - a*eft_cache%EFTGamma4P -  2*a**2*eft_cache%EFTOmegaPP))/a**3)

        F1dotdot = (2*eft_cache%EFTOmegaP*eft_cache%Hdot + eft_cache%Hdot*(3*eft_cache%EFTGamma3P + eft_cache%EFTGamma4P) + a*eft_cache%adotoa**2*(2*eft_cache%EFTOmegaPP &
            &+ 3*EFTGamma3PP + eft_cache%EFTGamma4PP))/a

        F2dotdot = (a**2*eft_par_cache%h0_Mpc**2*eft_cache%Hdot*eft_cache%EFTGamma2P +  eft_cache%adotoa*eft_cache%Hdot*(-8 - 8*eft_cache%EFTOmegaV - 12*eft_cache%EFTGamma3V &
            &- 4*eft_cache%EFTGamma4V + 5*a*eft_cache%EFTOmegaP +  9*a*eft_cache%EFTGamma3P + 3*a*eft_cache%EFTGamma4P + 3*a**2*eft_cache%EFTOmegaPP) &
            &+ (2 + 2*eft_cache%EFTOmegaV + 3*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot &
            &+ a**3*eft_par_cache%h0_Mpc**2*eft_cache%adotoa**2*eft_cache%EFTGamma2PP +  eft_cache%adotoa**3*(4 + 4*eft_cache%EFTOmegaV &
            &+ 6*eft_cache%EFTGamma3V + 2*eft_cache%EFTGamma4V - 4*a*eft_cache%EFTOmegaP - 6*a*eft_cache%EFTGamma3P - 2*a*eft_cache%EFTGamma4P + 2*a**2*eft_cache%EFTOmegaPP &
            &+ 3*a**2*EFTGamma3PP + a**2*eft_cache%EFTGamma4PP + a**3*eft_cache%EFTOmegaPPP))/a**3

        F3dotdot = (4*eft_cache%EFTc*(3*eft_cache%adotoa**2 - eft_cache%Hdot) +  2*(a*EFTcP*eft_cache%Hdot - 3*(2 + 2*eft_cache%EFTOmegaV &
            &+ 3*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V + 2*a*eft_cache%EFTOmegaP)*eft_cache%Hdot**2 + 2*a**3*eft_par_cache%h0_Mpc**2*eft_cache%Hdot*eft_cache%EFTGamma1P &
            &-  3*a*eft_par_cache%h0_Mpc*eft_cache%EFTGamma2V*eft_cache%Hdotdot) -  6*eft_cache%adotoa*(a*eft_par_cache%h0_Mpc*eft_cache%Hdot*(&
            &-4*eft_cache%EFTGamma2V + 3*a*eft_cache%EFTGamma2P) +  (2 + 2*eft_cache%EFTOmegaV + 3*eft_cache%EFTGamma3V + eft_cache%EFTGamma4V &
            &+ 2*a*eft_cache%EFTOmegaP)*eft_cache%Hdotdot) + eft_cache%adotoa**2*(-8*a*EFTcP +  3*eft_cache%Hdot*(24&
            & + 24*eft_cache%EFTOmegaV + 36*eft_cache%EFTGamma3V + 12*eft_cache%EFTGamma4V + 4*a*eft_cache%EFTOmegaP -   15*a*eft_cache%EFTGamma3P &
            &- 5*a*eft_cache%EFTGamma4P - 10*a**2*eft_cache%EFTOmegaPP) + 2*a**2*(EFTcPP + 2*a**2*eft_par_cache%h0_Mpc**2*eft_cache%EFTGamma1PP)) &
            &- 6*a*eft_par_cache%h0_Mpc*eft_cache%adotoa**3*(2*eft_cache%EFTGamma2V + a*(-2*eft_cache%EFTGamma2P + a*eft_cache%EFTGamma2PP)) &
            &- 3*eft_cache%adotoa**4*(12 + 12*eft_cache%EFTOmegaV + 18*eft_cache%EFTGamma3V + 6*eft_cache%EFTGamma4V - 4*a*eft_cache%EFTOmegaP - 12*a*eft_cache%EFTGamma3P &
            &-  4*a*eft_cache%EFTGamma4P - 2*a**2*eft_cache%EFTOmegaPP + 3*a**2*EFTGamma3PP + a**2*eft_cache%EFTGamma4PP + 2*a**3*eft_cache%EFTOmegaPPP))/a**4

        ! !SP: old condition
        eft_cache%EFT_kinetic  = 9._dl*( 1._dl +eft_cache%EFTOmegaV -eft_cache%EFTGamma4V )*( 4._dl*eft_cache%EFTc*( 1._dl +eft_cache%EFTOmegaV -eft_cache%EFTGamma4V ) &
            & +3._dl*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*a**2 + a**2*eft_par_cache%h0_Mpc*( eft_par_cache%h0_Mpc*( 3._dl*eft_cache%EFTGamma2V**2 +8._dl*eft_cache%EFTGamma1V* &
            &( 1._dl +eft_cache%EFTOmegaV -eft_cache%EFTGamma4V ) +6._dl*eft_cache%adotoa*eft_cache%EFTGamma2V*eft_cache%EFTOmegaP ) ) )
        !
        ! !SP: new Kinetic term
        ! eft_cache%EFT_kinetic = (F1*F3 + 3._dl*F2**2)*F1*a**3 /(2._dl*F2**2)

        !SP just for check
        ! if (abs(eft_cache%EFT_kinetic)<1.d-8)eft_cache%EFT_kinetic = 0._dl
        ! eft_cache%EFT_kinetic  =(1._dl +eft_cache%EFTOmegaV)*( 4._dl*eft_cache%EFTc*(1._dl +eft_cache%EFTOmegaV) +3._dl*a**3*eft_cache%adotoa**2*eft_cache%EFTOmegaP*&
        !     &( 4._dl*(a-1._dl)*(1._dl +eft_cache%EFTOmegaV) +a*eft_cache%EFTOmegaP ))

        !SP: old condition
        eft_cache%EFT_gradient = 9._dl*(8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P - 16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2 + 16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V**2 &
            &- 2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP + 4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP -&
            &4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP - 4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP &
            &- 8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP + 3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2 +&
            &4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP**2 + 16._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV &
            &+ 16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV - 16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV -&
            &2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV + 4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV &
            &- 4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV +3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*eft_cache%EFTOmegaV &
            &+ 8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV**2 - a**2*eft_cache%EFTGamma2V**2*eft_par_cache%h0_mpc**2*(1 + eft_cache%EFTOmegaV) +&
            &4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaPP*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) + 4._dl*eft_cache%EFTc*(4._dl*eft_cache%EFTGamma5V &
            &+ (1._dl + eft_cache%EFTOmegaV)**2) -2._dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*(a*eft_cache%EFTGamma2P*(1._dl - eft_cache%EFTGamma4V + eft_cache%EFTOmegaV)*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) +&
            &eft_cache%EFTGamma2V*(eft_cache%EFTGamma4V*(-1._dl + 2._dl*a*eft_cache%EFTGamma5P + 2._dl*eft_cache%EFTGamma5V + a*eft_cache%EFTOmegaP - eft_cache%EFTOmegaV) +&
            &(1._dl + eft_cache%EFTOmegaV)*(1._dl + a*(eft_cache%EFTGamma4P - 2._dl*eft_cache%EFTGamma5P - eft_cache%EFTOmegaP) + eft_cache%EFTOmegaV) - 2._dl*eft_cache%EFTGamma5V*(1._dl - a*eft_cache%EFTGamma4P &
            &+ a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV))) +8._dl*eft_cache%EFTGamma5V*eft_cache%Hdot + 16._dl*eft_cache%EFTGamma5V**2*eft_cache%Hdot + 4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%Hdot &
            &+ 8._dl*a*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP*eft_cache%Hdot + 16._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV*eft_cache%Hdot +&
            &16._dl*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV*eft_cache%Hdot + 4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV*eft_cache%Hdot + 8._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV**2*eft_cache%Hdot +&
            &4._dl*eft_cache%EFTGamma4V**2*(eft_cache%adotoa**2*(1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl*eft_cache%EFTGamma5V + a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV) - (1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)*eft_cache%Hdot) +&
            &2._dl*eft_cache%EFTGamma4V*(eft_cache%adotoa**2*(-(a**2*eft_cache%EFTOmegaP**2) + a**2*eft_cache%EFTOmegaPP*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) -&
            &4._dl*(1._dl + eft_cache%EFTOmegaV)*(1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) - a*eft_cache%EFTOmegaP*(3._dl + 2._dl*a*eft_cache%EFTGamma5P &
            &+ 2._dl*eft_cache%EFTGamma5V + 3._dl*eft_cache%EFTOmegaV)) +(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)*(4._dl + a*eft_cache%EFTOmegaP + 4._dl*eft_cache%EFTOmegaV)*eft_cache%Hdot))
        !
        !SP: new gradient term
          ! eft_cache%EFT_gradient = 2._dl*( 2._dl*eft_cache%EFTc*F1**2 +2._dl*F1**2*(eft_cache%Hdot-eft_cache%adotoa**2)/a**2*(1._dl+eft_cache%EFTOmegaV) +F1**2*eft_cache%adotoa/a*( F2-eft_cache%adotoa*eft_cache%EFTOmegaP ) &
          ! & +F1**2*( eft_cache%EFTOmegaP/a*eft_cache%Hdot +eft_cache%adotoa**2*eft_cache%EFTOmegaPP ) -2._dl*F2**2*(1._dl+eft_cache%EFTOmegaV) -F1**2*F2dot +2._dl*F2*F1*F1dot )&
          ! &/( F1*( 3._dl*F2**2+F1*F3 ) )
        ! end if
        ! if (eft_cache%EFTOmegaV==0._dl .and. eft_cache%EFTGamma1V==0._dl .and. eft_cache%EFTGamma2V==0._dl .and. eft_cache%EFTGamma3V==0._dl )then   !SP: to avoid numerical issues
        !   eft_cache%EFT_gradient = 2._dl
        !   return
        ! end if
        ! if (3*F2**2 + 2*F1*F3==0._dl .or. F1**2*F2**4*F3**2==0._dl .or.3*F2**2 + F1*F3==0._dl) then
        !   return
        ! end if
        ! !
        !implement the mass conditions
        eft_cache%EFT_mu1 = ((-4*F1*F2**2*F3*(-2*F1*F3*F2dot +F2*(F3*F1dot + F1*F3dot))**2)/(3*F2**2 + F1*F3) +(4*F1*F3*(-2*F1*F3*F2dot +F2*(F3*F1dot + F1*F3dot))**2)/(3 &
            &+ (F1*F3)/F2**2) -(9*F2**6*(1 + 1/(3.*Sqrt(F2**4/(3*F2**2 + 2*F1*F3)**2)))*(6*F1*F2*F3*(3*F2**2 + F1*F3)*eft_cache%adotoa*(F2*F3*F1dot &
            &+F1*(F2*F3dot + 2*F3*(-F2dot + rho_d))) +a*(-3*F2**4*F3**2*F1dot**2 +2*F1*F2**2*F3**2*(-(F3*F1dot**2) + 3*F2**2*F1dotdot) &
            &+2*F1**3*F3*(-2*F2**2*F3dot**2 +F2*F3*(F2*F3dotdot - 2*F3dot*(-2*F2dot + rho_d)) -2*F3**2*(F2dot**2 - F2dot*rho_d + rho_d**2 &
            &+F2*(F2dotdot - rho_d_dot))) +F1**2*F2*(4*F3**3*F1dot*F2dot +F2**3*(-9*F3dot**2 + 6*F3*F3dotdot) +2*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot &
            &+6*(F2dot - rho_d)*rho_d) +12*F2**2*F3*(F3dot*(F2dot - rho_d) +F3*(-F2dotdot + rho_d_dot))))))/(a*(3*F2**2 + F1*F3)**2*(3*F2**2 + 2*F1*F3)) &
            &+(2*F1*F2**2*(6*F1*F2**2*F3*(3*F2**2 + F1*F3)*eft_cache%adotoa*(-2*F3**2*F1dot + 3*F2**2*F3dot +6*F2*F3*(-F2dot + rho_d)) &
            &+a*(6*F2**4*F3**3*F1dot**2 +2*F1**3*F3**2*(-2*F3*F2dot + F2*F3dot)**2 +F1**2*F2*F3*(-8*F3**3*F1dot*F2dot +F2**3*(-9*F3dot**2 + 6*F3*F3dotdot) &
            &+4*F2*F3**2*(F1dot*F3dot - F3*F1dotdot +3*(F2dot - rho_d)*rho_d) +12*F2**2*F3*(F3dot*(F2dot - rho_d) +F3*(-F2dotdot + rho_d_dot))) &
            &+F1*(4*F2**2*F3**4*F1dot**2 +9*F2**6*(-3*F3dot**2 + 2*F3*F3dotdot) -12*F2**4*F3**2*(F3*F1dotdot +3*rho_d*(-F2dot + rho_d)) &
            &+36*F2**5*F3*(F3dot*(F2dot - rho_d) +F3*(-F2dotdot + rho_d_dot))))))/(a*(3*F2**2 + F1*F3)*(3*F2**2 + 2*F1*F3)))/(16.*F1**2*F2**4*F3**2)

        eft_cache%EFT_mu2 = ((-4*F1*F2**2*F3*(-2*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2)/(3*F2**2 + F1*F3) + (4*F1*F3*(-2*F1*F3*F2dot + F2*(F3*F1dot &
            &+ F1*F3dot))**2)/ (3 + (F1*F3)/F2**2) - (2*F2**2*F3*(6*F1*F2**2*F3*(3*F2**2 + F1*F3)*eft_cache%adotoa*(9*F2**4*F1dot + 6*F1*F2**2*F3*F1dot &
            &+ F1**2*(-3*F2**2*F3dot + 2*F3*(F3*F1dot + 3*F2*(F2dot - rho_d)))) + a*(-27*F2**8*F3*F1dot**2 - 2*F1**5*F3**2*(-2*F3*F2dot + F2*F3dot)**2 &
            &+ 18*F1*F2**6*F3*(-2*F3*F1dot**2 + 3*F2**2*F1dotdot) + 18*F1**2*F2**4*F3*(-(F1dot*(F3**2*F1dot - 2*F2*F3*F2dot + F2**2*F3dot)) &
            &+ 3*F2**2*F3*F1dotdot) + 2*F1**3*F2**2*(-2*F3**4*F1dot**2 + 12*F2*F3**3*F1dot*F2dot + 9*F2**4*(F3dot**2 - F3*F3dotdot) + 6*F2**2*F3**2*(&
            &-3*F2dot**2 - F1dot*F3dot + 2*F3*F1dotdot - 3*F2dot*rho_d + 3*rho_d**2) + 18*F2**3*F3*(F3dot*rho_d + F3*(F2dotdot - rho_d_dot))) &
            &+ F1**4*F2*F3*(8*F3**3*F1dot*F2dot + 3*F2**3*(F3dot**2 - 2*F3*F3dotdot) + 4*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot + 3*(-2*F2dot + rho_d)*(F2dot &
            &+ rho_d)) + 12*F2**2*F3*(F3dot*(F2dot + rho_d) + F3*(F2dotdot - rho_d_dot))))))/(a*(3*F2**2 + F1*F3)**2*(3*F2**2 + 2*F1*F3)) &
            &+ (9*F2**6*(1 + 1/(3.*Sqrt(F2**4/(3*F2**2 + 2*F1*F3)**2)))*(6*F1*F2*F3*(3*F2**2 + F1*F3)*eft_cache%adotoa*(F2*F3*F1dot + F1*(F2*F3dot &
            &+ 2*F3*(-F2dot + rho_d))) + a*(-3*F2**4*F3**2*F1dot**2 + 2*F1*F2**2*F3**2*(-(F3*F1dot**2) + 3*F2**2*F1dotdot) + 2*F1**3*F3*(-2*F2**2*F3dot**2 &
            &+ F2*F3*(F2*F3dotdot - 2*F3dot*(-2*F2dot + rho_d)) - 2*F3**2*(F2dot**2 - F2dot*rho_d + rho_d**2 + F2*(F2dotdot - rho_d_dot))) &
            &+ F1**2*F2*(4*F3**3*F1dot*F2dot + F2**3*(-9*F3dot**2 + 6*F3*F3dotdot) + 2*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot + 6*(F2dot - rho_d)*rho_d) &
            &+ 12*F2**2*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/(a*(3*F2**2 + F1*F3)**2*(3*F2**2 + 2*F1*F3)))/(16.*F1**2*F2**4*F3**2)

! if (IsNaN(eft_cache%EFT_mu1) .or. IsNaN(eft_cache%EFT_mu2)) then
!   print*, eft_cache%EFT_mu1, eft_cache%EFT_mu2
!   call self%parameter_values( 1, value1 )
!   call self%parameter_values( 2, value2 )
!   print*, value1, value2
!   print*, '1', 3*F2**2 + 2*F1*F3
!   print*, '2', F1**2*F2**4*F3**2
!   print*, '3', 3*F2**2 + F1*F3
!   print*, F1,F2,F3
! end if

    end subroutine EFTCAMBModelComputeStabilityFactors

    ! ---------------------------------------------------------------------------------------------
    !> Function that computes model specific stability requirements.
    function EFTCAMBModelAdditionalModelStability( self, a, eft_par_cache, eft_cache )

        implicit none

        class(EFTCAMB_model)                         :: self          !< the base class
        real(dl), intent(in)                         :: a             !< the input scale factor.
        type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

        logical :: EFTCAMBModelAdditionalModelStability               !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

        EFTCAMBModelAdditionalModelStability = .True.

    end function EFTCAMBModelAdditionalModelStability

    ! ---------------------------------------------------------------------------------------------

end module EFTCAMB_abstract_model

!----------------------------------------------------------------------------------------
