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

!> @file 11_EFTCAMB_stability.f90
!! This file contains the stability detection algorithm of EFTCAMB.


!----------------------------------------------------------------------------------------
!> This module contains the stability detection algorithm of EFTCAMB.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module EFTCAMB_stability

    use precision
    use EFT_def
    use EFTCAMB_cache
    use EFTCAMB_abstract_model
    use EFTCAMB_abstract_model_full
    use EFTCAMB_abstract_model_designer
    use EFTCAMB_main

    implicit none

    private

    public EFTCAMB_Stability_Check, EFTTestStability, EFTStability_cleanup, EFTStabilityComputation

    ! storage for some utility values that needs to be stored for the stability check.
    real(dl), save :: PastA1 = 0._dl
    real(dl), save :: PastAT = 0._dl

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that tests the stability of a theory in a time range.
    !! To ensure best time coverage scans with three different strategies.
    subroutine EFTCAMB_Stability_Check( success, input_EFTCAMB, params_cache, astart, aend, k_max )

        implicit none

        logical                      , intent(out)   :: success         !< Output of the subroutine. Tells whether the model is found stable or not.
        class(EFTCAMB)               , intent(in)    :: input_EFTCAMB   !< the EFTCAMB object for which the code is computing stability.
        type(EFTCAMB_parameter_cache), intent(inout) :: params_cache    !< the EFTCAMB parameter cache that contains all the physical parameters.
        real(dl)                     , intent(in)    :: astart          !< Initial scale factor.
        real(dl)                     , intent(in)    :: aend            !< Final scale factor.
        real(dl)                     , intent(inout) :: k_max           !< the input maximum k mode at which stability is computed.

        ! parameters of the stability sampler:
        integer , parameter :: indMax            = 10000   ! Number of points sampled.
        real(dl), parameter :: LogSamplingScale  = -10._dl ! Where to start with the log sampling

        real(dl) :: Atest, y
        integer  :: ind
        type(EFTCAMB_timestep_cache ) :: eft_cache

        ! 0) initial feedback:
        if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) then
            write(*,'(a)') '***************************************************************'
            write(*,'(a)') ' EFTCAMB: checking stability of the theory'
        end if
        if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
            write(*,'(a)')
        end if

        ! debug open cache files:
        if ( DebugEFTCAMB ) then
            call eft_cache%open_cache_files( input_EFTCAMB%outroot )
        end if

        ! 1) stability code:
        success = .true.

        !    - linear sampling:
        if ( success ) then
            call EFTStability_cleanup()
            do ind=1, indMax
                Atest = astart + REAL(ind-1)*(aend-astart)/REAL(indMax-1)
                success = EFTTestStability( Atest, k_max, input_EFTCAMB, params_cache, eft_cache )
                if ( .not. success ) then
                    if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
                        write(*,*)
                        write(*,'(a,E14.4)') '   Instability detected at a =', Atest
                    end if
                    exit
                end if
            end do
        end if

        !    - log sampling close to astart:
        if ( success ) then
            call EFTStability_cleanup()
            do ind=1, indMax
                y = LogSamplingScale + REAL(ind-1)*(0._dl-LogSamplingScale)/REAL(indMax-1)
                Atest = astart +(aend-astart)*10._dl**y
                success = EFTTestStability( Atest, k_max, input_EFTCAMB, params_cache, eft_cache )
                if ( .not. success ) then
                    if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
                        write(*,*)
                        write(*,'(a,E14.4)') '   Instability detected at a =', Atest
                    end if
                    exit
                end if
            end do
        end if

        !    - log sampling close to aend:
        if ( success ) then
            call EFTStability_cleanup()
            do ind=1, indMax
                Atest = aend +(astart-aend)*10._dl**y
                success = EFTTestStability( Atest, k_max, input_EFTCAMB, params_cache, eft_cache )
                if ( .not. success ) then
                    if ( input_EFTCAMB%EFTCAMB_feedback_level > 2 ) then
                        write(*,*)
                        write(*,'(a,E14.4)') '   Instability detected at a =', Atest
                    end if
                    exit
                end if
            end do
        end if

        ! debug close cache files:
        if ( DebugEFTCAMB ) then
            call eft_cache%close_cache_files( )
        end if

        ! 2) final feedback:
        if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) then
            if ( success ) then
                write(*,'(a)') ' EFTCAMB: theory stable'
            else
                write(*,'(a)') ' EFTCAMB: theory unstable'
            end if
        end if

    end subroutine EFTCAMB_Stability_Check

    ! ---------------------------------------------------------------------------------------------
    !> Function that fills the caches to check the stability of the theory.
    function EFTTestStability( a, k_max, input_EFTCAMB, params_cache, eft_cache )

        implicit none

        real(dl)                     , intent(in)    :: a                       !< the input scale factor.
        real(dl)                     , intent(inout) :: k_max                   !< the input maximum k mode at which stability is computed.
        class(EFTCAMB)               , intent(in)    :: input_EFTCAMB           !< the EFTCAMB object for which the code is computing stability.
        type(EFTCAMB_parameter_cache), intent(inout) :: params_cache            !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache               !< the EFTCAMB timestep cache that contains all the physical values.
        logical                                      :: EFTTestStability        !< Logical value returned by the function. If the model is stable this is True, otherwise False.

        ! Definitions of variables:
        logical  :: EFT_HaveNan_parameter, EFT_HaveNan_timestep
        real(dl) :: EFT_instability_rate, tempk, temp1, temp2, temp3, temp4, temp5
        integer  :: ind_max, ind
        real(dl) :: dtauda, test_dtauda
        external :: dtauda


	real(dl) :: EFTc, EFTcdot, adotoa, Hdot, Hdotdot, EFTOmegaV, EFTOmegaP, EFTOmegaPP, EFTOmegaPPP, EFTGamma1V, EFTGamma1P, EFTGamma2V, EFTGamma2P, EFTGamma3V, EFTGamma3P, EFTGamma4V, EFTGamma4P, EFTtemp_H0
	real(dl) :: F1, F2, F3, F1dot, F2dot, F3dot, F1dotdot, F2dotdot, F3dotdot, EFT_mu1, EFT_mu2
        real(dl) :: EFTcdotdot,EFTOmegaPPPP,  w_m, gpinudotdot_tot,EFTGamma4PP, EFTGamma3PP,EFTw_0, EFTw_1, EFTw_2, EFTGamma1PP, EFTGamma2PP, rho_d, rho_d_dot, EFTcP, EFTcPP, Hddd, EFT_GPINU_TOT

        ! Stability check initialization
        EFTTestStability = .true.
        ! reset the time-step cache:
        call eft_cache%initialize()
        ! fill it:
        call EFTStabilityComputation( a, input_EFTCAMB%model, params_cache, eft_cache )
        ! protect against k_max too small:
        if ( k_max < 0.1_dl ) k_max = 0.1_dl

        ! check stability of the theory:

        ! 0) dtauda should be finite:
        test_dtauda = dtauda(a)
        if ( test_dtauda > HUGE(test_dtauda) .or. IsNaN(test_dtauda) ) then
            EFTTestStability = .false.
            if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a)') '   Model dtauda is Nan'
            return
        end if
        ! 1) everything inside the parameter cache should not be a NaN:
        call params_cache%is_nan( EFT_HaveNan_parameter )
        if ( EFT_HaveNan_parameter ) then
            EFTTestStability = .false.
            if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a)') '   Model has Nan in the parameter cache'
            return
        end if
        ! 2) everything inside the time-step cache should not be a NaN:
        call eft_cache%is_nan( EFT_HaveNan_timestep )
        if ( EFT_HaveNan_timestep ) then
            EFTTestStability = .false.
            if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a)') '   Model has Nan in the timestep cache'
            return
        end if
        ! 3) enforce mathematical stability:
        if ( input_EFTCAMB%EFT_mathematical_stability ) then

            ! 1- the A coefficient should not change sign in time and in k, i.e. it shall not be zero.
            !    This is the strongest stability constraint since violating it would violate the mathematical
            !    consistency of the pi field equation.
            !    The first condition is A1/=0. Implemented by detecting sign changes in A1.
            if ( eft_cache%EFTpiA1*PastA1 < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a)') '   Mathematical instability: A is zero in time'
            end if
            PastA1 = eft_cache%EFTpiA1
            !    The second one is the condition on k.
            if ( (eft_cache%EFTpiA1 > 0 .and. eft_cache%EFTpiA1 + k_max**2*eft_cache%EFTpiA2 < 0) .or. &
                &(eft_cache%EFTpiA1 < 0 .and. eft_cache%EFTpiA1 + k_max**2*eft_cache%EFTpiA2 > 0) ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a)') '   Mathematical instability: A is zero in k'
            end if

            ! 2- the AT coefficient should not change sign in time, i.e. it shall not be zero.
            !    This is the second strongest stability constraint since violating it would
            !    violate the mathematical consistency of the tensor perturbation equation.
            !    Implemented by detecting sign changes in AT.
            if ( eft_cache%EFTAT*PastAT < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a)') '   Mathematical instability:  AT is zero in time'
            end if
            PastAT = eft_cache%EFTAT

            ! 3- we do not want (fast) growing exponential modes.
            !    This condition prevents the pi field from growing exponentially and destroying everything.
            !    Even though this condition is neither completely related to physics nor mathematics,
            !    violating it would completely mess up cosmological observables.

            !    This is the maximum allowed rate of instability. Units shall be Mpc^-1.
            EFT_instability_rate = 0._dl

            !    This condition needs to be tested in k. Sample in k.
            ind_max = 10

            do ind = 1, ind_max
                ! kmode to test. Linear sampling. Should suffice... (??)
                tempk = 0._dl + REAL(ind-1)*(k_max)/REAL(ind_max-1)
                ! vaule that discriminates between different cases:
                temp1 = (eft_cache%EFTpiB1 +eft_cache%EFTpiB2*tempk**2)
                temp2 = (eft_cache%EFTpiA1 +eft_cache%EFTpiA2*tempk**2)
                temp3 = temp1**2 -4._dl*temp2*(eft_cache%EFTpiC +eft_cache%EFTpiD1*tempk**2 + eft_cache%EFTpiD2*tempk**4)

                ! case 1:
                if ( temp3 > 0._dl .and. temp2 /= 0._dl ) then
                    temp4 = +0.5_dl*(-temp1 +sqrt(temp3))/temp2
                    temp5 = +0.5_dl*(-temp1 -sqrt(temp3))/temp2
                    if ( temp4>EFT_instability_rate .or. temp5>EFT_instability_rate ) then
                        EFTTestStability = .false.
                        if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) then
                            write(*,'(a,E11.4)')       '   Mathematical instability: growing exponential at k =', tempk
                            write(*,'(a,E11.4,E11.4)') '      Rate of instability: ', temp4, temp5
                        end if
                        exit
                    end if
                ! case 2:
                else if ( temp2 /= 0._dl ) then
                    temp4 = -0.5_dl*temp1/temp2
                    if ( temp4>EFT_instability_rate ) then
                        EFTTestStability = .false.
                        if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) then
                            write(*,'(a,E11.4)')       '   Mathematical instability: growing exponential at k =', tempk
                            write(*,'(a,E11.4,E11.4)') '      Rate of instability: ', temp4
                        end if
                        exit
                    end if
                end if

            end do

        end if


	!444) Mass instability
	if ( input_EFTCAMB%EFT_mass_stability ) then


	EFTc	= eft_cache%EFTc
	EFTcdot = eft_cache%EFTcdot
	adotoa = eft_cache%adotoa
	Hdot = eft_cache%Hdot
	Hdotdot = eft_cache%Hdotdot
	EFTOmegaV = eft_cache%EFTOmegaV
	EFTOmegaP = eft_cache%EFTOmegaP
	EFTOmegaPP = eft_cache%EFTOmegaPP
	EFTOmegaPPP = eft_cache%EFTOmegaPPP
	EFTGamma1V = eft_cache%EFTGamma1V
	EFTGamma1P = eft_cache%EFTGamma1P
	EFTGamma2V = eft_cache%EFTGamma2V
        EFTGamma2P = eft_cache%EFTGamma2P
	EFTGamma3V = eft_cache%EFTGamma3V
        EFTGamma3P = eft_cache%EFTGamma3P
        EFTGamma3PP = 0._dl!eft_cache%EFTGamma3PP
        EFTGamma1PP = 0._dl
        EFTGamma2PP = 0._dl
	EFTGamma4V = eft_cache%EFTGamma4V
        EFTGamma4P = eft_cache%EFTGamma4P
        EFTGamma4PP = eft_cache%EFTGamma4PP
	EFTtemp_H0 = params_cache%h0_Mpc
        EFTOmegaPPPP = 0._dl

	EFTcP = EFTcdot/(a*adotoa) + 2._dl*EFTc/a
        rho_d = (eft_cache%grhob_t+ eft_cache%grhoc_t)/a**2
        rho_d_dot = -3._dl*adotoa/a*rho_d

        gpinudotdot_tot= 0._dl !SP not massive nu implemented yet!
        EFTw_0 = -1._dl
        EFTw_1 = 0._dl
        EFTw_2 = 0._dl
        w_m = 1._dl/3._dl

        Hddd = (( + eft_cache%grhov_t*Hdot +6*EFTw_0*eft_cache%grhov_t*Hdot &
            &+ 9*EFTw_0**2*eft_cache%grhov_t*Hdot +( eft_cache%grhob_t +eft_cache%grhoc_t)*Hdot &
            &+ 6*w_m*( eft_cache%grhob_t +eft_cache%grhoc_t)*Hdot +9*w_m**2*( eft_cache%grhob_t +eft_cache%grhoc_t)*Hdot &
            &+ eft_cache%grhonu_tot*Hdot -3*a*eft_cache%grhov_t*Hdot*EFTw_1) +adotoa*(&
            &-9*eft_cache%gpinudot_tot -3*adotoa*(eft_cache%gpinu_tot + eft_cache%grhonu_tot)) - 3*gpinudotdot_tot &
            &-adotoa**2*(6*eft_cache%gpinu_tot + ( eft_cache%grhob_t +eft_cache%grhoc_t) + 9*w_m*( eft_cache%grhob_t +eft_cache%grhoc_t)&
            & + 27*w_m**2*( eft_cache%grhob_t +eft_cache%grhoc_t) +27*w_m**3*( eft_cache%grhob_t +eft_cache%grhoc_t) - 2*eft_cache%grhonu_tot &
            &+eft_cache%grhov_t*(1 + 27*EFTw_0**2 + 27*EFTw_0**3 - 6*a*EFTw_1) &
            &+EFTw_0*(9 - 27*a*EFTw_1)) + 3*a**2*EFTw_2)/(6._dl)

        EFTcdotdot = (a**2*eft_cache%grhov_t*Hdot*(-3*(1 + EFTw_0)**2 + a*EFTw_1) &
            &+  (EFTOmegaV*(8*Hdot**2 - 2*Hddd) - a*(-(Hdot**2*(EFTOmegaP &
            &- 3*a*EFTOmegaPP)) + adotoa*Hdotdot*(EFTOmegaP + a*EFTOmegaPP) + EFTOmegaP*Hddd)) &
            &+ a*adotoa*adotoa**3* (EFTOmegaP + a*(3*EFTOmegaPP + a*EFTOmegaPPP)) - adotoa* (4*EFTcdot &
            &+ Hdotdot*(-8*EFTOmegaV + a*(EFTOmegaP + 3*a*EFTOmegaPP)) + a*adotoa*Hdot*(-EFTOmegaP &
            &+ a*(5*EFTOmegaPP + 3*a*EFTOmegaPPP))) + adotoa**2*(a**2*eft_cache%grhov_t*(3 - 6*a*EFTw_1 &
            &+ 3*EFTw_0*(5 + EFTw_0*(7 + 3*EFTw_0) - 3*a*EFTw_1) + a**2*EFTw_2) &
            &+  Hdot* (-12*EFTOmegaV + a*(11*EFTOmegaP + 3*a*(EFTOmegaPP - a*EFTOmegaPPP)))) &
            &-  a*adotoa**4* (4*EFTOmegaP + a**2*(3*EFTOmegaPPP + a*EFTOmegaPPPP)))/(2.)

        EFTcPP = ( EFTcdotdot/adotoa +(3._dl-Hdot/adotoa)*EFTcdot +2._dl*adotoa*EFTc )/(a**2*adotoa**2)

        !< general expressions for Fs
        F1 = 2*(1 + EFTOmegaV) + 3*EFTGamma3V + EFTGamma4V

        F2 = EFTtemp_H0**2*EFTGamma2V + (adotoa*(2*(1 + EFTOmegaV) + 3*EFTGamma3V &
            &+ EFTGamma4V))/a + adotoa*EFTOmegaP

        F3 = (2*EFTc)/a**2 + 4*EFTtemp_H0**2*EFTGamma1V - (6*EFTtemp_H0*adotoa*EFTGamma2V)/a &
            &- (3*adotoa**2*(2*(1 + EFTOmegaV) + 3*EFTGamma3V + EFTGamma4V))/a**2 &
            &- (6*adotoa**2*EFTOmegaP)/a

        F1dot = adotoa*(2*EFTOmegaP + 3*EFTGamma3P + EFTGamma4P)

        F2dot = ((2 + 2*EFTOmegaV + 3*EFTGamma3V + EFTGamma4V + a*EFTOmegaP)*Hdot &
            &+  a**2*EFTtemp_H0**2*adotoa*EFTGamma2P + adotoa**2*(-2 - 2*EFTOmegaV &
            &- 3*EFTGamma3V - EFTGamma4V + 2*a*EFTOmegaP + 3*a*EFTGamma3P + a*EFTGamma4P &
            &+ a**2*EFTOmegaPP))/a**2

        F3dot = -((4*EFTc*adotoa + 6*a*EFTtemp_H0*EFTGamma2V*Hdot - 2*adotoa*(a*EFTcP &
            &- 3*(2 + 2*EFTOmegaV + 3*EFTGamma3V + EFTGamma4V + 2*a*EFTOmegaP)*Hdot &
            &+  2*a**3*EFTtemp_H0**2*EFTGamma1P) +  6*a*EFTtemp_H0*adotoa**2*(-EFTGamma2V &
            &+ a*EFTGamma2P) -  3*adotoa**3*(4 + 4*EFTOmegaV + 6*EFTGamma3V + 2*EFTGamma4V &
            &- 3*a*EFTGamma3P - a*EFTGamma4P -  2*a**2*EFTOmegaPP))/a**3)

        F1dotdot = (2*EFTOmegaP*Hdot + Hdot*(3*EFTGamma3P + EFTGamma4P) + a*adotoa**2*(2*EFTOmegaPP &
            &+ 3*EFTGamma3PP + EFTGamma4PP))/a

        F2dotdot = (a**2*EFTtemp_H0**2*Hdot*EFTGamma2P +  adotoa*Hdot*(-8 - 8*EFTOmegaV - 12*EFTGamma3V &
            &- 4*EFTGamma4V + 5*a*EFTOmegaP +  9*a*EFTGamma3P + 3*a*EFTGamma4P + 3*a**2*EFTOmegaPP) &
            &+ (2 + 2*EFTOmegaV + 3*EFTGamma3V + EFTGamma4V + a*EFTOmegaP)*Hdotdot &
            &+ a**3*EFTtemp_H0**2*adotoa**2*EFTGamma2PP +  adotoa**3*(4 + 4*EFTOmegaV &
            &+ 6*EFTGamma3V + 2*EFTGamma4V - 4*a*EFTOmegaP - 6*a*EFTGamma3P - 2*a*EFTGamma4P + 2*a**2*EFTOmegaPP &
            &+ 3*a**2*EFTGamma3PP + a**2*EFTGamma4PP + a**3*EFTOmegaPPP))/a**3
        F3dotdot = (4*EFTc*(3*adotoa**2 - Hdot) +  2*(a*EFTcP*Hdot - 3*(2 + 2*EFTOmegaV &
            &+ 3*EFTGamma3V + EFTGamma4V + 2*a*EFTOmegaP)*Hdot**2 + 2*a**3*EFTtemp_H0**2*Hdot*EFTGamma1P &
            &-  3*a*EFTtemp_H0*EFTGamma2V*Hdotdot) -  6*adotoa*(a*EFTtemp_H0*Hdot*(&
            &-4*EFTGamma2V + 3*a*EFTGamma2P) +  (2 + 2*EFTOmegaV + 3*EFTGamma3V + EFTGamma4V &
            &+ 2*a*EFTOmegaP)*Hdotdot) + adotoa**2*(-8*a*EFTcP +  3*Hdot*(24&
            & + 24*EFTOmegaV + 36*EFTGamma3V + 12*EFTGamma4V + 4*a*EFTOmegaP -   15*a*EFTGamma3P &
            &- 5*a*EFTGamma4P - 10*a**2*EFTOmegaPP) + 2*a**2*(EFTcPP + 2*a**2*EFTtemp_H0**2*EFTGamma1PP)) &
            &- 6*a*EFTtemp_H0*adotoa**3*(2*EFTGamma2V + a*(-2*EFTGamma2P + a*EFTGamma2PP)) &
            &- 3*adotoa**4*(12 + 12*EFTOmegaV + 18*EFTGamma3V + 6*EFTGamma4V - 4*a*EFTOmegaP - 12*a*EFTGamma3P &
            &-  4*a*EFTGamma4P - 2*a**2*EFTOmegaPP + 3*a**2*EFTGamma3PP + a**2*EFTGamma4PP + 2*a**3*EFTOmegaPPP))/a**4

        EFT_mu1 = ((-4*F1*F2**2*F3*(-2*F1*F3*F2dot +F2*(F3*F1dot + F1*F3dot))**2)/(3*F2**2 + F1*F3) +(4*F1*F3*(-2*F1*F3*F2dot +F2*(F3*F1dot + F1*F3dot))**2)/(3 &
            &+ (F1*F3)/F2**2) -(9*F2**6*(1 + 1/(3.*Sqrt(F2**4/(3*F2**2 + 2*F1*F3)**2)))*(6*F1*F2*F3*(3*F2**2 + F1*F3)*adotoa*(F2*F3*F1dot &
            &+F1*(F2*F3dot + 2*F3*(-F2dot + rho_d))) +a*(-3*F2**4*F3**2*F1dot**2 +2*F1*F2**2*F3**2*(-(F3*F1dot**2) + 3*F2**2*F1dotdot) &
            &+2*F1**3*F3*(-2*F2**2*F3dot**2 +F2*F3*(F2*F3dotdot - 2*F3dot*(-2*F2dot + rho_d)) -2*F3**2*(F2dot**2 - F2dot*rho_d + rho_d**2 &
            &+F2*(F2dotdot - rho_d_dot))) +F1**2*F2*(4*F3**3*F1dot*F2dot +F2**3*(-9*F3dot**2 + 6*F3*F3dotdot) +2*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot &
            &+6*(F2dot - rho_d)*rho_d) +12*F2**2*F3*(F3dot*(F2dot - rho_d) +F3*(-F2dotdot + rho_d_dot))))))/(a*(3*F2**2 + F1*F3)**2*(3*F2**2 + 2*F1*F3)) &
            &+(2*F1*F2**2*(6*F1*F2**2*F3*(3*F2**2 + F1*F3)*adotoa*(-2*F3**2*F1dot + 3*F2**2*F3dot +6*F2*F3*(-F2dot + rho_d)) &
            &+a*(6*F2**4*F3**3*F1dot**2 +2*F1**3*F3**2*(-2*F3*F2dot + F2*F3dot)**2 +F1**2*F2*F3*(-8*F3**3*F1dot*F2dot +F2**3*(-9*F3dot**2 + 6*F3*F3dotdot) &
            &+4*F2*F3**2*(F1dot*F3dot - F3*F1dotdot +3*(F2dot - rho_d)*rho_d) +12*F2**2*F3*(F3dot*(F2dot - rho_d) +F3*(-F2dotdot + rho_d_dot))) &
            &+F1*(4*F2**2*F3**4*F1dot**2 +9*F2**6*(-3*F3dot**2 + 2*F3*F3dotdot) -12*F2**4*F3**2*(F3*F1dotdot +3*rho_d*(-F2dot + rho_d)) &
            &+36*F2**5*F3*(F3dot*(F2dot - rho_d) +F3*(-F2dotdot + rho_d_dot))))))/(a*(3*F2**2 + F1*F3)*(3*F2**2 + 2*F1*F3)))/(16.*F1**2*F2**4*F3**2)

        EFT_mu2 = ((-4*F1*F2**2*F3*(-2*F1*F3*F2dot + F2*(F3*F1dot + F1*F3dot))**2)/(3*F2**2 + F1*F3) + (4*F1*F3*(-2*F1*F3*F2dot + F2*(F3*F1dot &
            &+ F1*F3dot))**2)/ (3 + (F1*F3)/F2**2) - (2*F2**2*F3*(6*F1*F2**2*F3*(3*F2**2 + F1*F3)*adotoa*(9*F2**4*F1dot + 6*F1*F2**2*F3*F1dot &
            &+ F1**2*(-3*F2**2*F3dot + 2*F3*(F3*F1dot + 3*F2*(F2dot - rho_d)))) + a*(-27*F2**8*F3*F1dot**2 - 2*F1**5*F3**2*(-2*F3*F2dot + F2*F3dot)**2 &
            &+ 18*F1*F2**6*F3*(-2*F3*F1dot**2 + 3*F2**2*F1dotdot) + 18*F1**2*F2**4*F3*(-(F1dot*(F3**2*F1dot - 2*F2*F3*F2dot + F2**2*F3dot)) &
            &+ 3*F2**2*F3*F1dotdot) + 2*F1**3*F2**2*(-2*F3**4*F1dot**2 + 12*F2*F3**3*F1dot*F2dot + 9*F2**4*(F3dot**2 - F3*F3dotdot) + 6*F2**2*F3**2*(&
            &-3*F2dot**2 - F1dot*F3dot + 2*F3*F1dotdot - 3*F2dot*rho_d + 3*rho_d**2) + 18*F2**3*F3*(F3dot*rho_d + F3*(F2dotdot - rho_d_dot))) &
            &+ F1**4*F2*F3*(8*F3**3*F1dot*F2dot + 3*F2**3*(F3dot**2 - 2*F3*F3dotdot) + 4*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot + 3*(-2*F2dot + rho_d)*(F2dot &
            &+ rho_d)) + 12*F2**2*F3*(F3dot*(F2dot + rho_d) + F3*(F2dotdot - rho_d_dot))))))/(a*(3*F2**2 + F1*F3)**2*(3*F2**2 + 2*F1*F3)) &
            &+ (9*F2**6*(1 + 1/(3.*Sqrt(F2**4/(3*F2**2 + 2*F1*F3)**2)))*(6*F1*F2*F3*(3*F2**2 + F1*F3)*adotoa*(F2*F3*F1dot + F1*(F2*F3dot &
            &+ 2*F3*(-F2dot + rho_d))) + a*(-3*F2**4*F3**2*F1dot**2 + 2*F1*F2**2*F3**2*(-(F3*F1dot**2) + 3*F2**2*F1dotdot) + 2*F1**3*F3*(-2*F2**2*F3dot**2 &
            &+ F2*F3*(F2*F3dotdot - 2*F3dot*(-2*F2dot + rho_d)) - 2*F3**2*(F2dot**2 - F2dot*rho_d + rho_d**2 + F2*(F2dotdot - rho_d_dot))) &
            &+ F1**2*F2*(4*F3**3*F1dot*F2dot + F2**3*(-9*F3dot**2 + 6*F3*F3dotdot) + 2*F2*F3**2*(-(F1dot*F3dot) + F3*F1dotdot + 6*(F2dot - rho_d)*rho_d) &
            &+ 12*F2**2*F3*(F3dot*(F2dot - rho_d) + F3*(-F2dotdot + rho_d_dot))))))/(a*(3*F2**2 + F1*F3)**2*(3*F2**2 + 2*F1*F3)))/(16.*F1**2*F2**4*F3**2)

	 if ( EFT_mu1 < -adotoa**2/a**2 ) then
                ! if ( eft_cache%EFT_mu1 < 0._dl .and. a>aRGR) then
                EFTTestStability = .false.
                if (  input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,*) '   Mass instability: mu_1 instability. mass term = ', EFT_mu1     
	 end if

         if ( EFT_mu2 < -adotoa**2/a**2 ) then
                !   if ( eft_cache%EFT_mu2 < 0._dl .and. a>aRGR) then
                EFTTestStability = .false.
                  if (  input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,*) '   Mass instability: mu_1 instability. mass term = ', EFT_mu2
         end if





	end if



        ! 4) enforce model specific priors:
        if ( input_EFTCAMB%EFT_additional_priors ) then
            EFTTestStability = input_EFTCAMB%model%additional_model_stability( a, params_cache, eft_cache )
            if ( .not. EFTTestStability ) then
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a)') '   Model specific stability criteria are not met'
            end if
        end if
        ! 5) enforce physical viability:
        if ( input_EFTCAMB%EFT_physical_stability ) then

            ! the present conditions extend up to Horndeski. Enforce that:
            if ( (eft_cache%EFTGamma6V /= 0._dl) .or.      &
                & ( (eft_cache%EFTGamma3V + eft_cache%EFTGamma4V) /= 0._dl) ) then
                write(*,'(a)') '   EFTCAMB WARNING: stability for model beyond GLPV has not been worked out.'
                write(*,'(a)') '      It will be added in a future release.'
                write(*,'(a)') '      If you want to run this model disable EFT_physical_stability.'
                EFTTestStability = .false.
                return
            end if

            ! 1- Positive gravitational constant:
            if ( 1._dl +eft_cache%EFTOmegaV <= 0 ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a,E11.4)') '   Physical instability: negative gravitational constant = ', 1._dl +eft_cache%EFTOmegaV
            end if

            ! 2- Ghost condition:
            if ( eft_cache%EFT_kinetic < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a,E11.4)') '   Physical instability: ghost instability. Kinetic term = ', eft_cache%EFT_kinetic
            end if

            ! 3- Gradient instability:
            if ( eft_cache%EFT_gradient < 0._dl ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a,E11.4)') '   Physical instability: gradient instability. Gradient term = ', eft_cache%EFT_gradient
            end if

            ! 4- No tensor ghosts:
            if ( eft_cache%EFTAT < 0 ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a,E11.4)') '   Physical instability: tensor ghost instability. Tensor kinetic term = ', eft_cache%EFTAT
            end if

            ! 5- No tensor gradient:
            if ( eft_cache%EFTDT < 0 ) then
                EFTTestStability = .false.
                if ( input_EFTCAMB%EFTCAMB_feedback_level > 1 ) write(*,'(a,E11.4)') '   Physical instability: tensor gradient instability. Tensor gradient term = ', eft_cache%EFTDT
            end if

        end if

    end function EFTTestStability

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that restores the values stored in the module to the default.
    !! Needed for successive calls.
    subroutine EFTStability_cleanup()

        implicit none

        PastA1  = 0._dl
        PastAT  = 0._dl

    end subroutine EFTStability_cleanup

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that fills the caches to check the stability of the theory.
    subroutine EFTStabilityComputation( a, input_model, params_cache, eft_cache )

        implicit none

        real(dl), intent(in)                         :: a                       !< the input scale factor.
        class(EFTCAMB_model)         , intent(in)    :: input_model             !< the EFTCAMB model for which the code is computing the RGR time.
        type(EFTCAMB_parameter_cache), intent(inout) :: params_cache            !< the EFTCAMB parameter cache that contains all the physical parameters.
        type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache               !< the EFTCAMB timestep cache that contains all the physical values.

        ! Definitions of variables:
        real(dl) :: grhonu, gpinu, grhormass_t, grhonudot, gpinudot
        integer  :: nu_i, ind, ind_max

        ! start filling the cache:
        eft_cache%a = a
        ! compute background densities of different species
        eft_cache%grhob_t = params_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
        eft_cache%grhoc_t = params_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
        eft_cache%grhor_t = params_cache%grhornomass/a/a ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
        eft_cache%grhog_t = params_cache%grhog/a/a       ! 8\pi G_N \rho_{\gamma} a^2: radiation background density
        ! Massive neutrinos terms:
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                grhonu      = 0._dl
                gpinu       = 0._dl
                grhormass_t = params_cache%grhormass(nu_i)/a**2
                call params_cache%Nu_background(a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                eft_cache%grhonu_tot = eft_cache%grhonu_tot +grhormass_t*grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  +grhormass_t*gpinu
            end do
        end if
        ! assemble total densities and pressure:
        eft_cache%grhom_t  = eft_cache%grhob_t +eft_cache%grhoc_t +eft_cache%grhor_t +eft_cache%grhog_t +eft_cache%grhonu_tot
        eft_cache%gpresm_t = (+eft_cache%grhor_t +eft_cache%grhog_t)/3._dl +eft_cache%gpinu_tot
        ! compute the other things:
        select type ( model_temp => input_model )
            ! compute the background and the background EFT functions.
            class is ( EFTCAMB_full_model )
            ! background for full models. Here the expansion history is computed from the
            ! EFT functions. Hence compute them first and then compute the expansion history.
            call input_model%compute_background_EFT_functions( a, params_cache , eft_cache )
            call input_model%compute_adotoa( a, params_cache , eft_cache )
            class is ( EFTCAMB_designer_model )
            ! background for designer models. Here the expansion history is parametrized
            ! and does not depend on the EFT functions. Hence compute first the expansion history
            ! and then the EFT functions.
            call input_model%compute_adotoa( a, params_cache , eft_cache )
        end select
        ! compute massive neutrinos stuff:
        ! Massive neutrinos mod:
        if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
                grhonu      = 0._dl
                gpinu       = 0._dl
                grhonudot   = 0._dl
                gpinudot    = 0._dl
                grhormass_t = params_cache%grhormass(nu_i)/a**2
                call params_cache%Nu_background( a*params_cache%nu_masses(nu_i), grhonu, gpinu )
                eft_cache%grhonu_tot = eft_cache%grhonu_tot + grhormass_t*grhonu
                eft_cache%gpinu_tot  = eft_cache%gpinu_tot  + grhormass_t*gpinu
                eft_cache%grhonudot_tot = eft_cache%grhonudot_tot + grhormass_t*( params_cache%Nu_drho(a*params_cache%nu_masses(nu_i) ,eft_cache%adotoa, grhonu)&
                    & -4._dl*eft_cache%adotoa*grhonu)
                eft_cache%gpinudot_tot  = eft_cache%gpinudot_tot  + grhormass_t*( params_cache%Nu_pidot(a*params_cache%nu_masses(nu_i),eft_cache%adotoa, gpinu )&
                    & -4._dl*eft_cache%adotoa*gpinu)
            end do
        end if
        ! compute pressure dot:
        eft_cache%gpresdotm_t = -4._dl*eft_cache%adotoa*( eft_cache%grhog_t +eft_cache%grhor_t )/3._dl +eft_cache%gpinudot_tot
        ! compute remaining quantities related to H:
        call input_model%compute_H_derivs( a, params_cache , eft_cache )
        ! compute backgrond EFT functions if model is designer:
        select type ( model_temp => input_model )
            class is ( EFTCAMB_designer_model )
            call input_model%compute_background_EFT_functions( a, params_cache , eft_cache )
        end select
        ! compute all other background stuff:
        call input_model%compute_rhoQPQ( a, params_cache , eft_cache )
        ! compute second order EFT functions:
        call input_model%compute_secondorder_EFT_functions( a, params_cache , eft_cache )
        ! Compute pi field equations factors:
        call input_model%compute_pi_factors( a, params_cache , eft_cache )
        ! Compute coefficients for the tensor propagation equation:
        call input_model%compute_tensor_factors( a, params_cache , eft_cache )
        ! Compute kinetic and gradient terms:
        call input_model%compute_stability_factors( a, params_cache , eft_cache )

        ! dump cache if in debug mode:
        if ( DebugEFTCAMB ) then
            call eft_cache%dump_cache_files()
        end if

    end subroutine EFTStabilityComputation

    !----------------------------------------------------------------------------------------

end module EFTCAMB_stability

!----------------------------------------------------------------------------------------
