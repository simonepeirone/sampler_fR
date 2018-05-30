module Pure_EFT_Sigma_Mu

    !use EFTDef
    use EFTCAMB_cache
    contains
    
    !---------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !< This function compute the value of the function Sigma, depending on the
    !< scale k and the scale factor a
    function EFTSigma( a, k )

      use precision
      
      implicit none

      real(dl)             :: EFTSigma
      real(dl), intent(in) ::   a, k

      real(dl) :: EFT_H0, EFT_cV, EFTOmegaV, EFTOmegaP, EFTgamma2V, EFTgamma2P
      real(dl) :: EFTgamma3V, EFTgamma3P, H, Hdot, Hdotdot,a2,EFTOmegaPP,grhov_t
      real(dl) :: C_3, C_pi

      a2=a**2

      EFT_H0     = eft_par_cache%h0_Mpc
      EFTOmegaV  = CP%eft_cache%EFTOmegaV
      EFTOmegaP  = EFTOmega( a, 1 )
      EFTOmegaPP  = EFTOmega( a, 2 )
      EFTgamma2V = EFTGamma2( a, 0 )
      EFTgamma2P = EFTGamma2( a, 1 )
      EFTgamma3V = EFTGamma3( a, 0 )
      EFTgamma3P = EFTGamma3( a, 1 )

      ! compute expansion history
      call Expansion_History(a, H, Hdot, Hdotdot)

      if ( CP%EFTflag == 1 ) then !SP: designer approach
        grhov_t=grhov*a**2  !DB
      EFT_cV = (H*H - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
          & - 0.5_dl*a2*H*H*EFTOmegaPP&
          & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))
      else if ( CP%EFTflag == 5) then !SP: designer approach
        EFT_cV     = Omega_Lambda_c( a, 0 )
      end if

      ! compute the coefficients C_3 and C_pi
      C_3  = EFT_cV - EFT_H0/2._dl *( a*H*EFTgamma2V +a**2*H*EFTgamma2P ) + ( 2._dl * H**2 -Hdot )*EFTgamma3V +a*H**2*EFTgamma3P
      C_pi = 1.5_dl *H/a*( Hdotdot -2._dl *H**3 )*EFTOmegaP -3._dl*EFT_cV/a**2*( Hdot -H**2 ) &
          &+1.5_dl*EFT_H0*( ( Hdotdot -H*Hdot -H**3)/a*EFTgamma2V +H*(Hdot -H**2)*EFTgamma2P ) +3._dl/a**2*( Hdot -H**2 )**2*EFTgamma3V

      ! espression for Sigma
      EFTSigma = ( (2._dl +EFTgamma3V/(EFTOmegaV +1._dl))*C_3 +1._dl/(EFTOmegaV +1._dl)*( a*H*EFTOmegaP +H*EFTgamma3V +a*H*EFTgamma3P )&
          &*( 1.5_dl*a*H*EFTOmegaP +H*EFTgamma3V +a/2._dl*EFT_H0*EFTgamma2V +a*H*EFTgamma3P ) &
          &+C_pi *( EFTgamma3V/(EFTOmegaV +1._dl) +2._dl )*a**2/k**2)/&
          &( 2._dl/(EFTOmegaV +1._dl)*( a**2*H*EFTOmegaP +a*H*EFTgamma3V +a**2*H*EFTgamma3P )&
          &*( H*EFTOmegaP +EFT_H0*EFTgamma2V )*( EFTOmegaV +1._dl +EFTgamma3V ) -a**2/2._dl*( H*EFTOmegaP +EFT_H0*EFTgamma2V )**2&
          &+2._dl/( EFTOmegaV +1._dl )*( EFTOmegaV +1._dl +EFTgamma3V )**2*( C_3 +C_pi*a**2/k**2 ))

    end function EFTSigma

    !---------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !< This function compute the value of the function Mu, depending on the
    !< scale k and the scale factor a
    function EFTMu( a, k )

      use precision

      implicit none

      real(dl)             :: EFTMu
      real(dl), intent(in) ::   a, k

      real(dl) :: EFT_H0, EFT_cV, EFTOmegaV, EFTOmegaP, EFTgamma2V, EFTgamma2P
      real(dl) :: EFTgamma3V, EFTgamma3P, H, Hdot, Hdotdot, a2, EFTOmegaPP,grhov_t
      real(dl) :: C_3, C_pi

      a2=a**2

      EFT_H0     = (CP%h0/c_EFT)*1000._dl
      EFTOmegaV  = EFTOmega( a, 0 )
      EFTOmegaP  = EFTOmega( a, 1 )
      EFTOmegaPP  = EFTOmega( a, 2 )
      EFTgamma2V = EFTGamma2( a, 0 )
      EFTgamma2P = EFTGamma2( a, 1 )
      EFTgamma3V = EFTGamma3( a, 0 )
      EFTgamma3P = EFTGamma3( a, 1 )

      ! compute expansion history
      call Expansion_History(a, H, Hdot, Hdotdot)

      if ( CP%EFTflag == 1 ) then !SP: designer approach
        grhov_t=grhov*a**2  !DB
      EFT_cV = (H*H - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
          & - 0.5_dl*a2*H*H*EFTOmegaPP&
          & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))
      else if ( CP%EFTflag == 5) then !SP: designer approach
        EFT_cV     = Omega_Lambda_c( a, 0 )
      end if

      ! compute the coefficients C_3 and C_pi
      C_3  = EFT_cV - EFT_H0/2._dl *( a*H*EFTgamma2V +a**2*H*EFTgamma2P ) + ( 2._dl * H**2 -Hdot )*EFTgamma3V +a*H**2*EFTgamma3P
      C_pi = 1.5_dl *H/a*( Hdotdot -2._dl *H**3 )*EFTOmegaP -3._dl*EFT_cV/a**2*( Hdot -H**2 ) &
          &+1.5_dl*EFT_H0*( ( Hdotdot -H*Hdot -H**3)/a*EFTgamma2V +H*(Hdot -H**2)*EFTgamma2P ) +3._dl/a**2*( Hdot -H**2 )**2*EFTgamma3V

      !expression for Mu
      EFTMu = ( 1._dl/(EFTOmegaV +1._dl)*(a*H*EFTOmegaP +H*EFTgamma3V +a*H*EFTgamma3P)**2 +C_3+C_pi*a**2/k**2)/&
          &( 1._dl/(EFTOmegaV +1._dl)*(a**2*H*EFTOmegaP +a*H*EFTgamma3V +a**2*H*EFTgamma3P)&
          &*(H*EFTOmegaP +EFT_H0*EFTgamma2V)*(EFTOmegaV+1._dl+EFTgamma3V)-a**2/4._dl*(H*EFTOmegaP +EFT_H0*EFTgamma2V)**2&
          &+1._dl/(EFTOmegaV +1._dl)*(EFTOmegaV+1._dl+EFTgamma3V)**2*(C_3 +C_pi*a**2/k**2))

    end function EFTMu
    !---------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !< This function compute the value of the sound speed square, depending on the
    !< scale factor a
    function EFTCs_square( a )

      use precision

      implicit none

      real(dl)             :: EFTCs_square
      real(dl), intent(in) ::   a

      real(dl) :: EFT_H0, EFT_cV, EFTOmegaV, EFTOmegaP, EFTgamma2V, EFTgamma2P,EFTgamma1V
      real(dl) :: EFTgamma3V, EFTgamma3P, H, Hdot, Hdotdot, a2, EFTOmegaPP
      real(dl) :: alpha, alphaB,alphaB_prime, alphaT, alphaK, alphaM
      real(dl) :: grhob_t, grhoc_t, grhor_t, grhog_t, grhov_t, grhormass_t
      real(dl) :: EFT_gpinu_tot, EFT_grhonu_tot, EFT_grhonu, EFT_gpinu, EFT_gpinudot_tot
      real(dl) :: gpres_matter, grho_matter, gpres_dot_matter,grho,gpres,adotoa

      integer :: nu_i

      a2=a**2

      !compute H_0 and the EFTfunctions
      EFT_H0     = (CP%h0/c_EFT)*1000._dl
      EFTOmegaV  = EFTOmega( a, 0 )
      EFTOmegaP  = EFTOmega( a, 1 )
      EFTOmegaPP  = EFTOmega( a, 2 )
      EFTgamma1V = EFTGamma1( a, 0 )
      EFTgamma2V = EFTGamma2( a, 0 )
      EFTgamma2P = EFTGamma2( a, 1 )
      EFTgamma3V = EFTGamma3( a, 0 )
      EFTgamma3P = EFTGamma3( a, 1 )

      if ( CP%EFTflag == 1 ) then !SP: designer approach
        grhov_t=grhov*a**2  !DB
      EFT_cV = (H*H - Hdot)*(EFTOmegaV + a*EFTOmegaP*0.5_dl) &
          & - 0.5_dl*a2*H*H*EFTOmegaPP&
          & + 0.5_dl*grhov_t*(1._dl+EFTw(a,0))
      else if ( CP%EFTflag == 5) then !SP: designer approach
        EFT_cV     = Omega_Lambda_c( a, 0 )
      end if

      ! compute total pressure and density
      grhob_t=grhob/a     !DB
      grhoc_t=grhoc/a     !DB
      grhov_t=grhov*a**2  !DB

      ! massive neutrinos:
      EFT_gpinu_tot     = 0._dl
      EFT_grhonu_tot    = 0._dl
      EFT_gpinudot_tot = 0._dl

      if (CP%Num_Nu_Massive /= 0) then
          do nu_i = 1, CP%Nu_mass_eigenstates
              EFT_grhonu    = 0._dl
              EFT_gpinu     = 0._dl
              grhormass_t=grhormass(nu_i)/a**2
              call Nu_background(a*nu_masses(nu_i),EFT_grhonu,EFT_gpinu)

              EFT_grhonu_tot = EFT_grhonu_tot + grhormass_t*EFT_grhonu
              EFT_gpinu_tot  = EFT_gpinu_tot + grhormass_t*EFT_gpinu
              EFT_gpinudot_tot  = EFT_gpinudot_tot + grhormass_t*(Nu_pidot(a*nu_masses(nu_i),adotoa, EFT_gpinu)&
                    & -4._dl*adotoa*EFT_gpinu)
          end do
      end if

      ! other species:
      grhor_t = grhornomass/a**2
      grhog_t = grhog/a**2
      ! collect all pressure terms:
      gpres_matter = EFT_gpinu_tot + (grhog_t+grhor_t)/3._dl
      grho_matter = grhob_t +grhoc_t +grhor_t +grhog_t +EFT_grhonu_tot

      ! compute expansion history
      call Expansion_History(a, H, Hdot, Hdotdot)

      !compute the alpha functions
      alphaT = -EFTgamma3V/(1._dl+EFTOmegaV+EFTgamma3V)

      alphaM = a*( EFTOmegaP +EFTGamma3P )/(1._dl +EFTOmegaV +EFTgamma3V)

      alphaK = ( 2._dl*EFT_cV +4._dl*EFT_H0**2*EFTgamma1V*a2 )/( (1._dl+EFTOmegaV+EFTgamma3V)*H**2 )

      alphaB = 0.5_dl*( a*EFTgamma2V*EFT_H0+a*H*EFTOmegaP )/( H*( 1._dl+EFTOmegaV+EFTgamma3V ) )

      alphaB_prime = 0.5_dl*( ( EFT_H0*(EFTgamma2V +a*EFTgamma2P) +(H +Hdot/H)*EFTOmegaP +a*H*EFTOmegaPP)/( H*(1._dl+EFTOmegaV+EFTgamma3V) ) &
                & -( (a*EFT_H0*EFTgamma2V +a*H*EFTOmegaV)*(Hdot/(a*H)*( 1._dl+EFTOmegaV+EFTgamma3V ) +H*(EFTOmegaP +EFTgamma3P) ) )/( &
                &H**2*( 1_dl+EFTOmegaV+EFTgamma3V )**2 )  )

      alpha = alphaK +3._dl/2._dl*alphaB**2

      !compute the sound speed
      !EFTCs_square = 2._dl/alpha*( ( 1._dl -0.5_dl*alphaB )*( alphaM -alphaT +0.5_dl*alphaB*( 1._dl +alphaT ) +( 1._dl -Hdot/H**2 ) )  &
       !             & +0.5_dl*a*alphaB_prime -(gpres_matter+grho_matter)/(2._dl*H*( 1._dl+EFTOmegaV +EFTgamma3V ) ) )
      EFTCs_square = 2._dl/alpha*( ( 1._dl -0.5_dl*alphaB )*( alphaM -alphaT +0.5_dl*alphaB*( 1._dl +alphaT ) +( 1._dl -Hdot/H**2 ) )  &
                    & +0.5_dl*a*alphaB_prime -(gpres_matter+grho_matter)/(2._dl*H**2*( 1._dl+EFTOmegaV +EFTgamma3V ) ) )

    end function EFTCs_square

  !----------------------------------------------------------------------------


end module Pure_EFT_Sigma_Mu
