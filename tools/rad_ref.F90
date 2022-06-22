!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module rad_ref_mod

    use constants_mod, only: grav, rdgas, pi => pi_8
    use fv_arrays_mod, only: fv_grid_bounds_type, r_grid
    use gfdl_mp_mod, only: do_hail, rhor, rhos, rhog, rhoh, rnzr, rnzs, rnzg, rnzh
    use gfdl_mp_mod, only: do_hail_inline => do_hail ! assuming same densities and numbers in both inline and traditional gfdl mp
    use module_mp_radar

   real :: missing_value = - 1.e10

   logical :: module_is_initialized = .false.
   logical :: qsmith_tables_initialized = .false.

   character (len = 17) :: mod_name = 'gfdl_cloud_microphys'

   real, parameter :: n0r = 8.0e6, n0s = 3.0e6, n0g = 4.0e6
   real,                 parameter :: rvgas = 461.50       !< gfs: gas constant for water vapor
   real,                 parameter :: cp_air = 1004.6      !< gfs: heat capacity of dry air at constant pressure
   real,                 parameter :: hlv = 2.5e6          !< gfs: latent heat of evaporation
   real,                 parameter :: hlf = 3.3358e5       !< gfs: latent heat of fusion

   ! real, parameter :: rdgas = 287.04                     !< gfdl: gas constant for dry air

   ! real, parameter :: cp_air = rdgas * 7. / 2.           ! 1004.675, heat capacity of dry air at constant pressure
   real, parameter :: cp_vap = 4.0 * rvgas                 !< 1846.0, heat capacity of water vapore at constnat pressure
   ! real, parameter :: cv_air = 717.56                    ! satoh value
   real, parameter :: cv_air = cp_air - rdgas              !< 717.55, heat capacity of dry air at constant volume
   ! real, parameter :: cv_vap = 1410.0                    ! emanuel value
   real, parameter :: cv_vap = 3.0 * rvgas                 !< 1384.5, heat capacity of water vapor at constant volume

   ! the following two are from emanuel's book "atmospheric convection"
   ! real, parameter :: c_ice = 2106.0                     ! heat capacity of ice at 0 deg c: c = c_ice + 7.3 * (t - tice)
   ! real, parameter :: c_liq = 4190.0                     ! heat capacity of water at 0 deg c

   real, parameter :: c_ice = 1972.0                       !< gfdl: heat capacity of ice at - 15 deg c
   real, parameter :: c_liq = 4185.5                       !< gfdl: heat capacity of water at 15 deg c
   ! real, parameter :: c_liq = 4218.0                     ! ifs: heat capacity of liquid at 0 deg c

   real, parameter :: eps = rdgas / rvgas                  !< 0.6219934995
   real, parameter :: zvir = rvgas / rdgas - 1.            !< 0.6077338443

   real, parameter :: t_ice = 273.16                       !< freezing temperature
   real, parameter :: table_ice = 273.16                   !< freezing point for qs table

   ! real, parameter :: e00 = 610.71                       ! gfdl: saturation vapor pressure at 0 deg c
   real, parameter :: e00 = 611.21                         !< ifs: saturation vapor pressure at 0 deg c

   real, parameter :: dc_vap = cp_vap - c_liq              !< - 2339.5, isobaric heating / cooling
   real, parameter :: dc_ice = c_liq - c_ice               !< 2213.5, isobaric heating / colling

   real, parameter :: hlv0 = hlv                           !< gfs: evaporation latent heat coefficient at 0 deg c
   ! real, parameter :: hlv0 = 2.501e6                     ! emanuel appendix - 2
   real, parameter :: hlf0 = hlf                           !< gfs: fussion latent heat coefficient at 0 deg c
   ! real, parameter :: hlf0 = 3.337e5                     ! emanuel

   real, parameter :: lv0 = hlv0 - dc_vap * t_ice          !< 3.13905782e6, evaporation latent heat coefficient at 0 deg k
   real, parameter :: li00 = hlf0 - dc_ice * t_ice         !< - 2.7105966e5, fusion latent heat coefficient at 0 deg k

   real, parameter :: d2ice = dc_vap + dc_ice              !< - 126, isobaric heating/cooling
   real, parameter :: li2 = lv0 + li00                     !< 2.86799816e6, sublimation latent heat coefficient at 0 deg k

   real, parameter :: qrmin = 1.e-8                        !< min value for rain
   real, parameter :: qvmin = 1.e-20                       !< min value for water vapor (treated as zero)
   real, parameter :: qcmin = 1.e-12                       !< min value for cloud condensates

   real, parameter :: vr_min = 1.e-3                       !< min fall speed for rain
   real, parameter :: vf_min = 1.e-5                       !< min fall speed for cloud ice, snow, graupel

   real, parameter :: dz_min = 1.e-2                       !< use for correcting flipped height

   real, parameter :: sfcrho = 1.2                         !< surface air density

    ! density parameters


   public rhor, rhos, rhog, rhoh, rnzr, rnzs, rnzg, rnzh
   real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw !< constants for accretions
   real :: acco (3, 4)                                     !< constants for accretions
   real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (5), cgmlt (5)

   real :: es0, ces0
   real :: pie, rgrav, fac_rc
   real :: c_air, c_vap

   real :: lati, latv, lats, lat2, lcp, icp, tcp           !< used in Bigg mechanism and wet bulk

   real :: d0_vap                                          !< the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
   real :: lv00                                            !< the same as lv0, except that cp_vap can be cp_vap or cv_vap
 
    ! cloud microphysics switchers

   integer :: icloud_f = 0                                 !< cloud scheme
   integer :: irain_f = 0                                  !< cloud water to rain auto conversion scheme

   logical :: de_ice = .false.                             !< to prevent excessive build - up of cloud ice from external sources
   logical :: sedi_transport = .true.                      !< transport of momentum in sedimentation
   logical :: do_sedi_w = .false.                          !< transport of vertical motion in sedimentation
   logical :: do_sedi_heat = .true.                        !< transport of heat in sedimentation
   logical :: prog_ccn = .false.                           !< do prognostic ccn (yi ming's method)
   logical :: do_qa = .true.                               !< do inline cloud fraction
   logical :: rad_snow = .true.                            !< consider snow in cloud fraciton calculation
   logical :: rad_graupel = .true.                         !< consider graupel in cloud fraction calculation
   logical :: rad_rain = .true.                            !< consider rain in cloud fraction calculation
   logical :: fix_negative = .false.                       !< fix negative water species
   logical :: do_setup = .true. !< setup constants and parameters
   logical :: p_nonhydro = .false.                         !< perform hydrosatic adjustment on air density

   real, allocatable :: table (:), table2 (:), table3 (:), tablew (:)
   real, allocatable :: des (:), des2 (:), des3 (:), desw (:)

    logical :: tables_are_initialized = .false.

   ! logical :: master
   ! integer :: id_rh, id_vtr, id_vts, id_vtg, id_vti, id_rain, id_snow, id_graupel, &
   ! id_ice, id_prec, id_cond, id_var, id_droplets
   real, parameter :: dt_fr = 8.                           !< homogeneous freezing of all cloud water at t_wfr - dt_fr
   ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
   ! dt_fr can be considered as the error bar

   real :: p_min = 100.                                    !< minimum pressure (pascal) for mp to operate

   ! slj, the following parameters are for cloud - resolving resolution: 1 - 5 km

   ! qi0_crt = 0.8e-4
   ! qs0_crt = 0.6e-3
   ! c_psaci = 0.1
   ! c_pgacs = 0.1

   ! -----------------------------------------------------------------------
   ! namelist parameters
   ! -----------------------------------------------------------------------

   real :: cld_min = 0.05                                  !< minimum cloud fraction
   real :: tice = 273.16                                   !< set tice = 165. to trun off ice - phase phys (kessler emulator)

   real :: t_min = 178.                                    !< min temp to freeze - dry all water vapor
   real :: t_sub = 184.                                    !< min temp for sublimation of cloud ice
   real :: mp_time = 150.                                  !< maximum micro - physics time step (sec)

   ! relative humidity increment

   real :: rh_inc = 0.25                                   !< rh increment for complete evaporation of cloud water and cloud ice
   real :: rh_inr = 0.25                                   !< rh increment for minimum evaporation of rain
   real :: rh_ins = 0.25                                   !< rh increment for sublimation of snow

   ! conversion time scale

   real :: tau_r2g = 900.                                  !< rain freezing during fast_sat
   real :: tau_smlt = 900.                                 !< snow melting
   real :: tau_g2r = 600.                                  !< graupel melting to rain
   real :: tau_imlt = 600.                                 !< cloud ice melting
   real :: tau_i2s = 1000.                                 !< cloud ice to snow auto-conversion
   real :: tau_l2r = 900.                                  !< cloud water to rain auto-conversion
   real :: tau_v2l = 150.                                  !< water vapor to cloud water (condensation)
   real :: tau_l2v = 300.                                  !< cloud water to water vapor (evaporation)
   real :: tau_g2v = 900.                                  !< graupel sublimation
   real :: tau_v2g = 21600.                                !< graupel deposition -- make it a slow process

   ! horizontal subgrid variability

   real :: dw_land = 0.20                                  !< base value for subgrid deviation / variability over land
   real :: dw_ocean = 0.10                                 !< base value for ocean

   ! prescribed ccn

   real :: ccn_o = 90.                                     !< ccn over ocean (cm^ - 3)
   real :: ccn_l = 270.                                    !< ccn over land (cm^ - 3)

   real :: rthresh = 10.0e-6                               !< critical cloud drop radius (micro m)
    ! -----------------------------------------------------------------------
   ! wrf / wsm6 scheme: qi_gen = 4.92e-11 * (1.e3 * exp (0.1 * tmp)) ** 1.33
   ! optimized: qi_gen = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
   ! qi_gen ~ 4.808e-7 at 0 c; 1.818e-6 at - 10 c, 9.82679e-5 at - 40c
   ! the following value is constructed such that qc_crt = 0 at zero c and @ - 10c matches
   ! wrf / wsm6 ice initiation scheme; qi_crt = qi_gen * min (qi_lim, 0.1 * tmp) / den
   ! -----------------------------------------------------------------------

   real :: sat_adj0 = 0.90                                 !< adjustment factor (0: no, 1: full) during fast_sat_adj

   real :: qc_crt = 5.0e-8                                 !< mini condensate mixing ratio to allow partial cloudiness

   real :: qi_lim = 1.                                     !< cloud ice limiter to prevent large ice build up

   real :: ql_mlt = 2.0e-3                                 !< max value of cloud water allowed from melted cloud ice
   real :: qs_mlt = 1.0e-6                                 !< max cloud water due to snow melt

   real :: ql_gen = 1.0e-3                                 !< max cloud water generation during remapping step if fast_sat_adj = .t.
   real :: qi_gen = 1.82e-6                                !< max cloud ice generation during remapping step

   ! cloud condensate upper bounds: "safety valves" for ql & qi

   real :: ql0_max = 2.0e-3                                !< max cloud water value (auto converted to rain)
   real :: qi0_max = 1.0e-4                                !< max cloud ice value (by other sources)

   real :: qi0_crt = 1.0e-4                                !< cloud ice to snow autoconversion threshold (was 1.e-4);
                                                           !! qi0_crt is highly dependent on horizontal resolution
   real :: qr0_crt = 1.0e-4                                !< rain to snow or graupel/hail threshold
                                                           ! lfo used * mixing ratio * = 1.e-4 (hail in lfo)
   real :: qs0_crt = 1.0e-3                                !< snow to graupel density threshold (0.6e-3 in purdue lin scheme)

   real :: c_paut = 0.55                                   !< autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
   real :: c_psaci = 0.02                                  !< accretion: cloud ice to snow (was 0.1 in zetac)
   real :: c_piacr = 5.0                                   !< accretion: rain to ice:
   real :: c_cracw = 0.9                                   !< rain accretion efficiency
   real :: c_pgacs = 2.0e-3                                !< snow to graupel "accretion" eff. (was 0.1 in zetac)

   ! decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)

   real :: alin = 842.0                                    !< "a" in lin1983
   real :: clin = 4.8                                      !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)

   ! fall velocity tuning constants:

   logical :: const_vi = .false.                           !< if .t. the constants are specified by v * _fac
   logical :: const_vs = .false.                           !< if .t. the constants are specified by v * _fac
   logical :: const_vg = .false.                           !< if .t. the constants are specified by v * _fac
   logical :: const_vr = .false.                           !< if .t. the constants are specified by v * _fac

   ! good values:

   real :: vi_fac = 1.                                     !< if const_vi: 1 / 3
   real :: vs_fac = 1.                                     !< if const_vs: 1.
   real :: vg_fac = 1.                                     !< if const_vg: 2.
   real :: vr_fac = 1.                                     !< if const_vr: 4.

   ! upper bounds of fall speed (with variable speed option)

   real :: vi_max = 0.5                                    !< max fall speed for ice
   real :: vs_max = 5.0                                    !< max fall speed for snow
   real :: vg_max = 8.0                                    !< max fall speed for graupel
   real :: vr_max = 12.                                    !< max fall speed for rain

   ! cloud microphysics switchers

   logical :: fast_sat_adj = .false.                       !< has fast saturation adjustments
   logical :: z_slope_liq = .true.                         !< use linear mono slope for autocconversions
   logical :: z_slope_ice = .false.                        !< use linear mono slope for autocconversions
   logical :: use_ccn = .false.                            !< must be true when prog_ccn is false
   logical :: use_ppm = .false.                            !< use ppm fall scheme
   logical :: mono_prof = .true.                           !< perform terminal fall with mono ppm scheme
   logical :: mp_print = .false.                           !< cloud microphysics debugging printout

   ! real :: global_area = - 1.

   real :: log_10, tice0, t_wfr

    integer :: reiflag = 1
    ! 1: Heymsfield and Mcfarquhar, 1996
    ! 2: Wyser, 1998

    logical :: tintqs = .false. !< use temperature in the saturation mixing in PDF

    real :: rewmin = 5.0, rewmax = 10.0
    real :: reimin = 10.0, reimax = 150.0
    real :: rermin = 10.0, rermax = 10000.0
    real :: resmin = 150.0, resmax = 10000.0
    real :: regmin = 300.0, regmax = 10000.0

contains

subroutine rad_ref (q, pt, delp, peln, delz, dbz, maxdbz, allmax, bd, &
        npz, ncnst, hydrostatic, zvir, in0r, in0s, in0g, iliqskin, do_inline_mp, &
        sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, mp_top)

    ! code from mark stoelinga's dbzcalc.f from the rip package.
    ! currently just using values taken directly from that code, which is
    ! consistent for the mm5 reisner - 2 microphysics. from that file:

    ! this routine computes equivalent reflectivity factor (in dbz) at
    ! each model grid point. in calculating ze, the rip algorithm makes
    ! assumptions consistent with those made in an early version
    ! (ca. 1996) of the bulk mixed - phase microphysical scheme in the mm5
    ! model (i.e., the scheme known as "resiner - 2") . for each species:
    !
    ! 1. particles are assumed to be spheres of constant density. the
    ! densities of rain drops, snow particles, and graupel particles are
    ! taken to be rho_r = rho_l = 1000 kg m^ - 3, rho_s = 100 kg m^ - 3, and
    ! rho_g = 400 kg m^ - 3, respectively. (l refers to the density of
    ! liquid water.)
    !
    ! 2. the size distribution (in terms of the actual diameter of the
    ! particles, rather than the melted diameter or the equivalent solid
    ! ice sphere diameter) is assumed to follow an exponential
    ! distribution of the form n (d) = n_0 * exp (lambda * d) .
    !
    ! 3. if in0x = 0, the intercept parameter is assumed constant (as in
    ! early reisner - 2), with values of 8x10^6, 2x10^7, and 4x10^6 m^ - 4,
    ! for rain, snow, and graupel, respectively. various choices of
    ! in0x are available (or can be added) . currently, in0x = 1 gives the
    ! variable intercept for each species that is consistent with
    ! thompson, rasmussen, and manning (2004, monthly weather review,
    ! vol. 132, no. 2, pp. 519 - 542.)
    !
    ! 4. if iliqskin = 1, frozen particles that are at a temperature above
    ! freezing are assumed to scatter as a liquid particle.
    !
    ! more information on the derivation of simulated reflectivity in rip
    ! can be found in stoelinga (2005, unpublished write - up) . contact
    ! mark stoelinga (stoeling@atmos.washington.edu) for a copy.

    ! 22sep16: modifying to use the gfdl mp parameters. if doing so remember
    ! that the gfdl mp assumes a constant intercept (in0x = .false.)
    ! ferrier - aligo has an option for fixed slope (rather than fixed intercept) .
    ! thompson presumably is an extension of reisner mp.

    implicit none

    type (fv_grid_bounds_type), intent (in) :: bd

    logical, intent (in) :: hydrostatic, in0r, in0s, in0g, iliqskin, do_inline_mp

    integer, intent (in) :: npz, ncnst, mp_top
    integer, intent (in) :: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel

    real, intent (in), dimension (bd%isd:bd%ied, bd%jsd:bd%jed, npz) :: pt, delp
    real, intent (in), dimension (bd%is:, bd%js:, 1:) :: delz
    real, intent (in), dimension (bd%isd:bd%ied, bd%jsd:bd%jed, npz, ncnst) :: q
    real, intent (in), dimension (bd%is :bd%ie, npz + 1, bd%js:bd%je) :: peln
    real, intent (out), dimension (bd%is :bd%ie, bd%js :bd%je, npz) :: dbz
    real, intent (out), dimension (bd%is :bd%ie, bd%js :bd%je) :: maxdbz

    real, intent (in) :: zvir
    real, intent (out) :: allmax

    ! parameters for constant intercepts (in0[rsg] = .false.)
    ! using gfdl mp values

    real (kind = r_grid), parameter :: vconr = 2503.23638966667
    real (kind = r_grid), parameter :: vcong = 87.2382675
    real (kind = r_grid), parameter :: vcons = 6.6280504
    real (kind = r_grid), parameter :: vconh = vcong
    real (kind = r_grid), parameter :: normr = 25132741228.7183
    real (kind = r_grid), parameter :: normg = 5026548245.74367
    real (kind = r_grid), parameter :: normh = pi * rhoh * rnzh
    real (kind = r_grid), parameter :: norms = 942477796.076938

    ! constants for variable intercepts
    ! will need to be changed based on mp scheme

    real, parameter :: r1 = 1.e-15
    real, parameter :: ron = 8.e6
    real, parameter :: ron2 = 1.e10
    real, parameter :: son = 2.e7
    real, parameter :: gon = 5.e7
    real, parameter :: ron_min = 8.e6
    real, parameter :: ron_qr0 = 0.00010
    real, parameter :: ron_delqr0 = 0.25 * ron_qr0
    real, parameter :: ron_const1r = (ron2 - ron_min) * 0.5
    real, parameter :: ron_const2r = (ron2 + ron_min) * 0.5

    ! other constants

    real, parameter :: gamma_seven = 720.
    real, parameter :: alpha = 0.224
    real (kind = r_grid), parameter :: factor_s = gamma_seven * 1.e18 * (1. / (pi * rhos)) ** 1.75 &
         * (rhos / rhor) ** 2 * alpha
    real, parameter :: qmin = 1.e-12
    real, parameter :: tice = 273.16

    ! double precision

    real (kind = r_grid), dimension (bd%is:bd%ie) :: rhoair, denfac, z_e
    real (kind = r_grid) :: qr1, qs1, qg1, t1, t2, t3, rwat, vtr, vtg, vts
    real (kind = r_grid) :: factorb_s, factorb_g
    real (kind = r_grid) :: temp_c, pres, sonv, gonv, ronv

    real :: rhogh, vcongh, normgh

    integer :: i, j, k
    integer :: is, ie, js, je

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    if (rainwat < 1) return

    dbz (:, :, 1:mp_top) = - 20.
    maxdbz (:, :) = - 20. ! minimum value
    allmax = - 20.

    if ((do_hail .and. .not. do_inline_mp) .or. (do_hail_inline .and. do_inline_mp)) then
        rhogh = rhoh
        vcongh = vconh
        normgh = normh
    else
        rhogh = rhog
        vcongh = vcong
        normgh = normg
    endif

    !$omp parallel do default (shared) private (rhoair, t1, t2, t3, denfac, vtr, vtg, vts, z_e)
    do k = mp_top + 1, npz
        do j = js, je
            if (hydrostatic) then
                do i = is, ie
                    rhoair (i) = delp (i, j, k) / ((peln (i, k + 1, j) - peln (i, k, j)) * &
                        rdgas * pt (i, j, k) * (1. + zvir * q (i, j, k, sphum)))
                    denfac (i) = sqrt (min (10., 1.2 / rhoair (i)))
                    z_e (i) = 0.
                enddo
            else
                do i = is, ie
                    rhoair (i) = - delp (i, j, k) / (grav * delz (i, j, k)) ! moist air density
                    denfac (i) = sqrt (min (10., 1.2 / rhoair (i)))
                    z_e (i) = 0.
                enddo
            endif
            if (rainwat > 0) then
                do i = is, ie
                    ! the following form vectorizes better & more consistent with gfdl_mp
                    ! sjl notes: marshall - palmer, dbz = 200 * precip ** 1.6, precip = 3.6e6 * t1 / rhor * vtr ! [mm / hr]
                    ! gfdl_mp terminal fall speeds are used
                    ! date modified 20170701
                    ! account for excessively high cloud water - > autoconvert (diag only) excess cloud water
                    t1 = rhoair (i) * max (qmin, q (i, j, k, rainwat) + dim (q (i, j, k, liq_wat), 1.0e-3))
                    vtr = max (1.e-3, vconr * denfac (i) * exp (0.2 * log (t1 / normr)))
                    z_e (i) = 200. * exp (1.6 * log (3.6e6 * t1 / rhor * vtr))
                    ! z_e = 200. * (exp (1.6 * log (3.6e6 * t1 / rhor * vtr)) + &
                    ! exp (1.6 * log (3.6e6 * t3 / rhogh * vtg)) + &
                    ! exp (1.6 * log (3.6e6 * t2 / rhos * vts)))
                enddo
            endif
            if (graupel > 0) then
                do i = is, ie
                    t3 = rhoair (i) * max (qmin, q (i, j, k, graupel))
                    vtg = max (1.e-3, vcongh * denfac (i) * exp (0.125 * log (t3 / normgh)))
                    z_e (i) = z_e (i) + 200. * exp (1.6 * log (3.6e6 * t3 / rhogh * vtg))
                enddo
            endif
            if (snowwat > 0) then
                do i = is, ie
                    t2 = rhoair (i) * max (qmin, q (i, j, k, snowwat))
                    ! vts = max (1.e-3, vcons * denfac * exp (0.0625 * log (t2 / norms)))
                    z_e (i) = z_e (i) + (factor_s / alpha) * t2 * exp (0.75 * log (t2 / rnzs))
                    ! z_e = 200. * (exp (1.6 * log (3.6e6 * t1 / rhor * vtr)) + &
                    ! exp (1.6 * log (3.6e6 * t3 / rhogh * vtg)) + &
                    ! exp (1.6 * log (3.6e6 * t2 / rhos * vts)))
                enddo
            endif
            do i = is, ie
                dbz (i, j, k) = 10. * log10 (max (0.01, z_e (i)))
            enddo
        enddo
    enddo

    !$omp parallel do default (shared)
    do j = js, je
        do k = mp_top + 1, npz
            do i = is, ie
                maxdbz (i, j) = max (dbz (i, j, k), maxdbz (i, j))
            enddo
        enddo
    enddo

    do j = js, je
        do i = is, ie
            allmax = max (maxdbz (i, j), allmax)
        enddo
    enddo

end subroutine rad_ref

!+---+-----------------------------------------------------------------+
!>\ingroup mod_gfdl_cloud_mp
!! This subroutine calculates radar reflectivity.
      subroutine refl10cm_gfdl (qv1d, qr1d, qs1d, qg1d,                 &
                       t1d, p1d, dBZ, kts, kte, ii, jj, melti)

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii,jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                      qv1d, qr1d, qs1d, qg1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ

!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL, DIMENSION(kts:kte):: rr, rs, rg
!      REAL:: temp_C

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilams, ilamg
      DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_s, N0_g
      DOUBLE PRECISION:: lamr, lams, lamg
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel
      DOUBLE PRECISION:: fmelt_s, fmelt_g

      INTEGER:: i, k, k_0, kbot, n
      LOGICAL, INTENT(IN):: melti
      DOUBLE PRECISION:: cback, x, eta, f_d
!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
!         temp_C = min(-0.001, temp(K)-273.15)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(rdgas*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            N0_r(k) = n0r
            lamr = (xam_r*xcrg(3)*N0_r(k)/rr(k))**(1./xcre(1))
            ilamr(k) = 1./lamr
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            L_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            N0_s(k) = n0s
            lams = (xam_s*xcsg(3)*N0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            L_qs(k) = .false.
         endif

         if (qg1d(k) .gt. 1.E-9) then
            rg(k) = qg1d(k)*rho(k)
            N0_g(k) = n0g
            lamg = (xam_g*xcgg(3)*N0_g(k)/rg(k))**(1./xcge(1))
            ilamg(k) = 1./lamg
            L_qg(k) = .true.
         else
            rg(k) = 1.E-12
            L_qg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      k_0 = kts
      K_LOOP:do k = kte-1, kts, -1
         if ( melti .and. (temp(k).gt.273.15) .and. L_qr(k)             &
              .and. (L_qs(k+1).or.L_qg(k+1)) ) then
            k_0 = MAX(k+1, k_0)
            EXIT K_LOOP
         endif
      enddo K_LOOP

 !+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (xam_g/900.0)*(xam_g/900.0)       &
                                    * N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
      enddo


!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.kts+1) then
       do k = k_0-1, kts, -1

!..Reflectivity contributed by melting snow
          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxDs(n)**xbm_s
              call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif


!..Reflectivity contributed by melting graupel

          if (L_qg(k) .and. L_qg(k_0) ) then
           fmelt_g = MAX(0.005d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = xam_g * xxDg(n)**xbm_g
              call rayleigh_soak_wetgraupel (x,DBLE(xocmg),DBLE(xobmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo


      end subroutine refl10cm_gfdl
!+---+-----------------------------------------------------------------+
!! @}
!! @}
end module rad_ref_mod
