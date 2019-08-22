*     declaration of variables and Parameters of input
*     ------------------------------------------------------------------

      implicit none

      integer i,       ! number of line
     >        landcov, !land cover type
     >        ntl      ! total number of input lines

      real a,        ! attenuation coefficient
     >     a11,      ! work variable
     >     a13,      ! work variable
     >     a22,      ! work variable
     >     a23,      ! work variable
     >     a33,      ! work variable
     >     a34,      ! work variable
     >     a35,      ! work variable
     >     a41,      ! work variable
     >     a42,      ! work variable
     >     a43,      ! work variable
     >     a44,      ! work variable
     >     a55,      ! work variable
     >     au,       ! allocation fraction of total photosynthesis in leaves
     >     as,       ! allocation fraction of total photosynthesis in stems
     >     af,       ! allocation fraction of total photosynthesis in fine roots
     >     ar,       ! allocation fraction of total photosynthesis in dense roots
     >     apar,     ! absorbed photosynthetically active radiation by upper canopy (MJ m-2 s-1)
     >     abu,      ! absortivity of leaves (adimensional)
     >     abss,     ! absortivity of stems (adimensional)
     >     alpha,    ! intrinsic quantum efficiency for c3 and c4 plant
     >     agub,     ! gross photosynthesis rate (mol co2 m-2 s-1)
     >     anub,     ! net photosynthesis rate (mol co2 m-2 s-1)
     >     avmuir,   ! average diffuse optical depth
     >     awcg,     ! available water content at layer g [0..1]
     >     awcd,     ! available water content at layer d [0..1]
     >     b,        ! exponent of the moisture release equation (adimensional)
     >     b1,       ! work variable
     >     b2,       ! work variable
     >     b3,       ! work variable
     >     b4,       ! work variable
     >     b5,       ! work variable
     >     bsw,      ! short wave balance at layers u and s (W m-2)
     >     bsr,      ! solar radiation balance at layers u and s (W m-2)
     >     c1,       ! work variable
     >     c2,       ! work variable
     >     c3,       ! work variable
     >     cimax,    ! maximum values for ci (for model stability)
     >     ciub,     ! intercellular co2 concentration - broadleaf (mol co2/mol air)
     >     csub,     ! leaf boundary layer co2 concentration - broadleaf (mol co2/mol air)
     >     cw,       ! specifc heat of water (J kg-1 K-1)
     >     cp,       ! specific heat of air at za (allowing for h2o vapor) (J kg-1 K-1)
     >     cu,       ! heat capacity of upper canopy leaves per unit leaf area (J m-2 K-1)
     >     cs,       ! heat capacity of upper canopy stems per unit stem area (J m-2 K-1)
     >     cg,	     ! heat capacity of soil (J m-2 k-1)
     >     cd,       ! heat capacity of soil at center of topmost layer (J m-2 K-1)
     >     cm,       ! specific heat of soil minerals (J kg-1 K-1)
     >     co,       ! specific heat of organic matter(J kg-1 K-1)
     >     co2a,     ! co2 flux at level a  (mol mol-1)
     >     carbu,    ! carbon allocation in leaves (kg C m-2)
     >     carbs,    ! carbon allocation in stems (kg C m-2)
     >     carbf,    ! carbon allocation in fine roots (kg C m-2)
     >     carbr,    ! carbon allocation in dense roots (kg C m-2)
     >     coefmub,  ! coefficients for stomatal conductance - broadleaf trees
     >     coefbub,  ! minimum conductance when net photosynthesis is zero broadleaf trees
     >     csia1,    ! atmospheric stability parameter
     >     csi2g,    ! atmospheric stability parameter
     >     d,        ! deep drainage (mm s-1)
     >     d13,      ! work variable
     >     d16,      ! work variable
     >     d23,      ! work variable
     >     d26,      ! work variable
     >     d34,      ! work variable
     >     d35,      ! work variable
     >     d36,      ! work variable
     >     d45,      ! work variable
     >     d46,      ! work variable
     >     d56,      ! work variable
     >     dtairmax, ! maximum delta t allowed
     >     dtvegmax, ! maximum delta t allowed
     >     dqmax,    ! maximum delta q allowed
     >     dtime,    ! time step (seconds)
     >     da,       ! zero-plane displacement height for upper canopy (m)
     >     date,     ! julian day  data was collected
     >     dripu,    ! drip of water from leaves (kg h2o m-2 s-1)
     >     drips,    ! drip of water from stems  (kg h2o m-2 s-1)
     >     densair,  ! air density   (kg m-3)
     >     densw,    ! water density (kg m-3)
     >     densm,    ! density of soil minerals (kg m-3)
     >     df,	     ! fine roots
     >     dr,	     ! dead dense roots
     >     du,       ! tipical dimension (m)
     >     ds,       ! tipical dimension (m)
     >     e12,      ! work variable
     >     ea,       ! work variable
     >     emu,      ! ir emissivity of upper-leaves veg plane
     >     ems,      ! ir emissivity of upper-stems veg plane
     >     epsilonn,  ! porosity (m-3 m-3)
     >     eg,	     ! thickness of soil layer (m)
     >     ed,       ! thickness of soil layer at center of topmost layer (m)
     >     etu,      ! leaf transpiration (kg H2O m-2 s-1)
     >	   eiu,      ! evaporation of intercepted water at leaves
     >     esat12,   ! saturation presion of air at level 12 (Pa)
     >     esatu,    ! saturation presion of air at leaf temperature (Pa)
     >     esats,    ! saturation presion of air at stem temperature (Pa)
     >     esatg,    ! saturation presion of air at soil temperature (Pa)
     >     ew,       ! weighted factor of water
     >     eair,     ! weighted factor of air
     >     em,       ! weighted factor of soil minerals
     >     fapar,    ! fraction of photosynthetically active radiation absorbed by upper canopy
     >     fwetu,    ! fraction of upper canopy leaf area wetted by intercepted liquid
     >     fwets,    ! fraction of upper canopy stem area wetted by intercepted liquid
     >     feu,	     ! leaf evapotranspiration flux (kg H2O m-2 s-1)
     >     fes,	     ! stem evapotranspiration flux (kg H2O m-2 s-1)
     >     feg,      ! soil evaporation flux (kg H2O m-2 s-1)
     >     fet,	     ! total evaporation flux (kg H2O m-2 s-1)
     >     fea,      ! atmospheric evaporation flux (kg H2O m-2 s-1)
     >     fhs,      ! heat flux in the soil
     >     flh,      ! total flux of latent heat (W m-2)
     >     fsh,      ! total flux of sensible heat (W m-2)
     >     fseng,    ! sensible heat flux from layer g to layer 12(W m-2)
     >     fsens,    ! sensible heat flux from stems to layer 12 (W m-2)
     >     fsenu,    ! sensible heat flux from leaves to layer 12 (W m-2)
     >     fiairg,   ! volume fraction of soil surface
     >     fiaird,   ! volume fraction of soil at center of topmost layer
     >     fim,	     ! volume fraction of soil minerals
     >     fio,	     ! volume fraction of organic matter
     >     fdown,    ! downward ir flux below tree level per overall area
     >     fup,      ! upward ir flux below tree level per overall area
     >     firg,     ! ir fluxes absorbed by upper soil(W m-2)
     >     firu,     ! ir fluxes absorbed by upper leaves(W m-2)
     >     firs,     ! ir fluxes absorbed by upper stems(W m-2)
     >     firatm,   ! downward incident ir flux (W m-2)
     >     fira,     ! ir balance calculate at layers u and s (W m-2)
     >     fsa,      ! downward incident solar flux (W m-2)
     >     fsu,      ! solar flux absorved by leaves (W m-2)
     >     fss,      ! solar flux absorved by stems (W m-2)
     >     fsg,      ! solar flux absorved by soil (W m-2)
     >     fsr,      ! solar flux reflected by leaves (W m-2)
     >     fg,       ! flux of water into the ground layer (kg m-2 s-1)
     >     fd ,      ! flux of water into the deep layer (kg H2O m-2 s-1)
     >     ftg,      ! function of soil temperature - ground layer
     >     ftd,      ! function of soil temperature - deep layer
     >     gg,	     ! function of soil moisture
     >     gd,	     ! function of soil moisture at center of topmost layer
     >     gmax,     ! maximum heat flux in the soil  (W m-2)
     >     gmin,     ! minimum heat flux in the soil  (W m-2)
     >     gpp,      ! gross primary productivity (g C m-2 s-1)
     >     grav,     ! gravitational acceleration (m s-2)
     >     gsub,     ! upper canopy stomatal conductance - broadleaf  (mol H2O m-2 s-1)
     >     gsubmin,  ! absolute minimum stomatal conductances - broadleaf trees
     >     gamaub,   ! leaf respiration coefficient - broadleaf trees
     >     gamap,    ! compensation point for gross photosynthesis (mol mol-1)
     >     gbco2u,   ! boundary layer conductance for co2 (mol co2 m-2 s-1)
     >     hvap,     ! latent heat of vaporization of water(J kg-1)
     >     hmin,     ! minimum sensible heat flux  (W m-2)
     >     hu,       ! respiration rate of soil litter (leaves)
     >     hs,       ! respiration rate of soil litter (stems)
     >     hf,       ! respiration rate of soil litter (fine roots)
     >     hr,       ! respiration rate of soil litter (dense roots)
     >     imax,     ! maximum infiltration capacity
     >     je,	     ! 'light limited' rate of photosynthesis (mol m-2 s-1)
     >     jc, 	     ! 'rubisco limited' rate of photosynthesis (mol m-2 s-1)
     >     kappag,   ! heat conductivities of soil (kg s-1 m-3)
     >     kappad,   ! heat conductivities of soil at center of topmost layer	(kg s-1 m-3)
     >     kb,       ! hydraulic conductivity of the bottom layer of soil (mm s-1)
     >     kg,	     ! hydraulic conductivity of ground layer of soil (mm s-1)
     >     kd,	     ! hydraulic conductivity of deep layer of soil (mm s-1)
     >     ks,       ! saturation hydraulic conductivity (mm s-1)
     >     kw,       ! thermal conductivity of water (W m-1 K-1)
     >     kair,     ! thermal conductivity of air (W m-1 K-1)
     >     km,       ! thermal conductivity of soil minerals (W m-1 K-1)
     >     ko,       ! o2/co2 kinetic parameters (mol mol-1)
     >     kc,       ! o2/co2 kinetic parameters (mol mol-1)
     >     ko15,     ! o2/co2 kinetic parameters at 15 degrees C (mol mol-1)
     >     kc15,     ! o2/co2 kinetic parameters at 15 degrees C (mol mol-1)
     >     lambdasap,! sapwood fraction of the total stem and coarse root biomass
     >     lue,	     ! light use efficiency (kg C MJ-1)
     >     lai,      ! canopy single-sided leaf area index (area leaf/area veg)
     >     lu,	     ! leaf litter carbon
     >     ls,	     ! stem litter carbon
     >     lm,       ! mean distance between leaves in the canopy
     >     mdens,    ! molar density of the air (mol m-3)
     >     mcng,     ! coreection of mass of water for layer g (when soil moisture goes below the minimum limit)
     >     mcxg,     ! correction of mass of water for layer g (when soil moisture goes above the maximum limit)
     >     mcnd,     ! coreection of mass of water for layer d (when soil moisture goes below the minimum limit)
     >     mcxd,     ! correction of mass of water for layer d (when soil moisture goes above the maximum limit)
     >     mo,       ! mass of organic matter (kg o.m.)
     >     n,        ! fraction of carbon lost due to growth
     >     npp,      ! net primary productivity (kg C m-2 s-1)
     >     nee,	     ! net ecosystem exchange (kg C m-2 s-1)
     >     o2conc,   ! o2 concentration (mol mol-1)
     >     par,      ! photosynthetically active radiation incident at upper canopy (w m-2)
     >     pi,       ! 3.141593
     >     pre,      ! surface pressure (Pa)
     >     prp,      ! precipitation rate (mm/hour)
     >     prpu,     ! precipitation intercepted by leaves (kg m-2 s-1)
     >     prps,     ! precipitation intercepted by stems (kg m-2 s-1)
     >     prpg,     ! precipitation intercepted by soil  (kg m-2 s-1)
     >     psiha1,   ! profile diabatic correction factors
     >     psima1,   ! profile diabatic correction factors
     >     psih2g,   ! profile diabatic correction factors
     >     psim2g,   ! profile diabatic correction factors
     >     psie,     ! air entry water potential (mm)
     >     psig,     ! matric potential - ground layer (mm)
     >     psid,     ! matric potential - deep layer (mm)
     >     qsat12,   ! saturation specific humidity of air at z12 (kg h2o/kg air)
     >     qsatg,    ! saturation specific humidity of soil surface (kg h2o/kg air)
     >     q12,      ! specific humidity of air at z12 (kg h2o/kg air)
     >     qa,       ! specific humidity of air at za (kg h2o/kg air)
     >     qu,	     ! specific humidity of air near the leaves (kg h2o/kg air)
     >     qs,	     ! specific humidity of air near the stems (kg h2o/kg air)
     >     qg,	     ! specific humidity of air in the soil (kg h2o/kg air)
     >     r,        ! gas constant for dry air (J kg-1 K-1)
     >     rf,	     ! respiration of fine root (mol CO2 m-2 s-1)
     >     rs,       ! respiration of stems     (mol CO2 m-2 s-1)
     >     rr,       ! respiration of roots     (mol CO2 m-2 s-1)
     >     rsoil,    ! respiration of soil (kg C m-2 s-1)
     >     reflu,    ! shortwave radiation reflectivity for leaves (adimensional)
     >     rh12,     ! relative humidity of air at z12 (%)
     >     rub,	     ! leaf respiration (mol CO2 m-2 s-1)
     >     runoff,   ! superficial drainage (kg H2O m-2 s-1)
     >     rwork,    ! work variable
     >     rwork1,   ! work variable
     >     rwork2,   ! work variable
     >     sl,       ! specific leaf area (m2 leaf per kg C)
     >     ssh,      ! air-vegetation sensible heat transfer coefficients for upper canopy stems (m s-1)
     >     ssv,      ! air-vegetation water vapor   transfer coefficients for upper canopy stems (m s-1)
     >     suh,      ! air-vegetation sensible heat transfer coefficients for upper canopy leaves (m s-1)
     >     suv,      ! air-vegetation water vapor   transfer coefficients for upper canopy leaves (m s-1)
     >     sai,      ! current single-sided stem area index
     >     stresstu, ! soil moisture stress factor for the upper canopy
     >     stressfac,! coefficient for soil moisture stress
     >     stef,     ! stefan-boltzmann constant (W m-2 K-4)
     >     sigah,    ! heat transfer coefficient between air-air (m s-1)
     >     sigav,    ! vapor transfer coefficient between air-air (m s-1)
     >     siggh,    ! heat transfer coefficient between soil-air (m s-1)
     >     siggv,    ! vapor transfer coefficient between soil-air (m s-1)
     >     ta,       ! reference level air temperature (K)
     >     ts,       ! upper canopy stem temperature (K)
     >     tu,       ! upper canopy leaf temperature (K)
     >     t12,      ! air temperature at z12 (K)
     >     tg,       ! soil surface temperature (K)
     >     tf,	     ! leaf temperature (K)
     >     td,       ! soil temperature at center of topmost layer (K)
     >     tga,      ! previous temperature of layer g
     >     tda,      ! previous temperature of layer d
     >     tava,     ! average temperature of air
     >     tavg,     ! average temperature of soil
     >     transmu,  ! shortwave radiation transmissivity for leaves (adimensional)
     >     transms,  ! shortwave radiation transmissivity for stems (adimensional)
     >     tempvm,   ! temperature stress factor for the upper canopy
     > 	   time,     ! time data was collected
     >     thetag,   ! fraction of soil moisture
     >     thetad,   ! fraction of soil moisture at center of topmost layer
     >     thetas,   ! saturation water content (m3 m-3)
     >     thetafc,  ! field capacity (m3 m-3)
     >     thetawilt,! wilting point of the soil (m3 m-3)
     >     tdrip,    ! time constant for leaf and stem liquid drip (seconds)
     >     tau,      ! co2/o2 ratio of kinetic parametrs
     >     tau15,    ! co2/o2 specificity ratio at 15 degrees C
     >     tauu,     ! residence time of carbon in leaf (seconds)
     >     taus,     ! residence time of carbon in stems (seconds)
     >     tauf,     ! residence time of carbon in fine roots (seconds)
     >     taur,     ! residence time of carbon in dense roots (seconds)
     >     ufrica1,  ! friction velocity
     >     ufric2g,  ! friction velocity
     >     u12,      ! wind speed at middle (m s-1)
     >	   ua,       ! W-E component of wind (m s-1)
     >     u1,	     ! wind speed at top (m s-1)
     >     u2,	     ! wind speed at bottom (m s-1)
     >     v,        ! parameter for calculating maximum infiltration
     >     vmax,     ! maximum rubisco capacity (mol co2 m-2 s-1)
     >     vmaxub,   ! maximum rubisco capacity of top leaf at 15 C (mol co2 m-2 s-1)
     >     wumax,    ! maximum intercepted water on upper canopy leaf area(kg h2o m-2)
     >     wsmax,    ! maximum intercepted water on upper canopy stem area(kg h2o m-2)
     >     wu,	     ! intercepted liquid h2o on upper canopy leaf area (kg m-2)
     >	   ws,	     ! intercepted liquid h2o on upper canopy stem area (kg m-2)
     >     wg,	     ! intercepted liquid h2o on soil surface (kg m-2)
     >     wd,	     ! intercepted liquid h2o on at center of topmost layer (kg m-2)
     >     wgmax,    ! maximum water stored at layer g (kg h2o m-2)
     >     wdmax,    ! maximum water stored at layer d (kg h2o m-2)
     >     wgmin,    ! minimum water stored at layer g (kg h2o m-2)
     >     wdmin,    ! minimum water stored at layer d (kg h2o m-2)
     >     wgfc,     ! water stored at field capacity at layer g  (kg h20 m-2)
     >     wdfc,     ! water stored at field capacity at layer d (kg h20 m-2)
     >     w,	     ! leaf width (m)
     >     year,     ! Year data was collected
     >     xg,	     ! relative soil saturation - layer g
     >     xd,	     ! relative soil saturation - layer d
     >     za,       ! height above a surface of atmospheric (m)
     >     z1,	     ! top heights of upper canopy (m)
     >     z2,       ! bottom heights of upper canopy (m)
     >     z12,	     ! middle heights of upper canopy	(m)
     >     zoa,	     ! upper canopy roughness length parameter (m)
     >     zog	     ! soil surface roughness length parameter (m)

      open(10,file ='ios/KM67/KM67.input',STATUS='OLD')
      open(11,file = 'ios/KM67/ir_balance.txt')
      open(12,file = 'ios/KM67/radiation_flux.txt')
      open(13,file = 'ios/KM67/aerodynamics.txt')
      open(14,file = 'ios/KM67/physiology.txt')
      open(15,file = 'ios/KM67/intercepted_water.txt')
      open(16,file = 'ios/KM67/mass_and_energy.txt')
      open(17,file = 'ios/KM67/soil_physics.txt')
      open(18,file = 'ios/KM67/carbon_balance.txt')
      open(19,file = 'ios/KM67output.txt')
     
c     time step in seconds
      dtime = 3600.
      
!     ------------------------------------------------------------------
!     Land Cover type: 0 to water, 1  to evergreen  needleleaf  forest,
!     2 to evergreen broadleaf forest, 3 to deciduous needleleaf forest,
!     4 to deciduous broadleaf forest, 5 to mixed forests, 6 to closed
!     shrubland, 7 to open shrublands, 8 to woody savannas,
!     9 to savannas, 10 to grasslands, 11 to  permanent wetlands,
!     12 to croplands, 13 to urban and built up,
!     14 to cropland natural vegetation mosaic, 15 to snow and ice and
!     16 to barren or sparsely vegetated
!     ------------------------------------------------------------------
      landcov   =       2    !land cover type


