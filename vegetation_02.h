*     initial parameters of the vegetation:
*     evergreen broadleaf forest
*     ------------------------------------------------------------------

c     total number of input lines
      ntl = 10464 !(km67_calib)

c     medium atmospheric pressure
      pre = 97585.13 !(km67_calib)

c     shortwave radiation reflectivity for leaves
      reflu = 0.090

c     Value of fraction of photosynthetically active radiation absorbed
c     by upper canopy
      fapar = 0.87   ! Senna et al. (2005)

c     shortwave radiation transmissivity for leaves and stems
c
      transmu = 0.16  !transmu*transms = 0.036,
      transms = 0.19  !according to Moura et al. (2000)

c     use uniform value 1.0 for average diffuse optical depth
      avmuir = 1.

c     initial values for canopy leaves and stems and soil temperature
      tu = 298.16
      ts = 298.16
      tg = 298.16
      td = 298.16

c     initial canopy air conditions
      t12 = 298.16
      q12 = 0.01

c     initial value relative humidity of air at z12 (%)
      rh12 = 0.30

c     physical constants (SI)
      stef = 5.67e-8
      hvap = 2.5104e+6
      cw  = 4.18e+3
      cp = 1004.64
      r = 287.04
      pi = 3.141593
      grav = 9.80616
      densw = 1000.
      mdens = 41.4

c     data collection height above the surface (m)
      za = 65.

c     top and bottom heights of upper canopy (m)
      z1  = 40.
      z2  = 30.

c     zero-plane displacement height for upper canopy (m)
      da =  30. ! Carswell et al. (submitted)

c     roughness length parameter (m)
      zoa = 2.35   ! Shuttleworth (1988)
      zog = 0.005

c     friction velocity
      ufrica1 = 0.01
      ufric2g = 0.01
      u12     = 0.01

c     typical dimension of leaves and stems (m)
      du = 0.072
      ds = 0.100

c     leaf width (m)
      w  = 0.10

c     specific leaf area (m2 leaf kg-1 C)
      sl = 13.  !Medina & Cuevas (1996); Roberts et al. (1996)

c     current single-sided stem area index
      sai = 1.

c     initial values of carbon stored in reservoirs (kg C m-2)
      carbu =  0.360
      carbs = 18.000
      carbf =  0.100
      carbr =  0.750

c     initial co2 and o2 concentration (mol mol-1)
      co2a   = 0.000365
      o2conc = 0.210000

c     intrinsic quantum efficiency for c3 and c4 plant
      alpha = 0.060 !c3 - broadleaf trees

c     co2/o2 specificity ratio at 15 degrees C
      tau15 = 4500.0

c     o2/co2 kinetic parameters at 15 degrees C (mol/mol)
      kc15 = 1.5e-04
      ko15 = 2.5e-01

c     initial value intercellular co2 concentration
      ciub = 0

c     initial value leaf boundary layer co2 concentration
      csub = 0

c     initial value upper canopy stomatal conductance
      gsub = 0

c     leaf respiration coefficients
      gamaub = 0.0150  !broadleaf trees

c     'm' coefficients for stomatal conductance
      coefmub = 7.0   ! broadleaf trees

c     'b' coefficients for stomatal conductance
c     (minimum conductance when net photosynthesis is zero)
      coefbub = 0.010   !broadleaf trees

c     absolute minimum stomatal conductances
      gsubmin = 0.00001 !broadleaf trees

c     maximum values for ci (to avoid numerical instability)
      cimax = 2000.e-06

c     nominal values for vmax of top leaf at 15 C (mol co2 m-2 s-1)
      vmaxub = 85.e-06 ! broadleaf trees

c     initial intercepted liquid h2o on upper canopy leaf and stem area
c     (kg h2o m-2)
      wu = 0.
      ws = 0.

c     time constant for leaf and stem liquid drip (day converted to
c     seconds)
      tdrip = 0.5 * 86400.

c     initial values of downward and upward sensible heat flux (W m-2)
      fsh = 20.2
      fseng =	-18.5

c     numerical stability parameters
      dtairmax = 1.0
      dtvegmax = 2.5
      dqmax = 0.003
      hmin = -50.
      gmax = 20.
      gmin = -60.

c     thickness of soil layers (m)
      eg = 0.100
      ed = 4.900

c     density of soil minerals (kg/m3)
      densm = 2650.

c     specific heat of soil minerals and  organic matter(J kg-1 K-1)
      cm = 0.87e+3
      co = 1.92e+3

c     volume fraction of air, organic matter, field capacity  and
c     wilting point
      epsilonn  = 0.477  ! porosity
      fio       = 0.5
      thetafc   = 0.36
      thetawilt = 0.23

c     coefficient for soil moisture stress
      stressfac = -5.0

c     initial fraction of soil moisture
      thetag = 0.36
      thetad = 0.36

c     exponent of the moisture release equation (adimensional)
      b = 9.066

c     air entry water potential (mm)
      psie = 223.9

c     saturation hydraulic conductivity (mm s-1)
      ks =  0.0032
      kb =  1.e-6

c     thermal conductivity of water, air and soil minerals (W m-1 K-1)
      kw   = 0.596
      kair = 0.025
      km   = 2.50

c     weighted factor of water, air and soil minerals
      ew   = 0.92
      eair = 1.75
      em   = 0.54

c     heat capacity of leaves and stems
      cu =  2.109e+3
      cs =  2.109e+4

c     allocation fraction of net primary productivity
      au = 0.35  !0.45(inicial)
      as = 0.25  !0.40(inicial)
      af = 0.35  !0.10(inicial)
      ar = 0.05  !0.05(inicial)

c     residence time of living carbon (years converted to seconds)
      tauu =  1.00*365.*86400.
      tauf =  1.00*365.*86400.
      taus = 25.00*365.*86400.
      taur = 25.00*365.*86400.

c     fraction of carbon lost due to growth
      n = 0.30

c     initial value of mass of litter (kg C m-2)
      lu = 0.50  !0.50
      ls = 0.25  !0.25
      df = 0.45  !0.44
      dr = 0.15  !0.135

c     respiration rate of soil litter (micromol CO2 kg C-1 s-1 converted
c     to kg C kg C-1 s-1)
      hu = 0.230 * 12.e-9 * 4.5
      hs = 0.230 * 12.e-9 * 4.5
      hf = 0.230 * 12.e-9 * 3.7
      hr = 0.230 * 12.e-9 * 3.7

c     sapwood fraction of the total stem and coarse root biomass
      lambdasap = 0.10

