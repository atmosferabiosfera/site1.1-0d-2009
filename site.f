c     ------------------------------------------------------------------
      program main
c     ------------------------------------------------------------------

c     ==================================================================
c             initialization of variables
c     ==================================================================
      include 'initsite.h'

      write(*,*) 'Executing SITE.....................'

      if (landcov .eq. 2) include 'vegetation_02.h'

c     ir emissivity and absortivity of upper-stems
      ems = 1.-exp(-sai/(2.*avmuir))
      abss = ems

c     maximum intercepted water on upper canopy stem area(kg h2o/m2)
      wsmax = 0.1*sai

      lai = carbu*sl
      z12 = da + zoa
      
c     volume fraction of air at layers g and d, and soil minerals
      fiairg    = epsilonn - thetag
      fiaird    = epsilonn - thetad
      fim       = 1. - epsilonn

c     initial value of deep drainage (mm s-1)
      d = kb

c     saturation water content (m3 m-3)
      thetas = epsilonn

c     maximum stored water in layers g and d (kg h2o m-2)
      wgmax = eg*epsilonn*densw
      wdmax = ed*epsilonn*densw

c     water stored at field capacity, layers g and d (kg h20 m-2)
      wgfc = thetafc * eg * densw
      wdfc = thetafc * ed * densw

c     minimum stored water in layers g and d (kg h2o m-2)
      wgmin = eg*thetawilt*densw
      wdmin = ed*thetawilt*densw

c     stored water in layers g and d (kg h2o m-2)
      wg = thetag*eg*densw
      wd = thetad*ed*densw

c     available water content [0..1]
      awcg = min(1.,(wg-wgmin)/(wgfc-wgmin))
      awcd = min(1.,(wd-wdmin)/(wdfc-wdmin))

c     parameter for calculating maximum infiltration
      v = -b*psie/(1000.*eg)

c     organic matter mass (kg C m-2 converted to kg o.m. m-2)
      mo = (lu+ls)*2.

c     heat capacity of ground soil layer
      cg = densm*eg*fim*cm + mo*fio*co + cw*wg

c     heat capacity of deep soil layer
      cd = densm*ed*fim*cm + cw*wd

c     calculate heat conductivities of ground soil layer
      kappag = (thetag*ew*kw + fiairg*eair*kair + fim*em*km)/
     >         (thetag*ew + fiairg*eair + fim*em)

c     calculate heat conductivities of deep soil layer
      kappad = (thetad*ew*kw + fiaird*eair*kair + fim*em*km)/
     >	       (thetad*ew + fiaird*eair + fim*em)

      write(11,'(a4,5x,a4,3x,a4,6x,a3,8x,a3,6x,a3,8x,a4,7x,a4,7x,a4,
     >         7x,a5,6x,a4)') 'year','date','time','emu','ems','fup',
     >         'fira','firu','firs','fdown','firg'
      write(12,'(a4,5x,a4,3x,a4,5x,a3,7x,a3,9x,a3,8x,a3,7x,a3,8x,a3)')
     >         'year','date','time','fsr','fsu','fss','fsg','bsw','par'
      write(13,'(a4,5x,a4,3x,a4,6x,a7,2x,a4,7x,a4,6x,a5,8x,a6,7x,a6,7x,
     >         a5,8x,a6,7x,a6,12x,a7,6x,a2,11x,a3,10x,a2,7x,a7,5x,a5,
     >         8x,a5,9x,a5,8x,a5,8x,a3,10x,a3,10x,a3,10x,a3)') 'year',
     >         'date','time','densair','tavg','tava','csia1','psiha1',
     >         'psima1','csi2g','psih2g','psim2g','ufrica1','u1','u12',
     >         'u2','ufric2g','sigah','sigav','siggh','siggv','suh',
     >         'suv','ssh','ssv'
      write(14,'(a4,5x,a4,3x,a4,5x,a2,8x,a5,9x,a3,9x,a2,13x,a2,9x,a6,
     >         3x,a3,14x,a8,1x,a5,8x,a4,11x,a6,3x,a4,9x,a3,10x,
     >         a2,11x,a2,11x,a4,8x,a4,10x,a4,9x,a4)') 'year','date',
     >         'time','tf','rwork','tau','kc','ko','tempvm','prp',
     >         'stresstu','gamap','ciub','gbco2u','vmax','rub','je',
     >         'jc','agub','anub','csub','gsub'
      write(15,'(a4,5x,a4,3x,a4,4x,a4,9x,a4,9x,a4,11x,a5,6x,a5,3x,a5,7x,
     >         a2,10x,a5,7x,a2,9x,a5,8x,a5,8x,a2,11x,a3,9x,a3,10x,a3,
     >         10x,a3,11x,a5,8x,a5,8x,a2,11x,a2,11x,a5,8x,a5)') 'year',
     >         'date','time','prpu','prps','prpg','wumax','wsmax',
     >         'esatu','qu','esats','qs','esatg','qsatg','qg','etu',
     >         'eiu','fes','feu','fwetu','fwets','wu','ws','dripu',
     >         'drips'
      write(16,'(a4,5x,a4,3x,a4,4x,a3,10x,a3,8x,a2,9x,a2,8x,a6,6x,a6,
     >         7x,a4,9x,a3,9x,a3,10x,a3,10x,a3,10x,a3,10x,a3,10x,a3,
     >         11x,a3,8x,a3)')'year','date','time','q12',
     >         't12','ts','tu','esat12','qsat12','rh12','etu','eiu',
     >         'feu','fes','feg','fet','fea','flh','fsh'
      write(17,'(a4,5x,a4,3x,a4,4x,a2,9x,a2,9x,a4,11x,a2,7x,a2,11x,a2,
     >         15x,a6,7x,a6,7x,a6,7x,a6,7x,a2,11x,a2,8x,a4,9x,a4,8x,a2,
     >         11x,a2,11x,a2,11x,a4,9x,a6,7x,a2,11x,a1,12x,a4,9x,a4,9x,
     >         a4,9x,a4,10x,a2,9x,a5,6x,a5,4x,a2,9x,a5,6x,a5,9x,a4,7x,
     >         a4,4x,a3)')'year','date','time','tg','td','prpg','mo',
     >         'cg','cd','kappag','kappad','thetag','thetad','xg','xd',
     >         'psig','psid','kg','kd','fg','imax','runoff','fd','d',
     >         'mcxg','mcng','mcxd','mcnd','wg','wgmax','wgmin','wd',
     >         'wdmax','wdmin','awcg','awcd','fhs'
      write(18,'(a4,5x,a4,3x,a4,4x,a2,11x,a2,11x,a2,13x,a3,8x,a3,8x,a2,
     >         9x,a2,9x,a5,5x,a5,7x,a5,6x,a5,6x,a2,9x,a2,9x,a2,9x,a2,
     >         7x,a5,7x,a3,10x,a3,13x,a3)') 'year','date','time',
     >         'rr','rf','rs','ftg','ftd','gg','gd','carbu','carbs',
     >         'carbf','carbr','lu','ls','df','dr','rsoil','npp',
     >         'nee','lai'
      write(19,'(a4,5x,a4,3x,a4,4x,a2,11x,a5,11x,a1,8x,a1,11x,a2,10x,a3,
     >         8x,a3,9x,a3,8x,a4,7x,a3)')'year','date','time','Rn',
     >         'PARin','G','H','LE','NEE','NPP','GPP','APAR','LAI'
      write(19,'(24x,a4,9x,a14,2x,a4,5x,a4,8x,a4,8x,a10,1x,a10,2x,a8,3x,
     >         a8,3x,a5)')'Wm-2','micromolm-2s-1','Wm-2','Wm-2','Wm-2',
     >         'KgCha-1h-1','KgCha-1h-1','gCm-2h-1','MJm-2h-1','m2m-2'

c     ==================================================================
c             Reading  data
c     ==================================================================
      read(10,*)
      do i = 2,(ntl+1)
        read(10,*) year,date,time,ta,qa,ua,fsa,firatm,prp

c       ================================================================
c             Adjusting units
c       ================================================================
        prp = prp/dtime

        tga = tg
c       ================================================================
c             Canopy ir balance
c       ================================================================
c       ir emissivity and absortivity of upper-leaves
        emu = 1.-exp(-lai/(2.*avmuir))
        abu = emu
 	
c       calculate upward ir fluxes below tree level per overall area
        fup = stef*(tg**4)

c       ir balance calculate at layers u and s (W m-2)
        fira = firatm - emu*stef*(tu**4) - ems*stef*(ts**4)*(1-abu) -
     >         (1-abu)*(1-abss)*fup
	
c       calculate ir fluxes absorbed by upper leaves
        firu = abu*ems*stef*(ts**4)
     >         + abu*(1.-ems)*fup
     >         + abu*firatm
     >         - 2.*emu*stef*(tu**4)	
	
c       calculate ir fluxes absorbed by upper stems
        firs = abss*emu*stef*(tu**4)
     >         + abss*fup
     >         + abss*(1.-emu)*firatm
     >         - 2.*ems*stef*(ts**4)
	
c       calculate downward ir fluxes below tree level per overall area
        fdown = (1.-abu)*(1.-abss)*firatm 
     >          + emu*(1.-abss)*stef*(tu**4)
     >          + ems*stef*(ts**4)
	
c       calculate ir fluxes absorbed by upper soil
        firg  = fdown-stef*(tg**4)
	
c       *************results print************
        write(11,100) year,date,time,emu,ems,fup,fira,firu,firs,fdown,
     >                firg
  100   format(1f5.0,2f8.0,8f11.4)

c       ================================================================
c             Canopy radiation flux
c       ================================================================
c       calculate solar flux absorved by leaves
        if (fsa.le.0) then
          fsu = 0
        else
          fsu = (1. - reflu - transmu)*fsa
        end if
	
c       calculate solar flux absorved by stems
        if (fsa.le.0) then
          fss = 0
        else
          fss = transmu*fsa*(1. - transms)
        end if

c       calculate solar flux absorved by the soil
        if (fsa.le.0) then
          fsg = 0
        else
          fsg = transmu*transms*fsa
        end if

c       calculate solar flux reflected by leaves
        fsr = reflu*fsa

c       short wave balance calculate at layers u and s (W m-2)
        if (fsa.le.0) then
          bsw = 0
        else
          bsw = (1. - reflu - transmu*transms)*fsa
        end if

c       solar radiation balance calculate at layers u and s (W m-2)
        bsr = bsw + fira

c       photosynthetically active radiation
        if (fsa.le.0) then
          par = 0
        else
          par = fsa*0.4432      ! 0.4432 from Oliveira (2000);
        end if

c       ************results print*************
        write(12,200) year,date,time,fsr,fsu,fss,fsg,bsw,par
  200   format(1f5.0,2f8.0,6f11.4)
	
c       ================================================================
c             aerodynamics 
c       ================================================================
c	density of near-surface air	(kg m-3)
	densair = pre/(r*t12)

c       average temperature (K)	
	tavg = (t12 + tg)/2.
	tava = (ta + t12)/2.

	csia1 = -(0.41*grav*(da+zoa)*fsh)/(mdens*cp*tava*ufrica1**3)
        csia1 = max(min(3.,csia1),-3.) ! to ensure numerical stability
	  if (csia1 .ge. 0) then
c 	      write(*,*) 'stable atmosphere csia1 =',csia1
	      psiha1 = 6*alog(1 + csia1)
	      psima1 = psiha1
	  else
c 	      write(*,*) 'unstable atmosphere csia1 =',csia1
	      psiha1 = -2*alog((1 + (1 - 16*csia1)**0.5)/2)
	      psima1 = 0.6*psiha1
	  end if
	
  	csi2g = -(0.41*grav*zog*fseng)/(mdens*cp*tavg*ufric2g**3)
        csi2g = max(min(3.,csi2g),-3.)  ! to ensure numerical stability
        if (csi2g .ge. 0) then
c 	      write(*,*) 'stable atmosphere csia2 =',csi2g
	      psih2g = 6*alog(1+csi2g)
	      psim2g = psih2g
	  else
c 	      write(*,*) 'unstable atmosphere csia2 =',csi2g
	      psih2g = -2*alog((1 + (1 - 16*csi2g)**0.5)/2)
	      psim2g = 0.6*psih2g
	  end if

c       calculation  of attenuation coefficient of wind
        lm = (6*(w**2)*z1/(pi*lai))**(0.333333)
        a  = sqrt(0.2*lai*z1/lm)
	
        ufrica1 = (0.41*ua)/(alog((za - da)/zoa) + psima1)
        u1 = max(0.1,ufrica1/0.41*(alog((z1 - da)/zoa) + psima1))
        u12 = max(0.1,u1*exp(a*(z12/z1-1)))
	u2 = max(0.1,u1*exp(a*(z2/z1-1)))
        ufric2g = (0.41*u2)/(alog(z2/zog) + psim2g)

c       calculate transfer coefficient between air-air, soil-air (m/s)
c       (Campbell e Norman, 1998) 
        ea = ((alog((za-da)/zoa)+psima1)*(alog((za-da)/zoa)+psiha1))
        if (ea .eq. 0.) ea = 0.1
        sigah = ua*(0.41**2)/ea
        sigav = 0.622*ua*(0.41**2)/ea
        e12 =  ((alog(z12/zog)+psim2g)*(alog(z12/zog)+psih2g))
        siggh = u12*(0.41**2)/e12
        siggv = 0.622*u12*(0.41**2)/e12

c       air-vegetation transfer coefficients for upper canopy
c       leaves (m s-1)
c	(0.029 or 0.018)/densair is a conversion factor from mol m-2 s-1
c       to m s-1
	suh = max(0.135*sqrt(u12/du),0.159)*0.029/densair
	suv = max(0.147*sqrt(u12/du),0.173)*0.018/densair
        ssh = max(0.135*sqrt(u12/ds),0.159)*0.029/densair
	ssv = max(0.147*sqrt(u12/ds),0.173)*0.018/densair

c       ************results print*************
         write(13,300) year,date,time,densair,tavg,tava,csia1,psiha1,
     >                 psima1,csi2g,psih2g,psim2g,ufrica1,u1,u12,u2,
     >                 ufric2g,sigah,sigav,siggh,siggv,suh,suv,ssh,ssv
  300   format(1f5.0,2f8.0,3f11.4,6E13.4,4f13.4,9E13.4)
	
c       ================================================================
c             canopy	physiology 
c       ================================================================
c       calculate physiological parameter values which are a function
c       of temperature
        tf = tu - 273.16
        tf = max(0.,tf)
        rwork = max(0.,3.47e-03 - 1./tf)
        tau = tau15*exp(-5000.0*rwork)
        kc  = kc15 *exp( 6000.0*rwork)
        ko  = ko15 *exp( 1400.0*rwork)
        tempvm = exp(3500.0*rwork)/
     >           ((1.0 + exp(0.40*(5.0 - tf)))*
     >           (1.0 + exp(0.40*(tf - 50.0))))

c       soil moisture stress factor for the upper canopy
        stresstu = (1.-exp(stressfac*awcd))/(1.-exp(stressfac))
	
c       upper canopy gamap values (mol/mol)
        gamap = o2conc/(2.*tau)
	
c       constrain ci values to acceptable bounds - to help ensure
c       numerical stability

        ciub = max (1.05*gamap, min(cimax,ciub))
	
c       calculate boundary layer parameters (mol CO2 m-2 s-1) 
        gbco2u = min (10.0, max(0.1,0.110*sqrt(u12/du))) 
	
c       vmax and dark respiration for current conditions
        vmax = vmaxub*tempvm*stresstu
        rub  = gamaub*vmax
	
c       'light limited' rate of photosynthesis (mol m-2 s-1)
        if (par.le.0) then
          je = 0
        else
          je = fapar*par*4.59e-06*alpha*(ciub - gamap)/
     >       (ciub + 2.*gamap)
        end if
	
c       'rubisco limited' rate of photosynthesis (mol m-2 s-1)
        jc = vmax*(ciub - gamap)/
     >       (ciub + kc*(1. + o2conc/ko))
	
c       calculate the net photosynthesis rate (mol m-2 s-1)
        agub = min(je,jc)
        anub = agub - rub
	
c       calculate co2 concentrations and stomatal condutance 
c	
c       calculate new value of cs using implicit scheme
        csub = 0.5*(csub+co2a-anub/gbco2u)
        csub = max(1.05*gamap,csub)
c
c       calculate new value of gs using implicit scheme
        gsub = 0.5*(gsub+(coefmub*anub*rh12/csub + coefbub * stresstu))
        gsub = max (gsubmin, coefbub*stresstu, gsub)
c
c       calculate new value of ci using implicit scheme
        ciub = 0.5 * (ciub + csub - 1.6 * anub/gsub)
        ciub = max (1.05 * gamap, min (cimax, ciub))
	
c       ************results print*************
        write(14,400) year,date,time,tf,rwork,tau,kc,ko,tempvm,prp,
     >                stresstu,gamap,ciub,gbco2u,vmax,rub,je,jc,agub,
     >                anub,csub,gsub

  400   format(1f5.0,2f8.0,1f11.4,1E13.4,1f13.4,1E13.4,2f11.4,1E13.4,
     >         1f13.4,2E13.4,1f11.4,8E13.4)

c       ================================================================
c             intercepted water
c       ================================================================
c       rates of precipitation intercepted by leaves and stems and
c       throughfall
        prpu = (prp*(1-exp(-0.50*lai)))
        prps = ((prp-prpu)*(1-exp(-0.50*sai)))
	prpg = prp - prpu - prps
	
c       maximum intercepted water on upper canopy leaf area(kg h2o/m2)
        wumax = 0.1*lai

c       fraction of upper canopy leaf area wetted by intercepted liquid
        fwetu = min(0.8,wu/wumax)

c       fraction of upper canopy stem area wetted by intercepted liquid
        fwets = min(0.8,ws/wsmax)
        
c       saturation specific humidities for leaves,stem and soil
        esatu = 610.78*10**((7.5*(tu-273.16))/(237.3 + (tu - 273.16)))
        qu    = 0.622*esatu/(pre - 0.378*esatu)
	esats = 610.78*10**((7.5*(ts-273.16))/(237.3 + (ts - 273.16)))
        qs    = 0.622*esats/(pre - 0.378*esats)
	esatg = 610.78*10**((7.5*(tg-273.16))/(237.3 + (tg - 273.16)))
	qsatg = 0.622*esatg/(pre - 0.378*esatg)

c	specific humidity of soil air
	qg = min(1.0,(1.-exp(-(thetag/thetas)**2.))/(1.-exp(-1.)))*qsatg

c	evapotranspiration of leaves and stems
c
c       leaf transpiration
        etu = max(0.,(qu-q12)*lai*densair*(1-fwetu)/
     >           (1/suv + densair/(0.029*gsub)))*awcd

c	evaporation of intercepted water of leaf and stem
        if (qu.ge.q12) then ! evaporation will happen
          eiu = min(fwetu*suv*(qu-q12)*lai*densair,0.98*wu/dtime)
        else                ! condensation will happen
          eiu = max(suv*(qu-q12)*lai*densair,(wu-wumax)/dtime)
        end if

        if (qs.ge.q12) then ! evaporation will happen
          fes = min(fwets*ssv*(qs-q12)*sai*densair,0.98*ws/dtime)
        else                ! condensation will happen
          fes = max(ssv*(qs-q12)*sai*densair,(ws-wsmax)/dtime)
        endif

        feu = etu + eiu
	
c	calculate intercepted liquid h2o on upper canopy leaf area
c       (kg m-2)
        wu = wu + (prpu - eiu)*dtime
	if (wu .gt. wumax) then
	  dripu = (wu - wumax)/dtime
	else
	  dripu = wu/tdrip
	end if
	wu = wu - dripu*dtime
	
c       calculate intercepted liquid h2o on upper canopy stem area
c       (kg m-2)
        ws = ws + (prps - fes)*dtime
	if (ws .gt. wsmax) then
	  drips = (ws - wsmax)/dtime
   	else
	  drips = ws/tdrip
	end if
        ws = ws - drips*dtime

c       fraction of upper canopy leaf and stem area wetted by
c       intercepted liquid
        fwetu = min(0.8,wu/wumax)
        fwets = min(0.8,ws/wsmax)

	if (fwetu.lt.0.1) then      ! to avoid underflow
	  fwetu = 0.00
	  dripu = wu/dtime
	  wu = 0.00
	endif
        if (fwets.lt.0.1) then      ! to avoid underflow
	  fwets = 0.00
	  drips = ws/dtime
	  ws = 0.00
	endif

c       ************results print*************
        write(15,500) year,date,time,prpu,prps,prpg,wumax,wsmax,esatu,
     >                qu,esats,qs,esatg,qsatg,qg,etu,eiu,fes,feu,fwetu,
     >                fwets,wu,ws,dripu,drips

  500   format(1f5.0,2f8.0,3E13.4,3f11.4,1E13.4,1f11.4,1E13.4,1f11.5,
     >	       12E13.4)
	
c       ================================================================
c	        mass and energy conservation
c       ================================================================
c       first solution of the system using latent heat fluxes estimates 
c       based on the temperatures of the previous time step
        a11 = 1+densair*cp*suh*2.*lai*dtime/(cu+cw*wu)
        a13 =  -densair*cp*suh*2.*lai*dtime/(cu+cw*wu)
        a22 = 1+densair*cp*ssh*2.*sai*dtime/(cs+cw*ws)
        a23 =  -densair*cp*ssh*2.*sai*dtime/(cs+cw*ws)
        a33 = -siggh*densair*cp*dtime/cg
        a34 = 1+(siggh*densair*cp*dtime/cg)+(kappag*dtime/cg)
        a35 = -siggh*densair*hvap*dtime/cg
        a41 = 2.*lai*suh
        a42 = 2.*sai*ssh
	a43 = sigah-siggh-(2.*lai*suh)-(2.*sai*ssh)
        a44 = siggh
        a55 = sigav+siggv
 
	b1 = tu + (fsu+firu-hvap*feu)*dtime/(cu + cw*wu)
        b2 = ts + (fss+firs-hvap*fes)*dtime/(cs + cw*ws)
        b3 = tg + (((fsg+firg) - (siggh*densair*hvap*qg) + 
     >             (kappag*td))*dtime)/cg
        b4 = sigah*ta
        b5 = fes + feu + siggv*qg + sigav*qa

	c1 = (a43/a41-a13/a11)*a41/a42-a23/a22
        if (c1 .eq. 0.) c1 = 0.1
        c2 = a44/a42
        c3 = (b4/a41-b1/a11)*a41/a42-b2/a22
	  
        d13 = a13/a11
	d16 = b1/a11
        d23 = a23/a22
        d26 = b2/a22
        d34 = a34/a33
        d35 = a35/a33
        d36 = b3/a33
        d45 = -a35/(a33*((c2/c1)-(a34/a33)))
        d46 = ((c3/c1)-(b3/a33))/((c2/c1)-(a34/a33))
        d56 = b5/a55
	   
c       upper canopy-air mass balance;
c       calculate upper air specific humidity
        q12 = d56
        if ((q12-qa) .gt. dqmax) then
          q12 = qa+dqmax
        else if (qa .gt. q12) then
          q12 = qa
        endif
	   
c       soil surface energy balance; calculate soil surface temperature
        tg = d46 - d45*q12
        if ((tg-ta) .gt. dtvegmax) then
          tg = ta+dtvegmax
        else if ((ta-tg) .gt. dtvegmax) then
          tg = ta-dtvegmax
        endif
	   
c       upper canopy-air energy balance; calculate upper air temperature
        t12 = d36 - d34*tg - d35*q12
        if ((t12-ta) .gt. dtairmax) then
          t12 = ta+dtairmax
        else if ((ta-t12) .gt. dtairmax) then
          t12 = ta-dtairmax
        endif
	   
c       upper stem energy balance; calculate upper stem temperature
        ts = d26 - d23*t12
        if ((ts-ta) .gt. dtvegmax) then
          ts = ta+dtvegmax
        else if ((ta-ts) .gt. dtvegmax) then
          ts = ta-dtvegmax
        endif
	   
c       upper leaf energy balance; calculate upper leaf temperature
        tu = d16 - d13*t12
        if ((tu-ta) .gt. dtvegmax) then
          tu = ta+dtvegmax
        else if ((ta-tu) .gt. dtvegmax) then
          tu = ta-dtvegmax
        endif

c       calculate the relative humidity in the canopy air space
c       correct q12 for eventual condensation
c       set a minimum value of 0.30 to avoid errors in the
c       physiological calculations
c       warning: some mass of water may be lost in this process
        esat12 = 610.78*10**((7.5*(t12-273.16))/
     >                       (237.3 + (t12 - 273.16)))
        qsat12 = 0.622*esat12/(pre - 0.378*esat12)

        rh12   = max(0.30, q12/qsat12)
        if (q12.gt.qsat12) then
          q12 = qsat12
          rh12 = 1.00
        endif

c       recalculation of latent heat fluxes leaf transpiration
        etu = max(0.,(qu-q12)*lai*densair*(1-fwetu)/
     >               (1/suv + densair/(0.029*gsub)))*awcd

c	evaporation of intercepted water of leaf and stem
        if (qu.ge.q12) then ! evaporation will happen
          eiu = min(fwetu*suv*(qu-q12)*lai*densair,0.98*wu/dtime)
        else                ! condensation will happen
          eiu = max(suv*(qu-q12)*lai*densair,(wu-wumax)/dtime)
        endif

        if (qs.ge.q12) then ! evaporation will happen
          fes = min(fwets*ssv*(qs-q12)*sai*densair,0.98*ws/dtime)
        else                ! condensation will happen
          fes = max(ssv*(qs-q12)*sai*densair,(ws-wsmax)/dtime)
        endif

        feu = etu + eiu
        
c       second solution of the system using latent heat fluxes estimates 
c       based on the temperatures of this time step
        a11 = 1+densair*cp*suh*2.*lai*dtime/(cu+cw*wu)
        a13 =  -densair*cp*suh*2.*lai*dtime/(cu+cw*wu)
        a22 = 1+densair*cp*ssh*2.*sai*dtime/(cs+cw*ws)
        a23 =  -densair*cp*ssh*2.*sai*dtime/(cs+cw*ws)
        a33 = -siggh*densair*cp*dtime/cg
        a34 = 1+(siggh*densair*cp*dtime/cg)+(kappag*dtime/cg)
        a35 = -siggh*densair*hvap*dtime/cg
        a41 = 2.*lai*suh
        a42 = 2.*sai*ssh
	a43 = sigah-siggh-(2.*lai*suh)-(2.*sai*ssh)
        a44 = siggh
        a55 = sigav+siggv
 
	b1 = tu + (fsu+firu-hvap*feu)*dtime/(cu + cw*wu)
        b2 = ts + (fss+firs-hvap*fes)*dtime/(cs + cw*ws)
        b3 = tg + (((fsg+firg) - (siggh*densair*hvap*qg) + 
     >             (kappag*td))*dtime)/cg
        b4 = sigah*ta
        b5 = fes + feu + siggv*qg + sigav*qa

	c1 = (a43/a41-a13/a11)*a41/a42-a23/a22
        if (c1 .eq. 0.) c1 = 0.1
        c2 = a44/a42
        c3 = (b4/a41-b1/a11)*a41/a42-b2/a22
	  
        d13 = a13/a11
	d16 = b1/a11
        d23 = a23/a22
        d26 = b2/a22
        d34 = a34/a33
        d35 = a35/a33
        d36 = b3/a33
        d45 = -a35/(a33*((c2/c1)-(a34/a33)))
        d46 = ((c3/c1)-(b3/a33))/((c2/c1)-(a34/a33))
        d56 = b5/a55
	   
c       upper canopy-air mass balance;
c       calculate upper air specific humidity
        q12 = d56
        if ((q12-qa) .gt. dqmax) then
          q12 = qa+dqmax
        else if (qa .gt. q12) then
          q12 = qa
        end if
	   
c       soil surface energy balance; calculate soil surface temperature
        tg = d46 - d45*q12
        if ((tg-ta) .gt. dtvegmax) then
          tg = ta+dtvegmax
        else if ((ta-tg) .gt. dtvegmax) then
          tg = ta-dtvegmax
        end if
	   
c       upper canopy-air energy balance; calculate upper air temperature
        t12 = d36 - d34*tg - d35*q12
        if ((t12-ta) .gt. dtairmax) then
          t12 = ta+dtairmax
        else if ((ta-t12) .gt. dtairmax) then
          t12 = ta-dtairmax
        end if
	   
c       upper stem energy balance; calculate upper stem temperature
        ts = d26 - d23*t12
        if ((ts-ta) .gt. dtvegmax) then
          ts = ta+dtvegmax
        else if ((ta-ts) .gt. dtvegmax) then
          ts = ta-dtvegmax
        endif
	   
c       upper leaf energy balance; calculate upper leaf temperature
        tu = d16 - d13*t12
        if ((tu-ta) .gt. dtvegmax) then
          tu = ta+dtvegmax
        else if ((ta-tu) .gt. dtvegmax) then
          tu = ta-dtvegmax
        endif

c       calculate the relative humidity in the canopy air space
c       correct q12 for eventual condensation
c       set a minimum value of 0.30 to avoid errors in the
c       physiological calculations
c       warning: some mass of water may be lost in this process
c       however, this two-step process may reduce the loss of mass
c       to a minimum
        esat12 = 610.78*10**((7.5*(t12-273.16))/
     >                       (237.3 + (t12 - 273.16)))
        qsat12 = 0.622*esat12/(pre - 0.378*esat12)
        rh12   = max(0.30, q12/qsat12)
        if (q12.gt.qsat12) then
          q12 = qsat12
          rh12 = 1.00
        endif

c       flux of water from the soil (kg h2o m-2 s-1)
	feg = siggv * densair * (qg - q12)*awcg

c       total water flux calculated from two different methods

c       total water flux (kg H2O m-2 s-1)
        fet   = feu + fes + feg

c       total water flux to atmosphere (kg H2O m-2 s-1)
        fea = sigav * densair * (q12 - qa)
        
c       total flux of latent heat (W m-2)
        flh = fet * hvap

c       flux of sensible heat from layer g to layer 12 (W m-2)
        fseng = densair*cp*siggh*(tg-t12)

c       flux of sensible heat from leaves and stems to layer 12 (W m-2)
        fsenu = densair*cp*suh*lai*(tu-t12)
        fsens = densair*cp*ssh*sai*(ts-t12)

c       total flux of sensible heat (W m-2)
        fsh = fsenu + fsens + fseng
        
        if (fsh .lt. hmin) then
          ts = ts - (fsh - hmin)*dtime/cs
          fsh = hmin
        end if
        
c       ************results print*************
        write(16,600) year,date,time,q12,t12,ts,tu,esat12,qsat12,
     >                rh12,etu,eiu,feu,fes,feg,fet,fea,flh,fsh
  600   format(1f5.0,2f8.0,1E13.4,4f11.4,9E13.4,2f11.4)
	   
c       ================================================================
c             soil physics
c       ================================================================

c	balance of water in the soil

	prpg = prpg + dripu + drips    ! (kg H2O m-2 s-1)

        runoff = 0.   !superficial drainage (kg H2O m-2 s-1)

        mcxg = 0.     !mass correction for layer g
        mcxd = 0.     !mass correction for layer d
        mcng = 0.     !mass correction for layer g
	mcnd = 0.     !mass correction for layer d
		  	
c	volumetric fraction of water and hydraulic conductivity
	thetag = wg/(eg*densw)
	thetad = wd/(ed*densw)

	xg = thetag/thetas
	xd = thetad/thetas
	
	psig = psie*xg**(-b)
	psid = psie*xd**(-b)
	
	kg = ks*xg**(2*b+3)
	kd = ks*xd**(2*b+3)
	
c	calculate maximum infiltration (mm/s)
        imax = ks*(v*xg+1-v)

c	calculate infiltration and runoff
	fg = min(imax, prpg) ! fg > 0 for downward flux
        if (prpg .gt. imax) then
	  runoff = runoff + (prpg - imax)
	end if

c       (Bonan, 1996) eg and ed converted from m to mm
c       fd > 0 for downward flux
        fd = (-(2*(psig - psid) + (eg + ed)*1000)/((eg*1000/kg) +
     >       (ed*1000/kd)))

c       calculate stored water in soil g layer
	wg = wg + dtime*(fg - fd - feg)

        if (wg .gt. wgmax) then
          mcxg = wg - wgmax
          if (wd .lt. wdmax) then
            if (mcxg .le. (wdmax-wd)) then
              wd = wd + mcxg
              mcxg = 0
            else
              wd = wdmax
              mcxg = mcxg - (wdmax-wd)
            end if
          end if
            wg = wgmax
	end if

	if (wg .lt. wgmin) then
	  mcng = wg - wgmin
	  wg = wgmin
        end if

        fd = fd + mcng/dtime

c       calculate stored water in soil d layer
        wd = wd + dtime*(fd-d-etu)

        if (wd .gt. wdmax) then
          mcxd = wd - wdmax
     	  if (mcxd .gt. (wgmax - wg)) then
            mcxd = mcxd - (wgmax - wg)
            wg = wgmax
          else
            wg = wg + mcxd
            mcxd = 0
          end if
            wd = wdmax
	end if

        if (wd .lt. wdmin) then
	  mcnd = wd - wdmin
	  wd = wdmin
        end if

c	flow across the bottom layer
        d = kd + mcnd/dtime

        runoff = runoff + ((mcxg+mcxd)/dtime)

        awcg = min(1.,(wg-wgmin)/(wgfc-wgmin))
        awcd = min(1.,(wd-wdmin)/(wdfc-wdmin))

c       mass of organic matter on the ground (kg m-2)
        mo = (lu+ls)*2.

c       volume fraction of air at layers g and d
        fiairg    = epsilonn - thetag
        fiaird    = epsilonn - thetad
	
c       heat capacity of ground soil layer
        cg = densm*eg*fim*cm + mo*fio*co + cw*wg
	
c       heat capacity of deep soil layer
        cd = densm*ed*fim*cm + cw*wd
	
c       calculate heat conductivities of ground soil layer
        kappag = (thetag*ew*kw + fiairg*eair*kair + fim*em*km)/
     >           (thetag*ew + fiairg*eair + fim*em)
	
c       calculate heat conductivities of deep soil layer
        kappad = (thetad*ew*kw + fiaird*eair*kair + fim*em*km)/
     >           (thetad*ew + fiaird*eair + fim*em)
	
c       calculate temperature of deep soil layer
        tda = td
        td = td + (dtime/cd)*(tg-td)/((eg/(2*kappag))+(ed/(2*kappad)))
          
c       calculate heat flux in the soil
        fhs = cg*(tga-tg)/dtime + cd*(tda-td)/dtime

        if (fhs .lt. gmin) then
          tg = tg + (fhs - gmin)*dtime/cg
          fhs = gmin
        end if
        
        if (fhs .gt. gmax) then
          tg = tg + (fhs - gmax)*dtime/cg
          fhs = gmax
        end if

c       ************results print*************
	write(17,700) year,date,time,tg,td,prpg,mo,cg,cd,kappag,kappad,
     >                thetag,thetad,xg,xd,psig,psid,kg,kd,fg,imax,
     >                runoff,fd,d,mcxg,mcng,mcxd,mcnd,wg,wgmax,wgmin,
     >                wd,wdmax,wdmin,awcg,awcd,fhs
  700   format(1f5.0,2f8.0,2f11.4,1E13.4,1f11.4,2E13.4,8f13.4,7E13.4,
     >         4E13.4,1f11.4,8f11.4 )
	
c       ================================================================
c             carbon balance
c       ================================================================
c       Calculate respiration terms a function of temperature at 15
c       degree C following respiration parametrization of Lloyd and
c       Taylor
        rwork1 = exp (3500.0 * (1./288.16 - 2./(tg+td)))
        rwork2 = exp (3500.0 * (3.47e-03 - 1./ts))
	
c       respiration of fine roots, stems and coarse roots
        rf = (3*1.00)/(365*86400)*(carbf*(1000.0/12.0))*rwork1
        rs = (3*0.01)/(365*86400)*(carbs*(1000.0/12.0))*rwork2
     >	     *lambdasap
        rr = (3*0.01)/(365*86400)*(carbr*(1000.0/12.0))*rwork1
     >	     *lambdasap
	  	
c       calculate of net primary productivity (NPP)
c       12e-3 is a conversion factor from mol CO2 to kg C
        npp = (1 - n)*(agub - rub - rs - rf - rr)*12.e-3

c       function of soil temperature
        ftg = 2**(((tg-273.16)-10.)/10.)
        ftd = 2**(((td-273.16)-10.)/10.)

c       function of soil moisture
        gg = max(min(0.25 + 0.75*(wg/wgmax),1.),0.)
        gd = max(min(0.25 + 0.75*(wd/wdmax),1.),0.)

c       carbon stored in leaves
        carbu = carbu + (au*npp - carbu/tauu)*dtime

c       new leaf area index
        lai = carbu*sl
    
c       carbon stored in stems
        carbs = carbs + (as*npp - carbs/taus)*dtime

c       carbon stored in fine roots
        carbf = carbf + (af*npp - carbf/tauf)*dtime

c       carbon stored in dense roots
        carbr = carbr + (ar*npp - carbr/taur)*dtime

c       soil carbon decomposition processes

c       soil leaf litter decomposition
        lu = lu + (carbu/tauu - ftg*gg*hu*lu)*dtime

c       soil stem litter decomposition
        ls = ls + (carbs/taus - ftg*gg*hs*ls)*dtime

c       soil fine root litter decomposition
        df = df + (carbf/tauf - (ftg + ftd)*(gg + gd)/4*hf*df)*dtime

c       soil dense root litter decomposition
        dr = dr + (carbr/taur - (ftg + ftd)*(gg + gd)/4*hr*dr)*dtime

c       respiration of soil
        rsoil= (ftg*gg*hu*lu) + (ftg*gg*hs*ls) +
     >         ((ftg + ftd)*(gg + gd)/4*hf*df) +
     >         ((ftg + ftd)*(gg + gd)/4*hr*dr)

c       net ecosystem exchange  (NEE)
        nee = rsoil - npp   ! nee < 0 during the day

c       ************results print*************
	write(18,800) year,date,time,rr,rf,rs,ftg,ftd,gg,gd,
     >                carbu,carbs,carbf,carbr,lu,ls,df,dr,rsoil,npp,
     >                nee,lai
  800   format(1f5.0,2f8.0,3E13.4,12f11.4,3E13.4,1f11.4)

C        if (par .le. 0.) then
C	  lue = -1.0
C	else
C	  lue = agub*12.e6/(fapar*par) !(g C MJ-1)
C	end if

	gpp  = agub*12.*dtime        !(g C m-2 h-1)
        apar = fapar*par*(dtime/1e6) !(MJ m-2 h-1)
        npp  = npp*dtime*1e4         !(kg C ha-1 h-1)
        nee  = nee*dtime*1e4         !(kg C ha-1 h-1)
        par  = par*10/2.35           !(micromol m-2 s-1)

        write(19,900) year,date,time,bsr,par,fhs,fsh,flh,nee,npp,gpp,
     >		      apar,lai
  900   format(1f5.0,2f8.0,2f11.4,6x,8f11.4)

      end do ! i

      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      
      write(*,*) '...................................'
      write(*,*) 'End of the process!'

      end
