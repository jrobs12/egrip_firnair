c==============================================================================
C	Model of firn gravitational and thermal fractionation
C       driven by seasonal temperature fluctuations or glacial terminations
C		with lock-in zone gas model 
C------------------------------------------------------------------------------
c     Modified by Anais 2010-01-13
c     firn_allgas8c.f identical to 8b except do the bb loop starting from the bottom.
c     do 1200 bb=jbmax,jbmin, -1
c     Also, get the run ready to do 1e-6 inst of 1e-8.
c    --- firn_allgas12.f---
c     change the formulation of uscale to an arc tan
c     firn_allgas15i.f for isotopes, works
c     firn_allgas15j.f try to use the same file, but for gases
c modified it to interpolate temperature and diffusivity profiles
c     firn_allgas16b.f : incorporates temperature model.
c     firn_allgas16iso.f : same, including the isotopes
c - make sure it doesnt' start interpolating temperature before the start of the temperature time series
c - fix the timestep in the temp model when we want to use ktemp: timestep needs to be ktemp*tstep.
c - re-initialise cprint=.true. inside the ng loop, otherwise, we don't get the csave routine to work for multiple gases.
c - adjust seasonal cycle
c -fiz csave: out of that kprint if statement that was making it useless
c firn_allgas16isot02_N8g_diff.f : 8g = with inert gases, and co2, ch4, ch3ccl3
c                                  diff: for the diffusivity optimisation
c                                       makes the amplitude proportional to diff0.
C---------------------------------------------------------------------
C       Jeff Severinghaus 12-22-98 Scripps Institution of Oceanography
C			Modified 10-31-2004 for lock-in zone evolution
C			Modified 11-18-2006 for 13CO2 correction with
c			  fixed	bottom box and linear taper of eddy diffusivity
c				near 25-40 m depth.
c			Modified 8-15-2007 to include exponential gravitational
c			  settling
c			Modified Nov-4-2008 to make dvbubble parameterization consistent
c			  with Severinghaus and Battle (2006).
c			Expulsion of air from lock-in disabled. 12-5-08
c			Modified 6-4-2009 to allow eddy diffusion in lock-in zone
c			Modified 7-22-2009 to include porosity gradient in lock-in zone
c			Modified 7-21-2010 Removed porosity gradient from lock-in zone and
c				changed air advection scheme in lock-in to airspeed-based temporal grid
c			Modified 7-23-2010 to run multiple gases, 2-D arrays
C     Explanation:
C       This is finite difference scheme for the time-dependent 
c       diffusion equation. It models effects of thermal diffusion and
C       gravity on firn air isotopic composition as modulated by temperature
c       transients during seasonal cycles or other temperature changes.
c       Fluxes are calculated according to a relaxation to the equilibrium
c       gradient (as in Bender et al., 1994): Flux = -D[d(C-Ceq)/dz] 
c       where D is diffusivity, C is concentration, z is depth, and Ceq is the 
c       equilibrium concentration including gravitational (Sowers et al.,1989)
c       and thermal diffusion (Grew et al., 1952) effects.      
c
c	There is a separate heat diffusion model, based on Alley and Koci (1979),
c	that sets the temperature profile before the gas model runs, at each time
c	step.  The models are not coupled; the output of the
c	gas model does not affect the heat model.
c
c	There is an uneven grid spacing in the lock-in zone, called the "temporal
c	grid", to improve the treatment of advection (advection is simulated by
c	moving the concentrations down by one grid cell per unit time, where unit
c	time corresponds to the temporal spacing of the grid points.  This unit
c   time is called "bundle" in the model and it has units of years.  It 
c	equals the time that the air takes to move from one gridpoint to the
c	next within the lock-in zone.  The grid spacing is calculated from the
c	air advection velocity, in turn derived from accumulation rate, density,
c	and the bubble closure parameterization)
c
C     References:
c	Alley, R.B., and Koci, B.R., Recent warming in central Greenland? Ann. 
c	  Glaciol. 14: (1990).
c       Bender, M.L., Sowers, T., Barnola, J.-M., and Chappellaz, J., 1994.
c         Changes in the O2/N2 ratio of the atmosphere during recent decades 
c         reflected in the composition of air in the firn at Vostok Station,
c         Antarctica, Geophysical Research Letters 21:189-192.
c       Grew, K.E., and Ibbs, T.L., 1952.  Thermal diffusion in gases, 
c         Cambridge University Press.
c       Reid, R.C., Prausnitz, J.M., and Sherwood, T.K., 1977. The properties
c         of Gases and Liquids, Third Edition, New York: McGraw-Hill.
c       Schwander, J., Stauffer, B., and Sigg, A., 1988. Air mixing in firn and
c         the age of the air at pore close-off, Annals of Glaciology 10:141-145
c       Schwander, J., Barnola, J.-M., Andrie, C., Leuenberger, M., Ludin, A.,
c         Raynaud, D., and Stauffer, B., 1993. The Age of the Air in the Firn 
c         and the Ice at Summit, Greenland, J. Geophys. Res. 98:2831-2838.
c       Sowers, T., Bender, M., and Raynaud, D., 1989. Elemental and isotopic
c         composition of occluded O2 and N2 in polar ice, JGR 94:5137-5150.
c       Yen, Y.-C., 1981. Review of thermal properties of snow, ice, and sea
c         ice, CRREL Report 81-10, US Army Corps of Engineers, Cold Regions 
c         Research and Engineering Laboratory, Hanover, N.H.
c
c     Accompanying files:
c	'inputdata_Megadunes_lockin.txt'  input file for parameters
c	'diffusivityMegadunes.txt'		input file for CO2-tuned effective diffusivity
c
c     Units: 
c	length        meters
c       concentration per mil with respect to atmosphere
c       time          s  (except sometimes I use years, as noted)
c       mass          kg (except mass difference which is in g/mol)
cc
c     Major variables (vector quantities positive downwards):
c       t(i)          temperature in box i, Celsius
c       c(n,i)        concentration or mixing ratio of gas n in box i
c		cl(n,i)		  concentration in the lock-in zone

c        [For example, c(5,7) would be the d15N (gas 5) in the sixth box
c         down from the surface (the atmosphere itself is box 1, with a depth
c         of zero.)]
c
c	  Other variables:
c	z(i)	      depth below surface, m
c	depth(i)	  depth below surface for lock-in zone model
c	deriv(i,j)    time derivative of gas concentration in box i
c	derivlock(i,j)  same for lock-in zone
c	diff(i)	      effective diffusivity of firn, found empirically from CO2 data (units m2/d in file)
c   D(j)		  free-air molecular diffusivity of gas j, expressed as a ratio to that of CO2
c	s(i)	      porosity of firn, m3/m3
c   stotal(i)	  same for lock-in zone
c	sc(i)	      closed porosity of firn, which is a parameterized function of density
c	sclosed(i)	  same for lock-in zone
c	openpor(i)	  open porosity of firn (s(i) - sc(i))
c	sopen(i)	  same for lock-in zone
c	vb(i)	      volume of air trapped in bubbles closed off in box i, m3/m2 (corrected for bubble compression)
c	vbl(i)		  same for lock-in zone
c	q	          concentration expressed as delta/1000 + 1
c	dvbubble      change in volume of trapped air in bubbles (normalized to surface pressure, i.e. corrected for bubble compression assuming bubbles shrink at same proportional rate as bulk firn specific volume)
c	dvbubblel	  same for lock-in zone
c
c     Parameters:
c		ng			  number of gases treated in one run
c       ldimlevels    number of diffusive zone gas domain boxes including surface box, also called dimlevels
c       levels        number of gas domain boxes, not including surface box
c		lt			  number of heat model boxes including surface box, also called tlevels (lt*dz must be larger than total firn thickness) 
c		layers		  number of gas boxes in the lock-in zone model (these boxes are not equally-spaced, but are set by air advection speed)
c		ldimtotal	  sum of dimlevels and layers (all the gas boxes in the whole model - this is used for the diffusivity input file)
c       dz            box size, meters
c		dzl			  tiny increment of depth used for locating temporal grid depths
c       timestep      the increment of time at each loop, s
c		samplestart   time in years after which time series begin to be collected
c
c     Constants:
c	u	      eddy diffusivity in convective zone, m2/s
c	ustart	      surface eddy diffusivity for convective zone, m2/s
c	uscale	      scale height of convective zone (exponential decay), m
c       Dm            molecular mass difference of gas pair, g/mol
c       g             acceleration of gravity, m/s^2
c       R             gas constant, J/mol K
c       omega         thermal diffusion sensitivity (permil per degree C)
c       pressure      atmospheric pressure, mbar
c	grav		Dm/1000*g*dz/R (has units of temperature K)
c	baro	molecular mass of air times g/RT (gradient of barometric pressure with depth)
c	aa		annual average temperature (C)
c	aaK		annual average temperature in Kelvin
c	newO2		initial value of atmospheric oxygen at surface
c	newO2date	date that applies to newO2 (keep it greater than length+1800 to avoid initiation
c			  of the reading of the surface concentration time series)
c	startseries	year CE in which program starts reading surface concentration time series (first
c			  date in time series file, equal to first value of newO2date if using timeseries)
c		
c   Selection of the appropriate time step:  (rule of thumb)
c       timestep should always be less than 1/3 of boxsize^2/diffusivity.
C       Approximate effective gas diffusivity = 2.3e-5 m2/s at surface
c
c 	HEAT COMPONENT OF GAS MODEL (courtesy of Richard Alley and Sigfus Johnsen)
c 		program calculates nonsteady temperature-depth profiles in firn
c
c 		program needs to know tlevels, the number of grid points (specify in declarations),
c  		then the distance step, accumulation rate, ice thickness and ice density,
c  		then the experiment length, time step, and number of output steps
c 		program produces temperature-depth curves for fixed bottom 
c  		and surface temperatures
c
c 		dkt(i)=temperature gradient of diffusivity (m2/s/C)
c 		dro(i)=depth gradient of density (kg/m3/m)
c 		ro(i)=density (kg/m3)
c 		sg(i)=load (kg/m2)
c 		t(i)=temperature (C)
c 		tnew(i)=temperatures at new time step
c 		z(i)=depth below surface (m)
c 		a=flow-law parameter (Pa-3s-2)
c 		bdot=accumulation rate of ice (m/s)
c 		cp=heat capacity (J/kg/C)
c 		dct=gradient of c with t (J/kg/C2)
c 		dkro=gradient of diffusivity with density (m2/s/kg/m3)
c 		dz=depth step (m)
c 		geotherm= initial temperature profile (C/m)
c 		g=gravitational acceleration (m/s2)
c 		tlevels=number of grid points in array
c 		jprint=test whether output this loop
c 		ntstep=number of time steps
c 		ntprint= number of times during experiment to print output
c 		qs=strain heating (C/s)
c 		rh=ice thickness (m)
c 		rhoi=ice density (kg/m3)
c 		rk=thermal diffusivity of firn (m2/s)
c 		spa=seconds per annum
c 		time=point in experiment (s)
c 		tmax=experiment length; read in as a and recalc as s
c 		tstep=time step; read in as steps/year and recalc as length in s
c 		vz=vertical velocity (m/s)
c 		w1,w2,w3,w4,w5,w6=computational intermediates
c
c
c*********************************************************
c   LOCKIN: A MODEL OF GAS FRACTIONATION IN THE FIRN LOCK-IN ZONE
c   Jeff Severinghaus, Scripps Inst Oceanography 9/20/2004
c
c	This model computes firn density with annual layers, 
c	summer layers having lower density than winter layers.  The lock-in zone begins
c	at a specified depth, the lock-in depth (LID).
c	Open porosity is decreased in the lock-in zone using an empirical fit 
c	to closed porosity data from Schwander et al. 1993.  This closing off
c	of open porosity results in size-fractionation of the remaining firn air.
c	This is for modeling small gases observed in the lock-in zone, such as Ne, Ar, and O2.
c
c     Variables (vector quantities positive downwards):
c       tl          		age of ice for lockin model, years
c       cl(i,j)          		concentration or "delta" of gas pair in box i for lockin model
c        [For example, cl(7,2) would be the dO2/N2 (gas 2) in the sixth box
c         down from the lock-in (the bottom firn box is lockinbox 1, with a depth
c         equal to the lock-in depth, lid.)]
c	rol				bulk density of firn in lockin
c	zl				depth below surface for lock-in model
c	sclosed(i)			closed porosity (from linear fit to Schwander 1993 data)	
c	dsclosed			temporary value of change in sclosed
c	q				temporary value of c(i)/1000 + 1
c	qtemp(i)			temporary value of q in the open porosity
c	porosity			total porosity (calculated from bulk density)
c	dporosity			annual change in total porosity due to densification
c	oporosityold			prior value of open porosity (in annual layer above)
c	oldporosity			prior value of total porosity (in annual layer above)
c
c     Parameters:
c	lid			lock-in depth (depth of the top of the lock-in zone)
c
c     Constants:
c	rhoice		real ice density from Bader
c
c initialize
      program firn
c      implicit real*8(a-h,o-z)
      implicit none
c These parameters determine the number of grid cells in the model.  The number of cells in the model is lt,
c	 and the thickness of the model is lt*dz - 1, where dz is the cell size.
c	 Extra cells are needed for computational reasons.  For example, ldimlevels should be h + 1, where
c	 h/dz is the thickness of the open firn column (and also the lock-in depth, or LID), and dz is the
c	 cell size.  For example, if the firn is 100 m thick, and the cells are 1 m in size, then dimlevels
c	 should be 101.  The thickness of the total firn column H must be <= (lt-1)*dz, so lt >= H/dz + 1. 
c	 Finding total firn thickness takes a little trial-and-error, because it depends on the closed porosity
c	 parameterization (by definition, total firn thickness is the depth where open porosity goes to zero).
c	 The number of separate layers in the lock-in zone is layers.
c	 ldimtotal is the sum of layers and ldimlevels - and is used for the empirical diffusivity file.
c --- improvement of the temperature model ----- 
c I did not change the equations, but separated it out in several depth intervals with varying resolutions:
c first step: same resolution as the gas, until the bottom of the gas model depth
c     for cod=85m and dz-0.2m, we need lclevel=426
c second step: 1m resolution until 500m, it is for the borehole temperature profile:
c     so 500-80 = 415 new boxes, so level1=lclevel+420=841
c third step: 25m resolution until the bottom. It's not a good model for that,
c     because of our vertical velocity parameterisation is not really valid until the bottom
c     for rh-3400m, 3400-500=2900, we need 116 depth, so lt=level1+116= 957
c define:
c     dz1= depth interval for the second step
c     dz2 =  depth interval for the third step
c     variabletemp = a boolean which decides whether we want to include temperature changes
c
c ----------------------------------------------------------------------
      integer ng, ldimlevels, lt,  layers, ldimtotal, level1,lclevel
      PARAMETER (ng = 8)
c      PARAMETER (ldimlevels = 336,lclevel=426, level1=841, lt = 957)
c      PARAMETER (ldimlevels = 68,lclevel=85, level1=500, lt = 616)
      PARAMETER (ldimlevels = 129,lclevel=171, level1=586, lt = 666)
c      PARAMETER (ldimlevels = 340,lclevel=441, level1=841, lt = 957)
      PARAMETER (layers = 5000, ldimtotal = 5340)
      real dkt(lt), dro(lt), ro(lt), sg(lt) 
      real t(lt), tnew(lt), z(lt), dvbubble(lt)      
      real s(lt), sc(lt),openpor(lt),w(lt),vb(lt),u(lt), du(lt)
      real deriv(ng,lt),c(ng,lt),cl(ng,layers),derivlock(ng,layers)
      real D(10),O2(ng),O2new(ng),O2date(ng),O2newdate(ng), Dm(ng)
      real buffer(ng),d1(ng,lt),d2(ng,lt),d0(ng,lt), decay(ng)
      real grav(ng,ldimlevels),dgravT(ng, ldimlevels)
      real dthermT(ng,ldimlevels), therm(ng,ldimlevels)
c      real year(layers)
      real depth(layers),wl(layers),pormidl(layers)
      real rolayer(layers),sclosed(layers),stotal(layers)
      real dvbubblel(layers),vbl(layers),diffmidl(layers)
      real thick(layers),sopen(layers), diff(ldimtotal)
      real ddiff(ldimtotal), dw(ldimtotal)
      real dporol(layers), ddiffl(layers), gradl(layers)
      real d2l(ng,layers), d1l(ng,layers), d0l(ng,layers)
      real newTsurfdate,newTsurf, startTseries, Tsurfdate, Tsurf
c      integer switch,switch1,
      integer dimlevels,dimtotal,tlevels, idl
      real*8 Fin, Fout
      character*160 fin0, fin33, fin34, fin5
      integer nfin, i, ii, jt, levels, indname

      real samplestart, aircontent, conveczone, bundle
      integer n, nlockstep
      real*8 dz, dzl, dz1, dz2
      real bdot, rh, aa, geotherm
      real length,tstep,ntprint,pressure,ustart,uscale
      real critdepth1, first0, first1, first2, first3
      real critdepth2, second0, second1, second2
c      real dummy, dummy1, dummy2
      real rhoice, spa, tmax
      integer klockstep, kbundle, ko, ntstep, kprint
      real g, R, dct, timestep, aaK, tmid, baro, pressurelock
      real   omega1(ng), omega2(ng), rhoclose, v
      real*8  dporo(lt)
      real airage, zl, wair
      integer j, jj, imax, jlock
      real*8 dlockz, dlockporin, dlockdiffin
      real time, ta, day
      real dt, d2t, rr, rk, dkro, cp, qs, vz, a1, a2, a3,a
c      double precision  gradout, gradin
      real grad(ldimlevels),  diff0(ldimtotal)
      real al, amp, jbb, pi, season, co2diff0
      integer bb, jdmax, jbmax,jbmin, jdmin
c integer nz, bmax
      character*160 ofolder, DiffFile, fileiniT, fThist
      integer charnum
c      real zi(lt),zdiff(ldimtotal)
      real zdmiddle(layers), diffmidl0(layers)
      real zdiff1, zdiff2, diff1, diff2, tp1, tp2, ztp1, ztp2
      real tp(lt)
      logical withtemp, iso, cprint
      integer ntemp, ktemp, tst
      real starttime, at, tendc
      real csave(ng,ldimlevels), clsave(ng, layers)
c --------------------------------------------------------------------

c
c	PARAMETERS TO BE SET
c       ----------------
c	lock-in depth (LID) in meters:
c	    lid = 67.
c max depth of the gas model (a bit more than cod
c            cod=79.

c     ratio of eddy diff to total diff in lock-in 0-1: .11
            al=0
c time in year A.D. when the simulation starts, supposes that concentration and
c  temperature time series are in year A.D.
       starttime=1800

c            print*, 'jmin=', jdmin, jdmax
c            print*, 'jbmin=', jbmin, jbmax
c	atmospheric gas mixing ratio or delta unit at surface used to force model
	   do 10 n = 1, ng
        O2new(n) = 1 
	O2(n) = 1
c	date to start reading surface gas history (10000 to skip history)
        O2date(n) = starttime
c     set decay to zero by default, and also grav
        decay(n)=0
        grav(n,1)=0
 10    continue
c if you want the temperature simulation either .true. or .false.
       withtemp=.true.
c if you want to update concentrations, set iso to true, if you are doing d15N, etc, set it to false
       iso=.true.
c enter the time when you want the concentrations to be saved (in run time)
c july 15 2008=day 196
       tendc=208.+196./365.
c logic input for the concentrations, initialised as .true.
       cprint=.true.
c	date to start reading surface temperature history (10000 to skip history)
	    newTsurfdate = 10000.
c	date to start sampling model to produce time series (10000 to skip history)
	    samplestart = 0.
c	total air content in m3STP/kg of ice. This quantity is produced by the model's
c	bubble closure parameterization but is entered here for the purpose of calculating
c	the advection velocity (if there is a disagreement then large spurious advection results).
c		aircontent = 0.000096037
c                aircontent = 8.6126d-5 ! for WAIS Divide
		aircontent = 9.4062E-005
c approximate depth of convective zone for model initialization
		conveczone = 4.
c number of times per year to calculate eddy diffusion in lock-in zone
		nlockstep = 1000
c number of years represented by each layer within the lock-in zone
		bundle = 0.5
c       ----------------
c for the temperature model:
            dz1=1.
            dz2=25.
c number of times per year to update the temperature model:
c need at least 3*diff(surface)/dz^2 steps per year, diff(surface)~1.76e-6 m2/s,
c  for 0.2m, we need 4166/year
c            ntemp=6000
c            ntemp=6000
c initialize site-specific parameters
      open (26,file='inputfiles/inputdata_NEEM_16.txt'
     $          ,status='old')
      read (26,*) dz, dzl, bdot, rh, aa, geotherm
      read (26,*) length,tstep,ntprint,pressure,ustart,uscale
c      read (26,*) critdepth1, first0, first1, first2
      read (26,*) critdepth1, first0, first1, first2, first3
      read (26,*) critdepth2, second0, second1, second2
      read(26,'(A80)') ofolder
      read(26,'(A80)') fileiniT
      read(26,'(A80)') DiffFile
      read(26,'(A80)') fThist
      read(26,*) jdmin, jdmax, amp, jbb
      close(26)
	print*, DiffFile

c     total depth in meters:
c            jdmax=67
c            jdmin=0
c     amplitude of diffusivity variations: 1e-6 - 1e-8
c            amp=1e-8
c     step of diffusivity variations:
c            jbb=0.5
            jbmin=int(jdmin/jbb)
            jbmax=int((jdmax)/jbb)
c for temperature variations
c        jbmin=jdmin
c        jbmax=jdmax

c --  not sure why we start here ---
c	  open (27,file=fThist,status='old')
c	  read (27,*) newTsurfdate,newTsurf
c	  startTseries = newTsurfdate
cc to disable the temperature history, set startTseries to 5000
cc		startTseries = 5000.
c 111         Tsurfdate=newTsurfdate
c             Tsurf=newTsurf
c	     read (27,*) newTsurfdate,newTsurf
c          if(starttime.ge.newTsurfdate) goto 111
 

		tlevels = lt
		dimtotal = ldimtotal
		dimlevels = ldimlevels
        levels = dimlevels - 1
        idl = tlevels - 1				

c	WHEN CHANGING GASES, YOU MUST CHANGE 5 DIFFERENT THINGS:
c	GAS-SPECIFIC CHANGE 1 (ATMOSPHERIC HISTORY - automatically zero if the next 14 lines are commented out)
c$$$	  open (41,file='SCENARIO_NEEM08_co2.txt',status='old')     
c$$$	  open (42,file='SCENARIO_NEEM08_ch4.txt',status='old')      
c$$$	  open (43,file='SCENARIO_NEEM08_ch3ccl3.txt',status='old')
c$$$	  open (44,file='SCENARIO_NEEM08_14co2.txt',status='old')
c$$$	  open (45,file='SCENARIO_NEEM08_hfc134a.txt',status='old')
c$$$	  open (46,file='SCENARIO_NEEM08_cfc113.txt',status='old')
c$$$	  open (47,file='SCENARIO_NEEM08_sf6.txt',status='old')
c$$$	  open (48,file='SCENARIO_NEEM08_cfc11.txt',status='old')
c$$$	  open (49,file='SCENARIO_NEEM08_cfc12.txt',status='old')
c$$$	  open (50,file='SCENARIO_NEEM08_ccl4.txt',status='old')

c free air diffusivity of CO2 (Matsunaga 1998-2009, Buizert 2011)
c$$$c 1.2059e-05 at -29C
c$$$        co2diff0=5.75d-10*(273.15-aa)**1.81
c$$$c      open (29,file='NEEMdiff_SIO_20cm_halfyr.txt',status='old')
c$$$        open(29, file=DiffFile,status='old' )
c$$$c        open(29, file='NEEMdiff_li6s_snr10tau10.txt')
c$$$c        open(29, file='diff_test03.txt')
c$$$      do 20 i = 1, dimtotal
c$$$c      read (29,*,end=20) zdiff(i),dummy1,diff(i)
c$$$      read (29,*,end=20) zdiff(i),diff(i)
c$$$c     verify it is not larger than the free air diffusivity:
c$$$      if(diff(i).ge.co2diff0) then
c$$$         print*, 'diffusivity input is too large at z=', zdiff(i),
c$$$     $    'free air diffusivity is', co2diff0, 'diff(i)=', diff(i)
c$$$         diff(i)=co2diff0
c$$$       end if
c$$$c     store the reference diffusivity profile
c$$$      diff0(i)=diff(i)
c$$$
c$$$  20  continue
c$$$c verify that ustart is not too large:
c$$$       if(ustart.ge.co2diff0) then
c$$$         print*, 'eddy diffusivity is too large',
c$$$     $    'free air diffusivity is', co2diff0, 'ustart', ustart
c$$$        ustart=co2diff0
c$$$       end if

c calculate real ice density (Bader, in Schwander et al. 1997)
      rhoice = 916.5 - 0.14438*aa - 0.00015175*aa**2
      spa = 3600.*24.*365.25
c speed up model by skipping steps in lock-in zone
	  klockstep = tstep/nlockstep
	  kbundle = tstep/aint(1/bundle)
      tmax = length*spa
c speed up model by skipping steps in the temp model:
c      ktemp=tstep/ntemp
      ktemp=10
c skip steps in printing timeseries output
      kprint =tstep/ntprint     

c convert timestep from frequency to time units
      tst = tstep
      tstep = spa/tstep
      ko = 1
      g = 9.82
      R = 8.314
      sg(1) = 0.
      dct = 7.755
      ntstep = tmax/tstep
c skip steps in printing timeseries output
c      kprint = ntstep/(ntprint-1)     
      aaK = aa + 273.15
	  tmid = aaK
	  baro = 0.02895*g/R/aak
	  pi = 3.1415927

c          baro=0
c calculate density at which air content becomes fixed (Martinerie et al., 1994)	  
	  pressurelock = 1/(1/rhoice+0.000000695*aaK-0.000043)
	  
c get the folder name for output files
      CALL system('mkdir ' // ofolder )
c find length of folder string ignoring trailing blanks:
      DO 15, I = LEN(ofolder), 1, -1 
         IF(ofolder(I:I) .NE. ' ') GO TO 21 
15    CONTINUE 
21    charnum = I 
            fin0 = ofolder(1:charnum) //
     $         "/NEEMfirnoutput16_test1_g0_1000.33"
          nfin=charnum+31


c****************************************************
c		GAS PARAMETERS TO BE SET
c	
c   GAS-SPECIFIC CHANGE 2 (Radioactive decay)
c	   decay = 1/8267./spa

c	GAS-SPECIFIC CHANGE 3 (DIFFUSIVITY IN AIR - "D" - RELATIVE TO CO2)
c	NEEM-recommended diffusivities, D/DCO2 (mostly Matsunaga):c	CO2-air
			D(6) = 1.000000
c	CH4-air
			D(7) = 1.367
c	N2O-air
c			D(n) = 0.981
c	SF6-air
c			D(7) = 0.554
c	H2-air
c			D(n) = 4.694
c	CFC-11
c			D(8) = 0.525
c	CFC-12
c			D(9) = 0.596
c	CFC-113
c			D(6) = 0.453
c	HFC-134a
c			D(5) = 0.630
c	CH3Br
c			D(n) = 0.736
c	CF4
c			D(n) = 0.823
c	CH3CCl3
			D(8) = 0.485
c	CCl4
c			D(10) = 0.470
c	CH3Cl (from Chen and Othmer)
c			D(n) = 0.789
c	CO	(from Chen and Othmer)
c			D(n) = 1.250
c	N2 (from Chen and Othmer)
c			D(1) = 1.275
c	O2 (from Chen and Othmer)
c			D(n) = 1.302
c	He	(from Chen and Othmer)
c			D(n) = 4.780
c	Ne (from Chen and Othmer)
c			D(n) = 2.140
c	Ar (from Chen and Othmer)
c			D(n) = 1.230
c	Kr (from Chen and Othmer)
c			D(n) = 0.962
c	Xe (from Chen and Othmer)
c			D(n) = 0.835
c	HD
c			D(n) = 3.895386
c	HH
c			D(n) = 4.693994
c	13CH4
c			D(n) = 1.340806
c	12CH4
c			D(n) = 1.366969
c	CDH3 (natural isotopic abundance 13C)
c			D(n) = 1.340469
c	CH4 (natural isotopic abundance 13C)
c			D(2) = 1.366676
c	14CO2 (natural isotopic abundance 17O and 18O)
c			D(4) = 0.991368
c	CO2 (natural isotopic abundance)
c			D(n) = 1.000000
c	13CO2 (natural isotopic abundance 17O and 18O)
c			D(n) = 0.995613
c	12CO2 (natural isotopic abundance 17O and 18O)
c			D(n) = 1.000048
c	12C18OO
c			D(n) = 0.991401
c	12C16OO
c			D(n) = 1.000089
c	12C18O
c			D(n) = 1.228754
c	12C17O
c			D(n) = 1.239117
c	12C16O
c			D(n) = 1.250172
c	15NN
			D(1) = 1.263893
c	14NN
c			D(1) = 1.275084
c                        D(2)=D(1)
c	18OO
c			D(n) = 1.283719
c	17OO
c			D(n) = 1.292637
c	16OO
c			D(n) = 1.302087
c	22Ne
c			D(n) = 2.087122
c	20Ne
c			D(n) = 2.145608
c	40Ar
			D(2) = 1.229952
c	38Ar
c			D(n) = 1.243488
c	36Ar
			D(3) = 1.258324
c	86Kr
			D(4) = 0.958741
c	84Kr
c			D(n) = 0.961616
c	82Kr
			D(5) = 0.964621
c	136Xe
c			D(n) = 0.832366
c	132Xe
c			D(n) = 0.834581
c	129Xe
c			D(n) = 0.836327
c		
		
c	GAS-SPECIFIC CHANGE 4 (MASS DIFFERENCE Dm, IF GRAVITATIONAL SETTLING IS INCLUDED)
c      For CO2, CH4, 134a, etc., Dm should be zero:
c		 Dm = 0.		 
c		 grav = Dm/1000*g*dz/R
		 Dm(1) = 1.
		 Dm(2) = 40.-28.
		 Dm(3) = 36.-28.
		 Dm(4) = 86.-28.
		 Dm(5) = 82.-28.
                 Dm(6) = 0.
                 Dm(7)=0.
                 Dm(8)=0.
c		 grav(1) = Dm(1)/1000*g/R/aaK
c                 print*, 'grav(1)=', grav(1), 'diff(1)*D(1)*grav(1)=',
c     $                diff(1)*D(1)*grav(1)
c		 grav(2) = Dm(2)/1000*g/R/aaK
               do 30 n=1,ng
c                  Dm(n)=0.
                  grav(n,1)=Dm(n)/1000.*g/R/aaK
c                  decay(n)=0.
c                  omega1(n)=0
c                  omega2(n)=0
 30            continue

c	GAS-SPECIFIC CHANGE 5 (THERMAL DIFFUSION SENSITIVITY - "omega1" and "omega2")
c
c		coefficients for quadratic T dependence of Omega-N15 (Omega = Omega1/T - Omega2/T^2)
		omega1(1) = 8.656
		omega2(1) = 1232.
c		coefficients for quadratic T dependence of Omega-Ar-40
c		omega1(2) = 26.08
c		omega2(2) = 3952.
c		coefficients for quadratic T dependence of Omega-O2/N2 
c		omega1 = 8.656*3.3
c		omega2 = 1232.*3.3
c		coefficients for quadratic T dependence of Omega-Ne/Ar 
c		omega1 = -8.656*48
c		omega2 = -1232.*48
c		coefficients for quadratic T dependence of Omega-Ne/N2 
c		omega1 = -8.656*13
c		omega2 = -1232.*13
c		coefficients for quadratic T dependence of Omega-Ar/N2 
		omega1(2) = 8.656*16
		omega2(2) = 1232.*16
c for 36/28 = 40/28 - 40/36
		omega1(3) = 8.656*16 -26.08
		omega2(3) = 1232.*16 -3952.
c		coefficients for quadratic T dependence of Omega-Kr/N2 
		omega1(5) = 8.656*32
		omega2(5) = 1232.*32
c for 86/28 = 86/82+82/28, not sure what Kr/N2 is
		omega1(4) = 8.656*32+5.05
		omega2(4) = 1232.*32+580.

c		coefficients for quadratic T dependence of Omega-Xe/N2 
c		omega1 = 8.656*35
c		omega2 = 1232.*35
c		coefficients for quadratic T dependence of Omega-Kr-86/82
c		omega1 = 5.05
c		omega2 = 580.
c		coefficients for quadratic T dependence of Omega-Xe-136/129
c		omega1 = 11.07
c		omega2 = 2000.
c
c     For CO2, CH4, 134a, etc., omega1 and omega2 should be zero:
c		omega1 = 0.
c                omega2 = 0.
                omega1(6)=0.
                omega2(6)=0.
                omega1(7)=0.
                omega2(7)=0.
                omega1(8)=0.
                omega2(8)=0.
c
c ----------------------------------------------------------------------
c *************** BEGIN OPEN FIRN COLUMN INITIALIZATION *************************
		vb(1) = 0.

c calculate depths and densities in open firn column using quadratic fit to data:
       rhoclose = second0 + second1*critdepth2 + second2*critdepth2**2	

c free air diffusivity of CO2 (Matsunaga 1998-2009, Buizert 2011)
c 1.2059e-05 at -29C
        co2diff0=5.75d-10*(273.15-aa)**1.81
c open diffusivity file:
        open(29, file=DiffFile,status='old' )
      read (29,*) zdiff1,diff1
      read (29,*) zdiff2,diff2


cc verify that ustart is not too large:
c       if(ustart.ge.co2diff0) then
c         print*, 'eddy diffusivity is too large',
c     $    'free air diffusivity is', co2diff0, 'ustart', ustart
c        ustart=co2diff0
c       end if

c open temperature profile:
      open(28, file = fileiniT, status = 'old')
      read(28,*) ztp1, tp1
      read(28,*) ztp2, tp2
c ------- initialisation loop --------------------
c      do 40 i=1,tlevels
       do 40 i=1,lclevel

c ------------------------------------------------
      z(i) = dz*(i-1)

c calculate diffusivity for the diffusive zone only:
      if(i.le.dimlevels) then 
 292     if(z(i).gt.zdiff2) then
            zdiff1=zdiff2
            diff1=diff2
            read (29,*) zdiff2,diff2
         end if
         if(z(i).gt.zdiff2) goto 292
         diff(i)=diff1+(diff2-diff1)/(zdiff2-zdiff1)*(z(i)-zdiff1)
c     verify it is not larger than the free air diffusivity:
         if(diff(i).ge.co2diff0) then
            print*, 'diffusivity input is too large at z=', z(i),
     $       'free air diffusivity is', co2diff0, 'diff(i)=', diff(i)
            diff(i)=co2diff0
         end if
         if(diff(i).le.0) then
            print*, 'error: negative diffusivity at z=', z(i),
     $           'diff(i)=', diff(i)
            diff(i)=0
            end if
c     store the reference diffusivity profile
         diff0(i)=diff(i)
      end if

c 
c calculate density
cc      if (z(i).lt.critdepth1) ro(i)=first0+ first1*z(i)+ 
cc     $      first2*exp(first3*(critdepth1-z(i)))	  
c	  if(z(i).lt.critdepth1)ro(i)=first0+ first1*z(i)+ first2*z(i)**2
c      if(z(i).ge.critdepth1)ro(i)=second0+second1*z(i)+second2*z(i)**2		
c      if (z(i).ge.critdepth2) ro(i) = rhoice-(rhoice-rhoclose)
c     $  *exp(-(z(i)-critdepth2)/(rhoice-rhoclose)*(second1+2*second2
c     $     *critdepth2))
c calculate depths and densities in open firn column using quadratic fit to data:
      if (z(i).lt.critdepth1) 
     $  ro(i)=first0+first1*z(i)+first2*exp(first3*(critdepth1-z(i)))
      if (z(i).ge.critdepth1) ro(i)=second0+second1*z(i)+second2*z(i)**2		
      if (z(i).ge.critdepth2) ro(i) = rhoice-(rhoice-rhoclose)
     $  *exp(-(z(i)-critdepth2)/(rhoice-rhoclose)*(second1+2*second2
     $     *critdepth2))
  					
c --------------------------------------------------------------------

c calculate total porosity in whole firn
          s(i) = 1-ro(i)/rhoice
			
c calculate closed porosity in whole firn (Goujon et al., 2003; Emmanuel Winant in prep)
	sc(i) = 0.37*s(i)*(s(i)/(1-831.2/rhoice))**(-7.6)
	if (sc(i).gt.s(i)) sc(i) = s(i)
		
c calculate open porosity in whole firn			
        openpor(i) = s(i) - sc(i)
	if (openpor(i).lt.0.) openpor(i) = 0.

c calculate volume of bubbles (corrected for compression) 
	if (i.eq.1) go to 39
	dvbubble(i) = 2*(sc(i) - sc(i-1)*(s(i)/s(i-1)
     $	+ ro(i)/ro(i-1)-1))/(s(i)/s(i-1)+1)
	vb(i) = dvbubble(i) + vb(i-1)*ro(i)/ro(i-1)
		
c calculate vertical snow (matrix) velocity
 39	v = bdot*rhoice/ro(i)

c calculate vertical air velocity w at box nodes:
 	w(i) = bdot*rhoice/openpor(i)*(aircontent*aaK/273.15*1013.25
     $   /pressure - vb(i)/ro(i))		
	 if (ro(i).gt.(pressurelock-10.)) then
            w(i) = (v-w(i))*(ro(i)-pressurelock+10.)/20.+w(i)
	 end if	
	 if (ro(i).gt.(pressurelock+10.)) w(i) = v
c WARNING: this air velocity isn't exactly the same as the air velocity that is calculated from air content; it tapers to the firn velocity.  It is just the air velocity used for tracer advection.
c It is important to use a value of aircontent that is consistent with the bubble close-off parameterization, or else spurious air velocity can result near the full close-off.
c -----------------------------------------------------------------
c ----- For convective zone (calculated at nodes): --------
    	u(i) = ustart*exp((1.-i)*dz/uscale)
        du(i)=-u(i)/uscale
c      2. other possibility: an arc tan
c        u(i)=ustart*(1/2.-(ATAN(((i-1.)*dz-uscale)/0.1))/3.1415927)
c        du(i)= -ustart/3.1415927/(1+(((i-1.)*dz-uscale)/0.1)**2)/0.1
c        u(i)=ustart*(1/2.-(ATAN(((i-1.)*dz/uscale)))/3.1415927)
c       du(i)= -ustart/3.1415927/(1+(((i-1.)*dz/uscale))**2)
c     3. other posibility piecewise linear:
c       if(i*dz.le.uscale) then
c           u(i)=ustart*(1.-(i-1.)*dz/uscale)
c           du(i) =-ustart/uscale
c          u(i)=ustart
c          du(i)=0
c           else 
c           u(i)=0
c           du(i) =0
c           end if
c      if(i*dz.eq.uscale) du(i)=ustart/uscale
c     4. Or  try a function of openpor: ---
c                 u(i)=-ustart*openpor(i)/openpor(1)
c                 du(i)=-ustart*(openpor(i+1)-openpor(i-1))/2/dz
c --------------------------------------------------------------------
	  		
c initialize variables useful in the temperature model
c	t(i) = aa + geotherm*dz*(i-1)
	if (i.gt.1) sg(i)=sg(i-1)+ro(i)*dz
	dkt(i)=1.084d-16*ro(i)**3-1.88d-13*ro(i)**2
     $  +8.291d-11*ro(i)-1.569d-8  
	
  40	continue

c -- continue initialisation for the temperature model
c --- stage 2, with 1m spacing:
        do 41 i=lclevel+1, level1
           z(i) = z(lclevel)+ dz1*(i-lclevel)
            ro(i) = rhoice-(rhoice-rhoclose)
     $  *exp(-(z(i)-critdepth2)/(rhoice-rhoclose)*(second1+2*second2
     $     *critdepth2))
	sg(i)=sg(i-1)+ro(i)*dz1
	dkt(i)=1.084d-16*ro(i)**3-1.88d-13*ro(i)**2
     $  +8.291d-11*ro(i)-1.569d-8  
 41        continue
c --- stage 3, with 25m spacing:
        do 42 i=level1+1, lt
           z(i) = z(level1)+ dz2*(i-level1)
            ro(i) = rhoice-(rhoice-rhoclose)
     $  *exp(-(z(i)-critdepth2)/(rhoice-rhoclose)*(second1+2*second2
     $     *critdepth2))
	sg(i)=sg(i-1)+ro(i)*dz2
	dkt(i)=1.084d-16*ro(i)**3-1.88d-13*ro(i)**2
     $  +8.291d-11*ro(i)-1.569d-8  
 42   continue

           print*, 'dz=', dz, 'until', z(lclevel),', then dz=', dz1, 
     $          'until', z(level1), 'then dz=', dz2,  'until', z(lt)

c --- calculate dro:
           dro(1)=(ro(2)-ro(1))/(z(2)-z(1))
           dro(lt)=(ro(lt)-ro(lt-1))/(z(lt)-z(lt-1))
           do 43, i=2,lt-1
              dro(i)=(ro(i+1)-ro(i-1))/(z(i+1)-z(i-1))
 43           continue

c ----  initialize temperatures with stationary profile ------
c	t(i) = aa + geotherm*dz*(i-1)
              do 44, i=1,lt
      if(z(i).gt.ztp2) then
         ztp1=ztp2
         tp1=tp2
         read (28,*) ztp2,tp2
      end if
      tp(i)=tp1+(tp2-tp1)/(ztp2-ztp1)*(z(i)-ztp1)
      t(i)=tp(i)
 44   continue
        close(28)
c        close (29)
        print*, 'done with diffusive firn'
	  open (49,file=ofolder(1:charnum) 
     $       //'Densityetc.d',status='replace')
       write(49,'(8f12.7)')z(1),ro(1),openpor(1),s(1),sc(1),0.,0.,0.
c ****************************************************
c calculate density, porosity, diffusivity gradients
        do 50 i=2,dimlevels-1

c        do 50 i=2,tlevels-1
c       		dro(i)=(ro(i+1)-ro(i-1))/(2.*dz)
c     change ddiff into ddiff(i) - Anais
c	ddiff(i) = (diff(i+1)-diff(i-1))/2/dz
	dporo(i) = (openpor(i+1)-openpor(i-1))/2/dz/openpor(i)
c add dw(i) first derivative of vertical velocity -Anais
                dw(i)=(w(i+1)-w(i-1))/2/dz
c			if (i.eq.dimlevels) dlockpor = dporo  
c--- no need for this, move it to the ng loop
c		 do 48 n = 1,ng	
c		d1(n,i) = ddiff*D(n)-u(i)/uscale + diff(i)*D(n)*
c     $       (dporo + baro) - w(i)

c ---anais full ---
c$$$        d2(n,i) = diff(i)*D(n) + u(i)
c$$$	d1(n,i) = (dporo(i)+baro)*(diff(i)*D(n)+u(i))+ddiff(i)*D(n)
c$$$     $          +du(i) -diff(i)*D(n)*grav(n)-w(i)
c$$$        d0(n,i)= (dporo(i)+baro)*(-diff(i)*D(n)*grav(n))
c$$$     $      -ddiff(i)*D(n)*grav(n) - decay(n) 
	 	
c --------- Jeff's expressions ----------------------------------------		
c		d1(n,i) = u(i)*baro - diff(i)*D(n)*(grav(n)-omega(n)*dt)
c     $	 +	  ddiff(i)*D(n) - u(i)/uscale + (diff(i)*D(n)+u(i))
c     $   *(dporo(i) + baro) - w(i)
c
c	    d0(n,i) = diff(i)*D(n)*(grav(n)*dt/aaK+omega(n)*d2t
c     $   + domega*dt*dt) - ddiff(i)*D(n)*(grav(n)-omega(n)*dt)
c     $	 - baro*u(i)*(1/uscale+dt/aak) 
c     $   +(u(i)*baro-diff(i)*D(n)*(grav(n)-omega(n)*dt))
c     $	 *(dporo(i)+baro) - decay(n)
c --------------------------------------------------------------------
c   48    continue
c---- end no need for this
          write(49, '(8g12.5)')z(i), ro(i), openpor(i), s(i), 
     $               sc(i),dporo(i),w(i),dw(i)

  50   continue
        print*, 'done with diffusive firn derivatives'
        close(49)
c ****************************************************
c *************** END OPEN FIRN COLUMN INITIALIZATION *************************

c       --------------------------------------------------------

c ************ INITIALIZE LOCK-IN ZONE ****************************************
c calculate depths of grid points in the lockin zone for air advection scheme

c     buffer = 0.
 		j = 1
		jj=dimlevels
	    airage = 0.	
		imax = aint(dz*(lt-dimlevels)/dzl)
		wl(1) = w(dimlevels)
		depth(1) = z(dimlevels)
		rolayer(1) = ro(dimlevels)
		stotal(1) = s(dimlevels)
		sclosed(1) = sc(dimlevels)
		sopen(1) = openpor(dimlevels)
c		diffmidl(1) = (diff(dimlevels)+diff(dimlevels+1))/2
                do 84 n=1,ng
		derivlock(n,1) = 0.
 84             continue
c begin loop to assign depths to lock-in layers
		do 81 i = 1, imax
	    zl = dzl*(i-1) + z(dimlevels)

c calculate lock-in "air advection age" by integrating air velocity
		if (zl-z(jj).gt.dz) jj = jj + 1
                if(jj.ge.lclevel) then
                   print*, 'lclevel too small'
                   call abort
                end if
		wair = (w(jj+1)-w(jj))/dz*(zl-z(jj))+w(jj)
		airage = airage + dzl/wair/spa
	
c identify layer with constant air age spacing equal to bundle
		if (airage.ge.bundle) then
			airage = 0.
			j = j+1
			
c record the air velocity at the layer
			wl(j) = wair
				
c find the depth of the firn layer, for each bundle	
			depth(j)=zl

c find the thickness of the layer
			thick(j-1) = depth(j) - depth(j-1)
			
c 	find the bulk density of the layer
		rolayer(j) = second0 + second1*zl + second2*zl**2		

c 	calculate the total porosity:		     
		stotal(j) = 1 - rolayer(j)/rhoice

c calculate closed porosity (Goujon et al., 2003; Emmanuel Winant, in prep.)
	sclosed(j)=0.37*stotal(j)*(stotal(j)/(1-831.2/rhoice))**(-7.6)
		sopen(j) = stotal(j) - sclosed(j)

c calculate diffusivities at the midpoints
             if(j.eq.2) then
                zdmiddle(1)=z(dimlevels)+thick(j-1)/2
 290            if(zdmiddle(1).gt.zdiff2) then
                   zdiff1=zdiff2
                   diff1=diff2
                   read (29,*) zdiff2,diff2
                end if
                 if(zdmiddle(1).gt.zdiff2) goto 290
c$$$                if(zmiddle(1).le.zdiff(1)) then
c$$$                   backspace 29
c$$$                   backspace 29
c$$$                   read (29,*) zdiff1,diff1
c$$$                   read (29,*) zdiff2,diff2
c$$$                   end if
                diffmidl(1)=diff1+(diff2-diff1)/(zdiff2-zdiff1)
     $               *(zdmiddle(1)-zdiff1)
         diffmidl0(1)=diffmidl(1)
         if(diffmidl(1).le.0) then
            print*, 'error: negative diffusivity at z=', zdmiddle(1),
     $           'diff(i)=',diffmidl(1)
           diffmidl(1)=0
         end if

             end if 

                zdmiddle(j)=depth(j-1)+thick(j-1)/2
 291     if(zdmiddle(j).gt.zdiff2) then
            zdiff1=zdiff2
            diff1=diff2
            read (29,*) zdiff2,diff2
         end if
         if(zdmiddle(j).gt.zdiff2) goto 291

         diffmidl(j)=diff1+(diff2-diff1)/(zdiff2-zdiff1)
     $        *(zdmiddle(j)-zdiff1)
c store this value for later:
         diffmidl0(j)=diffmidl(j)
c		diffmidl(j) = (diff(levels+j)+diff(levels+j+1))/2
		
c --------------------------------------------------------------------
         if(diffmidl(j).le.0) then
            print*, 'error: negative diffusivity at z=', zdmiddle(j),
     $           'diff(i)=',diffmidl(j)
           diffmidl(j)=0
           end if
      if (sopen(j).le.0.) go to 85
				
      end if
				
   81 continue
          print*, 'done with lid', zdiff1,diff1,zdiff2,diff2

c calculate composition of gases in remaining open pores 
   85	 jlock = j
        sopen(j) = 0.
	pormidl(1) = (sopen(1)+sopen(2))/2
	vbl(1) = vb(dimlevels)
        print*, 'done with lid', zdiff1,diff1,zdiff2,diff2

        print*,'nb of steps in lockin:jlock=', jlock
        print*, 'depth',depth(jlock)
c initialize firn air composition in open pores
c		do 89 n = 1, ng
c		    do 87 i = 1, j-1
c              cl(n,i) = c(n,dimlevels)	
c 87			continue 
c 89	    continue
           close(29)
		do 91 i = 2, jlock-1
		
c calculate each layer's thickness, and the volume increment of bubbles enclosed each year
     
      dvbubblel(i)=2*(sclosed(i)-sclosed(i-1)*(stotal(i)/stotal(i-1)
     $		+ rolayer(i)/rolayer(i-1)-1))/(stotal(i)/stotal(i-1)+1)
	vbl(i) = dvbubblel(i) + vbl(i-1)*rolayer(i)/rolayer(i-1)

c calculate the midpoint porosity
	   pormidl(i) = (sopen(i)+sopen(i+1))/2
c    --- calculate the intermediate for time dependant simulation:
c     review: thick(i) is the distance between node i and i-1.
           dporol(i)=(sopen(i+1)-sopen(i-1))/(thick(i)+thick(i+1))
     $          /sopen(i)
c           ddiffl(i)=(diff(levels+i+1)-diff(levels+i-1))
c     $          /(thick(i)+thick(i+1))
           ddiffl(i)=(diffmidl(i)-diffmidl(i-1))*2/(thick(i)-thick(i-1))
c$$$           do 92 n=1,ng
c$$$	d2l(n,i) = diff(levels+i)*(al+(1-al)*D(n))
c$$$	d1l(n,i)=(dporol(i)+baro)*(diff(levels+i)*(al+(1-al)*D(n)))+
c$$$     $       ddiffl(i)*(al+(1-al)*D(n))
c$$$        d0l(n,i)=  - decay(n) 
c$$$ 92     continue
   91 continue

c      jlock = j
        print*, 'vbl(jlock-1)=', vbl(jlock-1), 
     $             'aircontent=',
     $  vbl(jlock-1)/rolayer(jlock-1)*pressure/1013.25*273.15/aak
c     $ aircontent*aak/273.15*1013.25/pressure*rolayer(jlock-1)
        print*, 'vb(dimlevels)=', vb(dimlevels) 

c calculate thickness of lock-in box, diffusivities and porosities
	  dlockz = dz/2 + thick(1)/2
	  dlockporin = (openpor(levels)+openpor(dimlevels))/2
	  dlockdiffin = (diff(levels)+diff(dimlevels))/2
          print*, 'depth step in the lock in zone=', dlockz
c	  insure that the bottom of the model is impermeable
c	  diff(jlock-1+levels)=0.
          diffmidl(jlock-1)=0.
c     --- adjust the derivatives of poro and diff at dimlevels:
        dporo(dimlevels)= (pormidl(1)-dlockporin)/sopen(1)/dlockz
        ddiff(dimlevels)= (diffmidl(1)-dlockdiffin)/dlockz


c		
c ****************** END OF LOCK-IN INITIALIZATION ******************************		
c -------------------------------------------------------------------
c -------------------------------------------------------------------
c need to start loop over several diffusivity profiles
c       bb=0
c        do 1200 bb=jbmax,jbmin, -1
        do 1200 bb=jbmin,jbmax
c -------------------------------------------------------------------
c ---  update fin with bb:
        indname = 1000	  
        write(fin0(nfin-3:nfin), '(I4)') indname
       write(fin0(nfin-5:nfin-5), '(I1)') 0

        open(36, file=fin0(1:nfin) // ".36"
     $       , status='replace')
        open(37, file=fin0(1:nfin) // ".37"
     $       , status='replace')
        open(34, file=fin0(1:nfin) // ".34"
     $       , status='replace')
        n=1
 910     format(2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5,
     $        2x, g12.5,2x, g12.5,
     $        2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5)
 911     format(2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5,
     $        2x, g12.5,2x, g12.5,
     $        2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5,2x, g12.5)


c --- update  and save the diffusivity profile
           do 71 i=1, dimlevels
           diff(i)=diff0(i)
              if(z(i).ge.bb*jbb.and.z(i).lt.(bb+1)*jbb) then
           diff(i)=diff0(i)+(z(i)-bb*jbb)/jbb*amp*
     $                diff0(int((bb+1)*jbb/dz))
           end if
             if(z(i).ge.(bb+1)*jbb.and.
     $          z(i).lt.(bb+2)*jbb) then
           diff(i)=diff0(i)+((bb+2)*jbb-z(i))/jbb*amp*
     $                diff0(int((bb+1)*jbb/dz))

           end if
c save diffusivity profile
           write(34, FMT='(1x,f12.5,2x,E12.5,1x,E12.5,1x)')
     $          z(i),diff0(i),diff(i)
 71        continue

           do 72 j=1, jlock
           diffmidl(j)=diffmidl0(j)
              if(zdmiddle(j).ge.bb*jbb.and.
     $          zdmiddle(j).lt.(bb+1)*jbb) then
           diffmidl(j)=diffmidl0(j)+(zdmiddle(j)-bb*jbb)/jbb*amp
           end if
             if(zdmiddle(j).ge.(bb+1)*jbb.and.
     $          zdmiddle(j).lt.(bb+2)*jbb) then
           diffmidl(j)=diffmidl0(j)+((bb+2)*jbb-zdmiddle(j))/jbb*amp
           end if
c append to the diffusivity file
           write(34, FMT='(1x,f12.5,2x,E12.5,1x,E12.5,1x)')
     $          zdmiddle(j),diffmidl0(j),diffmidl(j)
 72        continue

        close(34)
c --------------------------------------------------------------------
c update cprint:
                cprint=.true.
c impermeability at the end of the domaine:
          diffmidl(jlock-1)=0.
          diffmidl(jlock)=0.

c Update diffusivity gradients:
        do 73 i=2,dimlevels-1
	ddiff(i) = (diff(i+1)-diff(i-1))/2/dz
 73     continue
c at the middle of dimlevel box:
	  dlockdiffin = (diff(levels)+diff(dimlevels))/2
c     --- adjust the derivatives at dimlevels:
        ddiff(dimlevels)= (diffmidl(1)-dlockdiffin)/dlockz
c -- in the lock-in zone
        do 74 i=2,jlock-1
        ddiffl(i)=(diffmidl(i)-diffmidl(i-1))*2/(thick(i)-thick(i-1))
 74      continue
c note: there is no ddiffl(1), it would be ddiff(dimlevels), 
c since it refers to the node dimlevels

c --- print all the terms in the equation, to see which one matters:

         do 75 i=1,dimlevels
         write(36,fmt=910) z(i),dporo(i),baro,diff(i), ddiff(i),
     $        ddiff(i)*D(n), du(i),-diff(i)*D(n)*grav(n,1),
     $        u(i)*baro, w(i),ddiff(i)*D(n)*grav(n,1),du(i)*baro

        write(37,fmt=911) z(i),openpor(i),dporo(i),ro(i), dro(i),
     $         diff(i),ddiff(i), s(i), sc(i), w(i), u(i), du(i)
 75      continue

         do 76 i=2,jlock-1
        write(36,fmt=910) depth(i),dporol(i),baro,diffmidl(i),ddiffl(i),
     $        ddiffl(i)*D(n), du(i),0,
     $        u(i)*baro, wl(i),ddiffl(i)*D(n)*grav(n,1),du(i)*baro

          write(37,fmt=911) depth(i),sopen(i),dporol(i),rolayer(i),
     $          0,diffmidl(i),ddiffl(i), stotal(i), sclosed(i),
     $        wl(i), al*diffmidl(i),al*ddiffl(i)

 76      continue
	close(36)
        close(37)

c -------------------------------------------------
c -------------------------------------------------
c grand loop through the different gases
		do 1100 n = 1, ng
c ------- derivatives initialisation-------------
       write(fin0(nfin-5:nfin-5), '(I1)') n
         open(38, file=fin0(1:nfin) // ".38"
     $       , status='replace')
               open(35, file=fin0(1:nfin) // ".35"
     $       , status='replace')
               open(39, file=fin0(1:nfin) // ".39"
     $       , status='replace')

c$$$        do 65 i=1,dimlevels-1
c$$$        d2(n,i) = diff(i)*D(n) + u(i)
c$$$	d1(n,i) = (dporo(i)+baro)*(diff(i)*D(n)+u(i))+ddiff(i)*D(n)
c$$$     $          +du(i) -diff(i)*D(n)*grav(n)-w(i)
c$$$        d0(n,i)= (dporo(i)+baro)*(-diff(i)*D(n)*grav(n))
c$$$     $      -ddiff(i)*D(n)*grav(n) - decay(n) 
c$$$       if(isnan(d2(n,i))) print*, 'NaN alert 918 d2', n,i
c$$$       if(isnan(d1(n,i))) print*, 'NaN alert 918 d1', n,i
c$$$       if(isnan(d0(n,i))) print*, 'NaN alert 918 d0', n,i
c$$$       write(38, '(4f12.7)')z(i), d2(n,i), d1(n,i), d0(n,i)
c$$$ 65     continue
c$$$        d2(n,dimlevels) = diff(dimlevels)*D(n)
c$$$	d1(n,dimlevels)=(dporo(dimlevels)+baro)*(diff(dimlevels)*D(n))+
c$$$     $       ddiff(dimlevels)*D(n)
c$$$        d0(n,dimlevels)=  - decay(n) 
c$$$        i=dimlevels
c$$$       write(38, '(4f12.7)')z(i), d2(n,i), d1(n,i), d0(n,i)
c$$$
c$$$        do 66 i=2,jlock-1
c$$$	d2l(n,i) = diffmidl(i)*(al+(1-al)*D(n))
c$$$	d1l(n,i)=(dporol(i)+baro)*(diffmidl(i)*(al+(1-al)*D(n)))+
c$$$     $       ddiffl(i)*(al+(1-al)*D(n))
c$$$        d0l(n,i)=  - decay(n) 
c$$$c	d2l(n,i) = diffmidl(i)*D(n)
c$$$c	d1l(n,i)=(dporol(i)+baro)*(diffmidl(i)*D(n))+
c$$$c     $       ddiffl(i)*D(n)
c$$$c        d0l(n,i)=  - decay(n) 
c$$$       write(38,'(4f12.7)')depth(i),d2l(n,i), d1l(n,i), d0l(n,i)
c$$$
c$$$       if(isnan(d2l(n,i))) print*, 'NaN alert 918 d2', n,i
c$$$       if(isnan(d1l(n,i))) print*, 'NaN alert 918 d1', n,i
c$$$       if(isnan(d0l(n,i))) print*, 'NaN alert 918 d0', n,i
c$$$
c$$$ 66   continue
c$$$      close(38)
        do 65 i=1,dimlevels-1
        d2(n,i) = diff(i)*D(n) + u(i)
	d1(n,i) = (dporo(i)+baro)*(diff(i)*D(n)+u(i))+ddiff(i)*D(n)
     $          +du(i) -diff(i)*D(n)*grav(n,1)-w(i)
        d0(n,i)= (dporo(i)+baro)*(-diff(i)*D(n)*grav(n,1))
     $      -ddiff(i)*D(n)*grav(n,1) - decay(n) 
c       write(38, '(4f12.7)')z(i), d2(n,i), d1(n,i), d0(n,i)
 65     continue
        d2(n,dimlevels) = diff(dimlevels)*D(n)
	d1(n,dimlevels)=(dporo(dimlevels)+baro)*(diff(dimlevels)*D(n))+
     $       ddiff(dimlevels)*D(n)
        d0(n,dimlevels)=  - decay(n) 
c       write(38, '(4f12.7)')z(i), d2(n,i), d1(n,i), d0(n,i)
        do 66 i=1,jlock-1
	d2l(n,i) = diffmidl(i)*(al+(1-al)*D(n))
	d1l(n,i)=(dporol(i)+baro)*(diffmidl(i)*(al+(1-al)*D(n)))+
     $       ddiffl(i)*(al+(1-al)*D(n))
        d0l(n,i)=  - decay(n) 
c       write(38,'(4f12.7)')depth(i),d2l(n,i), d1l(n,i), d0l(n,i)
c	d2l(n,i) = diffmidl(i)*D(n)
c	d1l(n,i)=(dporol(i)+baro)*(diffmidl(i)*D(n))+
c     $       ddiffl(i)*D(n)
c        d0l(n,i)=  - decay(n) 

 66   continue

c --- re-initialise time:
      ta=0
c ------------------------------------------------------------		  
c concentration initialisation
c ------------------------------------------------------------
c	GAS-SPECIFIC CHANGE 1 (ATMOSPHERIC HISTORY - automatically zero if the next 14 lines are commented out)
      if(iso) then
c	  open (41,file='co2_history2011.txt',status='old')     
c	  open (42,file='CH4_history.txt',status='old')      
          open (41,file='inputfiles/SCENARIO_isotopes.txt',status='old')
          open (42,file='inputfiles/SCENARIO_isotopes.txt',status='old')
          open (43,file='inputfiles/SCENARIO_isotopes.txt',status='old')
          open (44,file='inputfiles/SCENARIO_isotopes.txt',status='old')
          open (45,file='inputfiles/SCENARIO_isotopes.txt',status='old')

	  open (46,file='inputfiles/SCENARIO_NEEM08_CO2.txt',status='old')     
	  open (47,file='inputfiles/SCENARIO_NEEM08_CH4.txt',status='old')      
          open(48,file='inputfiles/SCENARIO_NEEM08_CH3CCL3.txt')
c	  open (44,file='SCENARIO_NEEM08_14CO2.txt',status='old')
c	  open (45,file='SCENARIO_NEEM08_HFC134a.txt',status='old')
c	  open (46,file='SCENARIO_NEEM08_CFC113.txt',status='old')
c	  open (47,file='SCENARIO_NEEM08_SF6.txt',status='old')
c	  open (48,file='SCENARIO_NEEM08_CFC12.txt',status='old')
c	  open (49,file='SCENARIO_NEEM08_CFC11.txt',status='old')
c	  open (50,file='SCENARIO_NEEM08_CCL4.txt',status='old')

c initialise concentrations:
	 read (40+n,*) O2newdate(n),O2new(n)
	 if (ta.ge.(O2newdate(n)-1800.)) then
 201        O2date(n) = O2newdate(n)
	     O2(n) = O2new(n)
	     read (40+n,*) O2newdate(n),O2new(n)
             if (ta.ge.(O2newdate(n)-1800.)) goto 201
c$$$                     print*, 'update concentrations',ta,O2date(n),
c$$$     $                    O2(n),O2new(n)
	  end if

         O2(n) = O2new(n)
         print*,'initialise concentrations',O2(n),'=',O2new(n),
     $        'date=', O2date(n)
         end if
c --  end just iso
         buffer(n) = 0.
         d1(n,1) = 0.
         d2(n,1) = 0.
         d0(n,1) = 0.
         deriv(n,1) = 0.
c set initial condition gas concentration profiles
         do 60 i=1,dimlevels
            c(n,i) = O2new(n)
c Initial condition for gravitationally & thermally fractionated gases
c         if (z(i).ge.conveczone) c(n,i) = (exp(grav(n)/dz*(z(i) - 
c     $	  conveczone)/aaK)-1)*1000 - 
c         if (z(i).ge.conveczone) c(n,i) = (exp(grav(n)*(z(i) - 
c     $	  conveczone))-1)*1000 - 
c     $	 geotherm*0.75*(i-1)*dz*(omega1/aaK-omega2/aaK**2) + O2new(n)

c--- for ratios --- 
         if (z(i).ge.conveczone) c(n,i) = O2new(n)*exp((grav(n,1) -
     $	   geotherm*(omega1(n)/aaK-omega2(n)/aaK**2)/1000)
     $	  *(z(i)-conveczone))

  60	 continue
c initialise for lock-in
        do 87 i = 1, jlock-1
           cl(n,i) = c(n,dimlevels)	
 87	continue 
c end concentration initialisation
c     Plot concentration
c        if(n.eq.1) then
c        open(51, file='Initialconcentration.51', status='replace')
c         do 93 jj=1,dimlevels
c            write (51, "(2g10.5)") z(jj),c(1,jj)
c 93      continue
c           do 92 i = 2, jlock-1
c               write (51, "(2g10.5)") depth(i),cl(1,i)
cc 920	format (f6.3,2x,f12.5,2x)
c 92         continue
c            end if
c update cprint:
                cprint=.true.

c ------------------------------------------------------------
c Surface temperature boundary condition:
	  open (27,file=fThist,status='old')
	  read (27,*) newTsurfdate,newTsurf
	  startTseries = newTsurfdate
c to disable the temperature history, set startTseries to 5000
	startTseries = 5000.
 112    Tsurfdate=newTsurfdate
             Tsurf=newTsurf
	     read (27,*) newTsurfdate,newTsurf
          if(starttime.ge.newTsurfdate) goto 112

c initialize temperature profile
c        do 88 i=1,tlevels
c	t(i) = aa + geotherm*dz*(i-1)
c 88     continue
c        open(27, file = fileiniT, status = 'old')
     	do 59 i=1,lt
c          read(27,*) zi(i), t(i)
        t(i)=tp(i)
 59     continue 
c       close(27)

c ------------------------------------------------------------
c       ========================================================
c       main timestepping loop
c       --------------------------------------------------------
c
c 		Top boundary condition simulation
c
      do 1000 jt = 1,ntstep
       time = jt*tstep
       ta=time/spa
       day=(ta-aint(ta))*365.25
			
c*********************************************************
c  calculate surface temperature over time (top boundary condition, model forcing)
       IF(withtemp.and.mod(jt,ktemp).eq.0) then

c	approximate generic Antarctic seasonal cycle (day 0 is Jan 0)
c	   t(1) = aa + 24.16
c	  if (day.gt.29..and.day.le.84.) t(1) =(aa+24.16) - 35./55.*(day-29)
c	  if (day.gt.84..and.day.le.279.) t(1) = aa - 10.84
c	  if (day.gt.279..and.day.le.329.) t(1)=(aa-10.84)+35./50.*(day-279)

c		perturbation to temperature
c	  if (ta.ge.100..and.ta.lt.200.) t(1) = t(1)-(ta-100)/100
c	  if (ta.ge.200..and.ta.lt.300.) t(1) = t(1)-1+(ta-200)/100

c      t(1) = aa + 17.1
c      if (day.gt.29..and.day.le.84.) t(1) = (aa+17.1) -25./55.*(day-29)
c      if (day.gt.84..and.day.le.279.) t(1) = aa - 7.9
c      if (day.gt.279..and.day.le.329.) t(1)= (aa-7.9) + 0.5*(day-279)

c     -or 2. sinusoidal seasonal cycle: add 0.0164 so that solstice 
c is on dec 25 and not jan 1
c      season=10.*(cos(2.*3.1415926*(ta-int(ta)+0.0164))+
c     $         0.3*cos(4.*3.1415926*(ta-int(ta)+0.0164)))
c          season=0
c		-or 3 NEEM seasonal cycle --- 
c	season= 13.5*sin((ta-floor(ta)+0.0274)*2.*3.1415927)
c     $         +7*sin((ta-floor(ta)+0.0274)*3.1415927)
	season=18.*(cos(2.*3.1415926*(ta-int(ta)+0.0164+0.5))+
     $         0.3*cos(4.*3.1415926*(ta-int(ta)+0.0164+0.5)))
c          season=-22*cos(2.*3.1415926*(ta-int(ta)))
			
c	Force model with surface temperature time series:
c			for Siple Dome, South Pole, or Megadunes
c	  if (ta.ge.startTseries-1800.) then
		  if (ta.ge.(newTsurfdate-starttime)) then
 200		 Tsurfdate = newTsurfdate
		     Tsurf = newTsurf
		     read (27,*) newTsurfdate,newTsurf
c			 if (ta.ge.(newTsurfdate-1800.)) go to 200
		  end if
c	      active interpolation:
c add this condition so that we dont extrapolate before we have any 
c data in the temperature time series file.
		if(ta.ge.Tsurfdate-starttime) then                  
		  aa = ((newTsurf-Tsurf)/(newTsurfdate-Tsurfdate)*
     $	    (ta-(Tsurfdate-starttime))+Tsurf) 	
                  endif
c	  end if
		
	
c   generic Greenland seasonal cycle (day 0 is Jan 0)
c      t(1) = aa - 10.
c      if (day.gt.208..and.day.le.271.) t(1)=(aa+21.3)-0.5*(day-208)
c      if (day.gt.150..and.day.le.208.) t(1) = aa + 21.3
c      if (day.gt.87..and.day.le.150.) t(1)=(aa-10.)+0.5*(day-87)
	
c$$$c --- add a triangle
c$$$                  at=0
c$$$              if(ta.ge.bb*jbb.and.ta.lt.(bb+1)*jbb) then
c$$$           at=(ta-bb*jbb)/jbb*amp
c$$$           end if
c$$$             if(ta.ge.(bb+1)*jbb.and.
c$$$     $         ta.lt.(bb+2)*jbb) then
c$$$           at=((bb+2)*jbb-ta)/jbb*amp
c$$$           end if
c$$$	
c$$$      t(1)=aa+season+at	
                  t(1)=aa+season

      do 301 i = 2, lclevel-1
       dt = (t(i+1)-t(i-1))/2/dz
       d2t = (t(i+1)-2.*t(i) + t(i-1))/dz/dz
       rr = ro(i)
       rk = (1-0.00882*(t(i)+30))*(-1.229d-14*rr**3
     $  + 2.1312d-11*rr**2-9.4d-9*rr + 1.779d-6)
       dkro = (1-0.00882*(t(i)+30))*(-3.687d-14*rr**2
     $  + 4.2624d-11*rr - 9.4d-9)
       a = 4.26d-13*exp(-7217/(273.15+t(i)))      
       cp = 2096.*(1+0.0037*t(i))
       qs = 2.*(bdot/rh)**(4./3.)/(rr*cp*a**(1./3.))
       vz = rhoice*bdot/rr*(1-sg(i)/rhoice/rh)
       a1 = (rk/rr+dkro)*dro(i)-vz    
       a2 = dkt(i) + rk*dct/cp
       a3 = sg(i)*bdot*rhoice/rr**3*dro(i)*g/cp + qs
       tnew(i) = t(i) + ktemp*tstep*(rk*d2t+a1*dt+a2*dt*dt+a3)
  301 continue
      i=lclevel
      dt = (t(i+1)-t(i-1))/(dz+dz1)
       d2t = ((t(i+1)-t(i))/dz1-(t(i)-t(i-1))/dz)/((dz+dz1)/2)
       rr = ro(i)
       rk = (1-0.00882*(t(i)+30))*(-1.229d-14*rr**3
     $  + 2.1312d-11*rr**2-9.4d-9*rr + 1.779d-6)
       dkro = (1-0.00882*(t(i)+30))*(-3.687d-14*rr**2
     $  + 4.2624d-11*rr - 9.4d-9)
       a = 4.26d-13*exp(-7217./(273.15+t(i)))      
       cp = 2096.*(1+0.0037*t(i))
       qs = 2.*(bdot/rh)**(4./3.)/(rr*cp*a**(1./3.))
       vz = rhoice*bdot/rr*(1-sg(i)/rhoice/rh)
       a1 = (rk/rr+dkro)*dro(i)-vz    
       a2 = dkt(i) + rk*dct/cp
       a3 = sg(i)*bdot*rhoice/rr**3*dro(i)*g/cp + qs
       tnew(i) = t(i) + ktemp*tstep*(rk*d2t+a1*dt+a2*dt*dt+a3)

      do 302 i = lclevel+1, level1-1
       dt = (t(i+1)-t(i-1))/2/dz1
       d2t = (t(i+1)-2.*t(i) + t(i-1))/dz1/dz1
       rr = ro(i)
       rk = (1-0.00882*(t(i)+30))*(-1.229d-14*rr**3
     $  + 2.1312d-11*rr**2-9.4d-9*rr + 1.779d-6)
       dkro = (1-0.00882*(t(i)+30))*(-3.687d-14*rr**2
     $  + 4.2624d-11*rr - 9.4d-9)
       a = 4.26d-13*exp(-7217/(273.15+t(i)))      
       cp = 2096.*(1+0.0037*t(i))
       qs = 2.*(bdot/rh)**(4./3.)/(rr*cp*a**(1./3.))
       vz = rhoice*bdot/rr*(1-sg(i)/rhoice/rh)
       a1 = (rk/rr+dkro)*dro(i)-vz    
       a2 = dkt(i) + rk*dct/cp
       a3 = sg(i)*bdot*rhoice/rr**3*dro(i)*g/cp + qs
       tnew(i) = t(i) + ktemp*tstep*(rk*d2t+a1*dt+a2*dt*dt+a3)
 302  continue

      i=level1
      dt = (t(i+1)-t(i-1))/(dz1+dz2)
       d2t = ((t(i+1)-t(i))/dz2-(t(i)-t(i-1))/dz1)/((dz1+dz2)/2)
       rr = ro(i)
       rk = (1-0.00882*(t(i)+30))*(-1.229d-14*rr**3
     $  + 2.1312d-11*rr**2-9.4d-9*rr + 1.779d-6)
       dkro = (1-0.00882*(t(i)+30))*(-3.687d-14*rr**2
     $  + 4.2624d-11*rr - 9.4d-9)
       a = 4.26d-13*exp(-7217./(273.15+t(i)))      
       cp = 2096.*(1+0.0037*t(i))
       qs = 2.*(bdot/rh)**(4./3.)/(rr*cp*a**(1./3.))
       vz = rhoice*bdot/rr*(1-sg(i)/rhoice/rh)
       a1 = (rk/rr+dkro)*dro(i)-vz    
       a2 = dkt(i) + rk*dct/cp
       a3 = sg(i)*bdot*rhoice/rr**3*dro(i)*g/cp + qs
       tnew(i) = t(i) + ktemp*tstep*(rk*d2t+a1*dt+a2*dt*dt+a3)


      do 303 i = level1+1, lt-1
       dt = (t(i+1)-t(i-1))/2/dz2
       d2t = (t(i+1)-2.*t(i) + t(i-1))/dz2/dz2
       rr = ro(i)
       rk = (1-0.00882*(t(i)+30))*(-1.229d-14*rr**3
     $  + 2.1312d-11*rr**2-9.4d-9*rr + 1.779d-6)
       dkro = (1-0.00882*(t(i)+30))*(-3.687d-14*rr**2
     $  + 4.2624d-11*rr - 9.4d-9)
       a = 4.26d-13*exp(-7217/(273.15+t(i)))      
       cp = 2096.*(1+0.0037*t(i))
       qs = 2.*(bdot/rh)**(4./3.)/(rr*cp*a**(1./3.))
       vz = rhoice*bdot/rr*(1-sg(i)/rhoice/rh)
       a1 = (rk/rr+dkro)*dro(i)-vz    
       a2 = dkt(i) + rk*dct/cp
       a3 = sg(i)*bdot*rhoice/rr**3*dro(i)*g/cp + qs
       tnew(i) = t(i) + ktemp*tstep*(rk*d2t+a1*dt+a2*dt*dt+a3)
 303  continue
		
c    calculate the temperature at the bottom of the first box
c not sure what this is used for -Anais
  		tmid = (t(1)+t(2))/2+273.1
c		calculate gas flux into first box down

c ---- assign new temperatures --- 

       do 550 i = 2, lt-1
        t(i) = tnew(i)
  550  continue
c       if(jt/ktemp.le.20) print*, t(1), t(2), t(3), t(4), t(5)
c --- update the derivatives ---- 
       i=1
       grav(n,i)=Dm(n)/1000*g/R/(273.16+t(i))
c        therm(n,i)=(omega1(n)/(273.16+t(i))-omega2(n)/(273.16+t(i))**2
c     $          -1/(273.16+t(i)))*(t(i+1)-t(i))/dz/1000
        therm(n,i)=(omega1(n)/(273.16+t(i))-omega2(n)/(273.16+t(i))**2
     $          )*(t(i+1)-t(i))/dz/1000
c       grav(n,i)=Dm(n)/1000*g/R/aaK
c        therm(n,i)=(omega1(n)/aaK-omega2(n)/(aaK)**2
c     $         -1/aaK )*(t(i+1)-t(i))/dz/1000
        dgravT(n,i) = (grav(n,i+1)-grav(n,i))/dz
        dthermT(n,i)=(therm(n,i+1)-therm(n,i))/dz
c        dthermT(n,i)=0

        d2(n,i) = diff(i)*D(n) + u(i)
	d1(n,i) = (dporo(i)+baro)*(diff(i)*D(n)+u(i))+ddiff(i)*D(n)
     $          +du(i) +diff(i)*D(n)*(-grav(n,i)+therm(n,i))-w(i)
        d0(n,i)= (dporo(i)+baro)*(diff(i)*D(n)*(-grav(n,i)+therm(n,i)))
     $    +ddiff(i)*D(n)*(-grav(n,i)+therm(n,i))
     $    +diff(i)*D(n)*(-dgravT(n,i)+dthermT(n,i))
     $       - decay(n) 

        do 365 i=2,dimlevels-1
        grav(n,i)=Dm(n)/1000*g/R/(273.16+t(i))
c        therm(n,i)=(omega1(n)/(273.16+t(i))-omega2(n)/(273.16+t(i))**2
c     $          -1/(273.16+t(i)))*(t(i+1)-t(i-1))/2/dz/1000
        therm(n,i)=(omega1(n)/(273.16+t(i))-omega2(n)/(273.16+t(i))**2
     $          )*(t(i+1)-t(i-1))/2/dz/1000
c       grav(n,i)=Dm(n)/1000*g/R/aaK
c        therm(n,i)=(omega1(n)/aaK-omega2(n)/(aaK)**2
c     $          +1/aaK)/1000*(t(i+1)-t(i-1))/dz/2
        dgravT(n,i) = (grav(n,i+1)-grav(n,i-1))/dz/2
        dthermT(n,i)=(therm(n,i+1)-therm(n,i-1))/2/dz
c        dthermT(n,i)=0
        d2(n,i) = diff(i)*D(n) + u(i)
	d1(n,i) = (dporo(i)+baro)*(diff(i)*D(n)+u(i))+ddiff(i)*D(n)
     $          +du(i) +diff(i)*D(n)*(-grav(n,i)+therm(n,i))-w(i)
        d0(n,i)= (dporo(i)+baro)*(diff(i)*D(n)*(-grav(n,i)+therm(n,i)))
     $    +ddiff(i)*D(n)*(-grav(n,i)+therm(n,i))
     $    +diff(i)*D(n)*(-dgravT(n,i)+dthermT(n,i))
     $       - decay(n) 
 365     continue
c at dimelvels, grav and therm go to zero, no more gravitational diffusion
        d2(n,dimlevels) = diff(dimlevels)*D(n)
	d1(n,dimlevels)=(dporo(dimlevels)+baro)*(diff(dimlevels)*D(n))+
     $       ddiff(dimlevels)*D(n)-w(i)
        d0(n,dimlevels)=  - decay(n) 

        do 366 i=1,jlock-1
	d2l(n,i) = diff(levels+i)*(al+(1-al)*D(n))
	d1l(n,i)=(dporol(i)+baro)*(diff(levels+i)*(al+(1-al)*D(n)))+
     $       ddiffl(i)*(al+(1-al)*D(n))
        d0l(n,i)=  - decay(n) 
c	d2l(n,i) = diffmidl(i)*D(n)
c	d1l(n,i)=(dporol(i)+baro)*(diffmidl(i)*D(n))+
c     $       ddiffl(i)*D(n)
c        d0l(n,i)=  - decay(n) 

 366   continue
c ---- output the temperature ---
c  .32: write a time series at the surface and lock-in depth:
	    if (int(jt/kprint)*kprint.eq.jt) then	
c             if(ta.le.samplestart+10) then
	  write (32,723) ta,aa,t(1),t(40),
     $              t(100),t(200),t(300),t(lt)
 723      format(2x,f10.4,2x,f10.4,2x,f10.5,2x,2x,f10.5,2x,
     $         f10.5,2x,f10.5,2x,f10.5,2x,f10.5)      
        end if

             end if 
c ------ end of the update of the temperature model ---

c	Force model with atmospheric concentration time series:
c			for NEEM, Siple Dome, South Pole, or Megadunes
         if(iso) then
		  if (ta.ge.(O2newdate(n)-starttime)) then
		     O2date(n) = O2newdate(n)
		     O2(n) = O2new(n)
		     read (40+n,*) O2newdate(n),O2new(n)
c$$$                     print*, 'update concentrations',ta,O2date(n),
c$$$     $                    O2(n),O2new(n)
		  end if

c	active interpolation:
		  c(n,1)=((O2new(n)-O2(n))/(O2newdate(n)-O2date(n))*
     $	    (ta-(O2date(n)-1800.))+O2(n))
        end if
c*********************************************************

c		calculate gas flux into first box down

c 310	gradout = (c(n,2) - c(n,1))/dz	
 310	grad(1) = (c(n,2) - c(n,1))/dz	
        if(isnan(grad(1))) print*, 'NaN alert grad(1)=',grad(1), jt
c        if(jt.le.10) print*,'jt=',jt, ' grad(1)=', grad(1)
c
c     -------------MAIN LOOP THROUGH OPEN FIRN COLUMN----------
        DO 500 i = 2, levels
c 
c    	gradin = gradout
		grad(i) = (c(n,i+1) - c(n,i))/dz
       if(isnan(grad(i))) then 
          print*, 'NaN alert 1121 grad', i, jt,'z=', z(i)
          print*, 'c(n,i+1)',c(n,i+1)
          print*, 'c(n,i)',c(n,i)
c          print*, 'dz',dz
          goto 1300
          end if
       if(grad(i).ge.1000) then 
          print*, ' alert 1121 grad',grad(i),'i=',i,'jt=',jt,'z=',z(i)
          print*, 'c(n,i+1)',c(n,i+1)
          print*, 'c(n,i)',c(n,i)
c          print*, 'dz',dz
c          goto 1300
          end if
        if(grad(i).le.-1000) then 
          print*, ' alert 1121 grad',grad(i),i,jt,'z=',z(i)
          print*, 'c(n,i+1)',c(n,i+1)
          print*, 'c(n,i)',c(n,i)
c          print*, 'dz',dz
c          goto 1300
          end if
 
 500      continue

          do 501 i=2,levels
c             if(jt.eq.660.and.i.eq.15) then
c          print*, 'grad(i-1)=', grad(i-1)
c          print*, 'grad(i)=', grad(i)
c          print*, 'c(n,i)=', c(n,i)
c          print*, 'c(n,i-1)',c(n,i-1)
c          print*, 'c(n,1)',c(n,1)
c          print*, 'c(n,2)',c(n,2)
c          print*, 'c(n,3)',c(n,3)
c          end if
          deriv(n,i) = d2(n,i)*(grad(i)-grad(i-1))/dz 
     $  + d1(n,i)*(grad(i-1) + grad(i))/2 +d0(n,i)*c(n,i)
       if(isnan(deriv(n,i))) then
          print*, 'NaN alert 1123 deriv', i, jt
          print*, 'grad(i-1)=', grad(i-1)
          print*, 'grad(i)=', grad(i)
          print*, 'grad(i+1)=', grad(i+1)
          print*, 'c(n,i)=', c(n,i)
          print*, 'c(n,i-1)',c(n,i-1)
c          print*, 'dz=', dz
          call abort
          end if

c                 if(i.le.20.and.jt.le.10) print*, grad(i)
c       if(jt.eq.1) write(52, '(3G15.5)') z(i), c(1,i), deriv(n,i)
c --------------------------------------------------------------------
				
 501    CONTINUE
c       ---------------------------------------------------------
c     assign new concentrations in open firn
c       ---------------------------------------------------------
c		Gases in open firn:
        do 600 i = 2, dimlevels
             c(n,i) = deriv(n,i)*tstep + c(n,i)
        if(c(n,i).ge.10000) then
           print*,'inf 1250 c', c(n,i),'deriv',
     $            deriv(n,i)*tstep,'tst',tstep, n,i,jt
           print*, 'd1=', d1(n,i)
           print*, 'd2=', d2(n,i)
           print*, 'd0=', d0(n,i)
          print*, 'grad(i-1)=', grad(i-1)
          print*, 'grad(i)=', grad(i)
          if(i.gt.dimlevels)print*, 'grad(i+1)=', grad(i+1)
          print*, 'c(n,i-1)=', c(n,i-1)
          print*, 'c(n,i)=', c(n,i)
          print*, 'c(n,i+1)=', c(n,i+1)

           call abort
           end if
        if(c(n,i).le.-10000) then
           print*,'inf 1250 c', c(n,i),'deriv',
     $            deriv(n,i)*tstep,'tst',tstep, n,i,jt
          print*, 'd1=', d1(n,i)
           print*, 'd2=', d2(n,i)
           print*, 'd0=', d0(n,i)
          print*, 'grad(i-1)=', grad(i-1)
          print*, 'grad(i)=', grad(i)
          print*, 'grad(i+1)=', grad(i+1)
          print*, 'c(n,i)=', c(n,i)

           call abort
           end if

c        if(jt.le.40.and.i.eq.15) then
c           print*,'c(n,i)', c(n,i),'deriv',
c     $            deriv(n,i)*tstep, n,i,jt
c           print*, 'd1=', d1(n,i)
c           print*, 'd2=', d2(n,i)
c           print*, 'd0=', d0(n,i)
c          print*, 'grad(i-1)=', grad(i-1)
c          print*, 'grad(i)=', grad(i)
c          print*, 'grad(i+1)=', grad(i+1)
c          print*, 'c(n,i)=', c(n,i)
c          print*, 'c(n,i-1)=', c(n,i-1)
c          print*, 'c(n,i+1)=', c(n,i+1)
c          end if

  600   continue
	
c --------------------------------------------------------------------
c     box at bottom of open firn domain (with quasi-impermeable layer and lock-in zone beneath it)

c	gradin = gradout
	grad(dimlevels) = (cl(n,2)-c(n,dimlevels))/thick(1)
	Fin = -dlockporin*dlockdiffin*(al+(1-al)*D(n))*grad(levels)
	Fout = -pormidl(1)*diffmidl(1)*(al+(1-al)*D(n))*grad(dimlevels)
	deriv(n,dimlevels) = (Fin - Fout)/sopen(1)/dlockz
     $   - w(dimlevels)*grad(levels)		
c with radioactive decay
c     $   - decay*c(n,dimlevels)	 		

c$$$c --- anais updates the expression for the lock in box --- 
c$$$        dporo(dimlevels)= (pormidl(1)-dlockporin)/sopen(1)/dlockz
c$$$c     dporo is equivalent to:
c$$$c     dporo(dimlevels)=(sopen(2)-poro(levels))/(dz+thick(1))/sopen(1)
c$$$        ddiff(dimlevels)= (diffmidl(1)-dlockdiffin)/dlockz
c$$$c     ddiff is equivalent to:
c$$$c       ddiff(dimlevels)= (diff(levels+2)-diff(levels))/(dz+thick(1))
c$$$        i=dimlevels
c$$$        d2(n,i) = diff(i)*D(n) + u(i)
c$$$	d1(n,i) = (dporo(i)+baro)*(diff(i)*D(n)+u(i))+ddiff(i)*D(n)
c$$$     $          +du(i) -diff(i)*D(n)*grav(n)-w(i)
c$$$        d0(n,i)= (dporo(i)+baro)*(-diff(i)*D(n)*grav(n))
c$$$     $      -ddiff(i)*D(n)*grav(n) - decay(n) 
c$$$
c$$$        do 510 i=2,dimlevels
c$$$        deriv(n,i) = d2(n,i)*(grad(i)-grad(i-1))/dz + d1(n,i)*
c$$$     $  (grad(i-1) + grad(i))/2 +d0(n,i)*c(n,i)
c$$$ 510    continue

c --------

c	  keep track of concentration in top box of lock-in zone, cl(1), which is
c		an average of c(dimlevels) over the year-bundle.  "Buffer" is just an
c		index for keeping track of this averaging, and it gets reset to
c		zero each year-bundle.
	cl(n,1) =(cl(n,1)*buffer(n)+c(n,dimlevels))/(buffer(n)+1.)
	buffer(n) = buffer(n) + 1.
c	 if(jt.le.10) print*,'jt=',jt, ' deriv(level)=', deriv(n,levels),
c     $       'deriv(n,dimlevels)=', deriv(n,dimlevels)		     
c *****************************************************************
c TIME-DEPENDENT SIMULATION OF LOCK-IN ZONE EVOLUTION
c Skip time steps to make model faster (skip klockstep-1 steps)
       if ((jt/klockstep)*klockstep.eq.jt) then		
				
c calculate firn air composition in open pores, with eddy diffusion term
c to simulate slow viscous mixing within the lock-in zone			
c --------------------------------------------------------------------
		
        do 610 i=2, jlock-1
           gradl(i)=(cl(n,i+1)-cl(n,i))/thick(i)
		Fin = Fout
        Fout=-pormidl(i)*diffmidl(i)*(al+(1-al)*D(n))*
     $               (cl(n,i+1)-cl(n,i))/thick(i)
        if(isnan(Fout)) print*,'NaN alert Fout 1185', i
        derivlock(n,i) = (Fin-Fout)/(thick(i)/2+thick(i-1)/2)/sopen(i)
        if(isnan(derivlock(n,i))) print*,'NaN alert derivlock 1185', i

c        derivlock(n,i) = d2l(n,i)*(gradl(i)-gradl(i-1))/(thick(i)/2.+
c     $   thick(i-1)/2.)+ d1l(n,i)*(gradl(i-1) + gradl(i))/2.
c     $          +d0l(n,i)*cl(n,i)
  610  continue	
	
c assign new concentrations in lock-in zone
	do 615 i=2, jlock-1
	cl(n,i) = cl(n,i) + derivlock(n,i)*tstep*klockstep
        if(cl(n,i).ge.10000)then
           print*,'inf 1222 cl', cl(n,i), n,i,jt, jt/klockstep, ta
           call abort
           end if
        if(cl(n,i).le.-10000) then
           print*,'inf 1222 cl', cl(n,i), n,i,jt, jt/klockstep, ta
           call abort
           end if
  615   continue

c Simluate advection due to bubble trapping by moving values down one temporal layer each year		
c calculate the time of year (a number between 0 and 1)
           if ((jt/kbundle)*kbundle.eq.jt) then 
	        buffer(n) = 0
c advect the gases downward by one layer to simulate the accumulation/bubble trapping process
     	     do 650 i = 1, jlock-2
                ii = jlock-i
                cl(n,ii) = cl(n,ii-1) 				
  650	     continue	 
	   end if

	end if
 			

c			go to 1200
c        end lock-in zone time-dependent model
c      ===============================================================
c       Output sampling routines for the open firn: 

	if(jt.eq.1) then
c            fin0 = ofolder(1:charnum) //
c     $         "/WAISfirnoutput10_test1_g0_1000.33"
c          nfin=charnum+31
c        bb=0
        indname = 1000   
        write(fin0(nfin-3:nfin), '(I4)') indname
       write(fin0(nfin-5:nfin-5), '(I1)') n
           fin5 = fin0(1:nfin) // ".5"//fin0(nfin-5:nfin-5)
         open(50+n, file=fin5, status = 'replace')
         print*, fin5
         open(32, file=fin0(1:nfin)//'.32')

         end if
c  write a time series at the surface and lock-in depth:
c	  if (ta.gt.samplestart) then
c
c --------------------------------------------------------------------
	    if ((jt/kprint)*kprint.eq.jt) then		
		
	  write(50+n,720)ta,c(n,1)-1, c(n,10)-1, c(n,40)-1,
     $              c(n,levels)-1,c(n,dimlevels)-1,cl(n,1)-1,cl(n,2)-1
 720  format(2x,f12.5,2x,g12.5,2x,g12.5,2x,g12.5,2x,g12.5,2x,g12.5,2x,
     $  g12.5,2x,g12.5,2x,g12.5) 

c$$$     
c$$$      write(60+n,721)ta,c(n,2)-c(n,1),
c$$$     $     c(n,dimlevels-4)-c(n,dimlevels-3),
c$$$     $     c(n,dimlevels-3)-c(n,dimlevels-2),
c$$$     $     c(n,dimlevels-2)-c(n,dimlevels-1),
c$$$     $              c(n,dimlevels)-c(n,dimlevels-1),cl(n,2)-cl(n,1),
c$$$     $     cl(n,6)-cl(n,5)
c$$$ 721  format(2x,f12.5,2x,g12.6,2x,g12.6,2x,g12.6,2x,g12.6,2x,g12.6,2x,
c$$$     $  g12.6,2x, g12.6,2x)  
c$$$        end if

	  end if

        if(ta.ge.tendc.and.cprint) then
          do  731 jj=1,dimlevels
           csave(n,jj)=c(n,jj)
 731      continue
          do 732 i=2,jlock
            clsave(n,i)=cl(n,i)
 732      continue
          cprint=.false.
          print*, 'printing csave, ta=', ta, 'tendc=', tendc
        end if

c$$$        if(jt.eq.2) then
c$$$           open(35, file='testjt2.35', status='replace')
c$$$           do 901 jj=1,dimlevels
c$$$            write (35,825) z(jj),c(1,jj),c(2,jj)
c$$$ 825        format (f6.3,2x,f12.5,2x,f10.4,2x)
c$$$ 901     continue
c$$$ 			
c$$$c complete the depth profile with the lock-in zone:
c$$$
c$$$           do 951 i = 2, jlock-1
c$$$	   write (35,925) depth(i),cl(1,i),cl(2,i)
c$$$ 925       format (f6.3,2x,f12.5,2x,f12.5,2x)
c$$$ 951    continue
c          end if

c write a depth profile at the end of the run, on July 16, 2008
c		 if (ta.gt.(length-0.462)) then
c       if(ta.gt.10) then
c old profile
c           fin0='NEEMfirnoutput.60'
c           write(fin0(16:17), '(I2)') 60+n 
c          open(60+n, file=fin0, status='replace')
c$$$         do 901 jj=1,dimlevels
c$$$	  write (60+n,801) z(jj),t(jj),openpor(jj),c(n,jj),vb(jj),w(jj)
c$$$ 801      format (2x,f7.2,2x,f9.4,2x,f13.6,2x,f12.5,2x,f13.6,2x,e10.3)
c$$$
c$$$ 901   continue
c$$$ 			
c$$$c complete the depth profile with the lock-in zone:
c$$$
c$$$           do 951 i = 2, jlock-1
c$$$	  write (60+n,921) depth(i), i*1.,sopen(i),cl(n,i),vbl(i),wl(i)
c$$$ 921      format (2x,f7.2,2x,f9.4,2x,f10.2,2x,f10.2,2x,f10.2,2x,e10.3)
c$$$
c$$$ 951   continue
c$$$
c        print*,'gas no,n,c(n,1),c(n, dimlevels), O2new(n),O2date(n), ta'	
c        print*, 'gas no', n, c(n,1), c(n, dimlevels), O2new(n)
c     $       ,o2date(n), ta+1800
c		 go to 1101	
c    	end if
c --- write a temperature time series:
c	    if (mod(jt, kprint*20).eq.0) then	
cc             if(ta.le.samplestart+10) then
c	  write (32,722) ta,t(1),t(2),t(3), t(4),t(40),
c     $              t(100),t(200),t(300)
c 722      format(2x,f10.4,2x,f10.5,2x,f10.5,2x,f10.5,2x,
c     $         f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)      
c        end if

c .35: write the annual average surface temperature	
c	    if (int(ta).eq.ta) then		
	    if (mod(jt,tst).eq.0) then		
	  write (35,723) ta,aa+at,aa,
     $              t(100),t(200),t(300),t(tlevels)
        end if


 1000   CONTINUE
c --------------------------------------------------------------------
c       end of main time stepping loop
 1101   close(40+n)
c 1101   close(41)
c        close(42)
        close(27)
        close(50+n)
        close(35)
        close(32)
        print*,'gas no,n,c(n,1),c(n, dimlevels), O2new(n),O2date(n), ta'	
        print*, 'gas no', n, c(n,1), c(n, dimlevels), O2(n)
     $       ,o2date(n), ta+1800

c .39: write a depth profile at the end of the run, 

      do 750 jj=1,lt
        write (39,810) z(jj), t(jj), ro(jj)
 810    format (2x,f7.2,2x,f10.5,2x,f10.4,2x,f10.5,2x,f10.7,2x,f10.5)
 750  continue
        close(39)
c .38: save the derivatives ar the end of the run (well, we just want ot check once what they're like)
        do 751 i=1,dimlevels
       write(38, '(4f12.7)')z(i), d2(n,i), d1(n,i), d0(n,i)
 751  continue
       do 752 i=1,jlock
       write(38,'(4f12.7)')depth(i),d2l(n,i), d1l(n,i), d0l(n,i)
 752  continue
       close (38)
c       =======================================================
 1100   CONTINUE
c		start again on the next gas
c       =======================================================

c     Output:
c        fin0='NEEMfirnoutput10_test1_1000.33'
c        nfin=27
c	if(jt.eq.1) then
c            fin0 = ofolder(1:charnum) //
c     $         "/WAISfirnoutput10_test1_1000.33"
c          nfin=charnum+28
c        bb=0
        indname = 1000 	  
        write(fin0(nfin-3:nfin), '(I4)') indname
           fin33 = fin0(1:nfin) // ".33"
c           fin34 = fin0(1:nfin) // ".34"
         open(33, file=fin33, status = 'replace')
c         open(34, file=fin34, status = 'replace')
         print*, fin33
c         end if
c     write depth profile in diffusive column
         do 900 jj=1,dimlevels
c            write (33,800) z(jj),c(1,jj),c(2,jj)
c     ,c(3,jj),c(4,jj),
c     $	 c(5,jj),c(6,jj) ,c(7,jj),c(8,jj), c(9,jj), c(10,jj)
c 800    format (f6.3,2x,f12.5,2x,f12.5,2x)
c     ,f12.5,2x,f13.6,2x,f13.6,
c     $   2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5,2x,f12.5,2x)
c            write (33,800) z(jj),(csave(1,jj)-1)*1000
            write (33,800) z(jj),(csave(1,jj)-1)*1000
     $          , (csave(2,jj)/csave(3,jj)-1)*1000
     $          , (csave(4,jj)/csave(5,jj)-1)*1000
     $          , (csave(2,jj)-1)*1000
     $          , (csave(4,jj)-1)*1000
     $          , csave(6,jj),csave(7,jj),csave(8,jj) 

c            write (33,800) z(jj),(c(1,jj)-1)*1000
c            write (33,800) z(jj),c(1,jj)
 800    format (f6.3,2x,8(f12.5,2x))
 900    continue
 			
c complete the depth profile with the lock-in zone:

           do 950 i = 2, jlock-1
c	   write (33,920) depth(i),cl(1,i),cl(2,i)
c     cl(3,i),cl(4,i),cl(5,i)
c      $	 ,cl(6,i) ,cl(7,i),cl(8,i), cl(9,i), cl(10,i)
c  920	format (f6.3,2x,f12.5,2x,f12.5,2x)
c ,f12.5,2x,f13.6,2x,f13.6,2x,
c     $   f9.5,2x,f9.5,2x,f12.5,2x,f12.5,2x,f12.5,2x)
               write (33,920) depth(i),
c     $             (cl(1,i)-1)*1000
     $             (clsave(1,i)-1)*1000
     $          , (clsave(2,i)/clsave(3,i)-1)*1000
     $          , (clsave(4,i)/clsave(5,i)-1)*1000
     $          , (clsave(2,i)-1)*1000
     $          , (clsave(4,i)-1)*1000
     $          , clsave(6,i), clsave(7,i)
     $             , clsave(8,i)

c               write (33,920) depth(i),cl(1,i)
 920	format (f6.3,2x,8(f12.5,2x))
  950	   continue
           close(33)
           print*, 'depth(jlock-1)', depth(jlock-1), 
     $          'depth(jlock)=', depth(jlock) 
c write the diffusivity used in this run, for sake of clarity
c           do 960 ii=1,ldimtotal
c              write(34, FMT='(1x,f12.5,2x,E12.5,1x)')zdiff(ii),diff(ii)
c 960          continue
c              close(34)
 1200      continue
      print *, char(7)
 1300	end
