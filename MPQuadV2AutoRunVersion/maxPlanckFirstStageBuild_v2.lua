--[[
This code will attempt to recreate the experiment done by the Max-Plack institute for quantum optics
(Rempe et al, Two-Dimensional trapping of dipolar molecules in time-varying electric fields).
Originally written by Dr. Demitri Balabanov, adapted for MPI simulation by Taylor Grubbs
This code creates a quadrupole geometry in simion and calculates the values for the electrostatic potential,  electric field magnitude, and stark potential.
--]]


-- ************************************************************
-- *********	Constants, Units and Conversions	***********
-- ************************************************************

-- MKS units
  m = 1.0 -- meter
  s = 1.0 -- second
  kg = 1.0 -- kilogram
  joule = kg * (m ^ 2.0) / (s ^ 2.0)
  kelvin = 1.0
  coulomb = 1.0
  volt = 1.0

-- Constants
  c = 299792458 * m / s
  h = 6.6260755 * 10.0 ^ -34.0 * joule * s
  kb = 1.380658 * 10.0 ^ -23.0 * joule / kelvin

-- Convenient units
  cm =  10.0 ^ -2.0 * m
  mm =  10.0 ^ -3.0 * m
  um = 10.0 ^ -6.0 * m
  amu = 1.6605402 * 10.0 ^ -27.0 * kg
  nsec = 10.0 ^ -9.0 * s
  kV = 10.0 ^ 3.0 * volt
  MHz = 10.0 ^ 6.0 * s ^ - 1.0
  e = 1.60217733 * 10.0 ^ -19.0 * coulomb
  eV = e * volt
  debye = (1.0 / c) * 10.0 ^ -21.0 * coulomb * (m ^ 2.0) / (s ^ 1.0)
  deb2wvn = (0.016791998 * cm ^ -1.0 / debye) * (kV/cm) ^ -1.0
  wavenumber = 1.24*10^-4 -- # of eV's in a wavenumber

--[[ notes:
  Debye to wavenumber conversion (deb2wvn) is for consistency with Alternate Grating Review paper
  Reference : 2006 J. Phys. B: At. Mol. Opt. Phys. 39 R263 (http://iopscience.iop.org/0953-4075/39/16/R01)
  (Debye / h*c)* kV / cm = 0.0167919918 cm^-1
--]]

-- ************************************************************
-- *****	Parameters of the Molecule(s) of Interest	*******
-- ************************************************************

  MolType = "ND3"	--Max Planck experiment uses ND3 molecule				-- a label to keep track of the type of molecule used
  mu = 1.47 * debye										-- permanent dipole moment [Debye] changed to 1.47 for ND3

  starkSlope = 1.44*10^-8 --assumed stark slope (1.2 cm^-1) for a certain ND3 state used by MPI. Converted into appropriate units for Simion eV/(V/mm) negative for high-field seeking, positive for low-field seeking

  mass = (20)*amu							-- [kg]
  constB = 44315.975 * MHz / c 								-- rotational B constant [m^-1]

  qnumJ = 1.0  												-- rotational J quantum number
  qnumM = 1.0  												-- M-sub-J quantum number; projection of total angular momentum (excluding nuclear spin) on space-fixed axis
  qnumK = 1.0  												-- quantum number of component of angular momentum along molecular axis

  vmax = 200 												-- [m/s]; initial (entering) speed of a synchronous particle in a molecular beam

--[[ notes:
The rotational constant for HCN is from JPL catalogue.
The observed transitions are from F. C. De Lucia and W. Gordy, 1969, Phys. Rev. 187, 58, and from F. C. Van den Heuvel, W. L. Meerts, and A. Dymanus, 1982, Chem. Phys. Lett. 92, 215.
The dipole moment is from G. Tomasevich, 1970, Thesis, Harvard University.

Note that |m-sub-j|=J states are always strong field seeking.
K=0 for linear (ex. diatomic) and symmetric top molecules.
--]]

-- ************************************************************
-- *************	Parameters of the Device	***************
-- ************************************************************

xsize = (300+1)					-- +1 b/c, actual PA length in WB is xsize-1, in addition, 0 is treated differently (LUA?) as far as the "fill" is
ysize = (300+1)					-- concerned (causing asymmetry in the created electrodes) hence, I added 1 grid cell everywhere along the PA boundary
zsize = 1

local quadVoltage = 5000 --  voltage magnitude of quadrupoles. ENTER THE EXACT VALUE THAT YOU WILL USE IN THE REAL WORLD EXPERIMENT NO SCALING NEEDED

scaleFactor = 100 --reduces PA size parameters by this factor; 100 converts each GU into 10 um which is needed for the Max-Planck experimental setup (want 1mm radius rods)
--voltage and other potential calculations will not be affected by this value.

local radius = 1 --radius of rods is 1 mm = 100 10^-5 m

-- ************************************************************
-- *************	Stark Energy Functions	*******************
-- ************************************************************

  function starkLinearShift(mu,Emag) 						-- stark energy shift linear with the applied electric field
	return -mu*Emag											-- [Joules]
  end

  function starkLinearShiftSlope(starkShift,Emag) --Calculates stark shift based on slope of state
      return starkShift*Emag
  end

  function stark1stOrder(mu,Emag,qnumJ,qnumM,qnumK)			-- stark energy from 1st order perturbation for symmetric-top wave function [Joules]; not including hyperfine structure
    if qnumJ == 0 then
		sp1 = 'Invalid: division by zero (J=0)'				-- status parameter #1
		return 0											-- invalid operation; division by 0
	else
		sp1 = 'OK' 											-- status parameter #1
		return -mu * Emag * qnumM * qnumK / (qnumJ*(qnumJ+1))	-- note that this vanishes for linear molecules (since K=0) and systems with
	end														-- nondegenerate energy levels (often lifted by inversion as in NH3)
  end

  function stark2ndOrder(mu,Emag,qnumJ,qnumM,qnumK,constB)	-- stark energy from 2nd order perturbation [Joules]; not including hyperfine structure
	if qnumJ == 0 and qnumK == 0 and qnumM == 0 then		-- relevant for linear molecules
		return -(mu^2) * (Emag^2) / (6 * h * constB * c)		-- see Townes and Schawlow
	else
		return ((mu^2) * (Emag^2) / (2 * h * constB*c))*((qnumJ^2 - qnumK^2)*(qnumJ^2 - qnumM^2) / ((qnumJ^3)*(2*qnumJ-1)*(2*qnumJ+1))-(((qnumJ+1)^2 - qnumK^2)*((qnumJ+1)^2-qnumM^2)/(((qnumJ+1)^3)*(2*qnumJ+1)*(2*qnumJ+3))))
	end
  end

  function starkOrder(mu,Emag,qnumJ,qnumM,qnumK,constB)		-- combined 1st+2nd order i.e. stark energy up to second order [Joules]
	if qnumJ == 0 and qnumK == 0 and qnumM == 0 then		-- second order; special case
		return -(mu^2) * (Emag^2) / (6 * h * constB*c)
	end
	if qnumJ ~= 0 and qnumK == 0 or qnumM == 0 then			-- second order only; first order effect = 0
		return ((mu^2) * (Emag^2) / (2 * h * constB*c))*((qnumJ^2 - qnumK^2)*(qnumJ^2 - qnumM^2) / ((qnumJ^3)*(2*qnumJ-1)*(2*qnumJ+1))-(((qnumJ+1)^2 - qnumK^2)*((qnumJ+1)^2-qnumM^2)/(((qnumJ+1)^3)*(2*qnumJ+1)*(2*qnumJ+3))))
	end
	if qnumJ ~= 0 and qnumK ~= 0 and qnumM ~= 0 then		-- first and second order; note that second order effect is typically 2 or 3 orders of magnitude smaller than first order effect
		return (-mu * Emag * qnumM * qnumK / (qnumJ*(qnumJ+1))) + ((mu^2) * (Emag^2) / (2 * h * constB*c))*((qnumJ^2 - qnumK^2)*(qnumJ^2 - qnumM^2) / ((qnumJ^3)*(2*qnumJ-1)*(2*qnumJ+1))-(((qnumJ+1)^2 - qnumK^2)*((qnumJ+1)^2-qnumM^2)/(((qnumJ+1)^3)*(2*qnumJ+1)*(2*qnumJ+3))))
	end
  end

  function starkPendular(mu,Emag,qnumJ,qnumM,constB) 		-- High field, pendular-state model [Joules]
		vp = 2.0 * qnumJ - math.abs(qnumM)					-- pendular-state quantum number (good quantum # in high field limit) [dimensionless]
		lambda = mu*Emag*1000/(constB*h*c)					-- dimensionless coupling parameter; lambda << 1 is weak-field limit; lambda >> 1 is strong-field limit; *1000 to work out the units
    return (-lambda+(vp + 1.0)*(2.0 * lambda)^0.5)*constB*h*c
  end

  function deltaKE(vmax,mass,nelect) 						-- kinetic energy to be remove by each field-stage [joules]
    return 0.5 * mass * (vmax ^ 2.0) / nelect
  end


--[[ notes:

--]]

-- ************************************************************
-- *********	Auxiliary Variables and Parameters	 **********
-- ************************************************************

local rscale = 10000					-- 10000; scaling used so that Stark Potentials are visible in SIMION WB (for visual purposes only).

-- ************************************************************
-- *****************	Program Code	 **********************
-- ************************************************************


local pa = simion.pas:open()

pa:size(xsize,ysize,zsize)

pa.dx_mm = 1/scaleFactor
pa.dy_mm = 1/scaleFactor
pa.dz_mm = 1/scaleFactor

--reduced values of xsize and ysize because grid unit size was decreased but
--simion still treats all numbers in terms of millimeters when creating the electrodes

local xSizemm = (xsize) / scaleFactor
local ySizemm = (ysize) / scaleFactor

--This actually creates the electrodes in the pa file
pa:fill{ function(xg,yg,zg)

	if math.sqrt((xg-xSizemm)^2+(ySizemm-yg)^2)<=radius then

		return 1,true


	elseif math.sqrt(xg^2+(ySizemm-yg)^2)<=radius then

		return 2,true


	elseif math.sqrt(xg^2+yg^2)<=radius then

		return 3,true


	elseif math.sqrt((xg-xSizemm)^2+yg^2)<=radius then

		return 4,true
	end
end,
 surface='fractional'
}

-- name of the PA containing electrostaticpotentials
pa:save('firstStage.pa#')
pa:refine{ convergence=1e-7}								-- , skipped_point = true

pa:fast_adjust{[1]=-quadVoltage, [2]=0, [3]=quadVoltage, [4]=0}


--- make PA containing Electric field magnitudes ---

pa2 = simion.pas:open()
pa2:load 'firstStage.pa0'

for xg, yg, zg in pa2:points() do
	if not pa2:electrode(xg,yg,zg) then

 		Ex,Ey,Ez = pa:field_vc(xg,yg,zg)					-- [V/gu]
		Emag = math.sqrt(Ex^2+Ey^2+Ez^2)*scaleFactor				-- E-field magnitude [V/mm] multiplied by 100 for scaling
 		pa2:potential(xg,yg,zg,Emag)
		end

	if pa2:electrode(xg,yg,zg) then
		Econ = 0								-- steady-state E-field inside a perfect conductor; can put in
		pa2:potential(xg,yg,zg,Econ)						-- skin depth and exp decay near the surface of a conductor
	end
end
pa2:save('firstStageEfield.pa0')

--- make PA containing Stark Potential ---

pa3 = simion.pas:open()
pa3:load 'firstStage.pa0'


for xg, yg, zg in pa3:points() do
	if not pa3:electrode(xg,yg,zg) then

 		Ex,Ey,Ez = pa:field_vc(xg,yg,zg)				-- [V/gu]
		Emag = math.sqrt(Ex^2+Ey^2+Ez^2)*scaleFactor			-- E-field magnitude [V/mm]
		W = starkLinearShiftSlope(starkSlope,Emag)	-- eV

		if xg==150 and yg==150 then

			print(W/wavenumber) --outputs actual value of stark energy at center of quadrupole in wavenumbers
			print(Emag)

		end

		pa3:potential(xg,yg,zg,W)


	end
end
pa3:save('firstStageStark.pa0')
