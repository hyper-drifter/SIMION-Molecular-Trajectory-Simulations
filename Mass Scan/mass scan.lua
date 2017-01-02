--This program is intended to investigate the use of Stark effect to filter molecules of different masses via time-varying electric fields.
--There is no way to change the voltages on the electrodes and then recalculate the new values of the stark effect in a fly program in SIMION (that I know of). So I have instead altered SIMION's charge parameter. Since we are assuming a Stark shift of alpha*Emag, the voltage and stark shifts are directly proportional. So changing the charge by a factor is like changing the voltage by the same factor.
--Taylor Grubbs 9-9-16


local isIn2 = false --determines which state the quadrupole is in
local halfCycleNum = 0
local timeAfterLastSwitch
local trapped = 0 --keeps track of whether or not a molecule was trapped

--this reads in velocities and positions from txt files.

local g = io.open('Bigvelocities2Sim.tsv','r') --need to specify velocity and position files
local velocityFile = io.open('Bigvelocities2Sim.tsv','r')
local positionFile = io.open('positions2Sim.tsv','r')
local massFile = io.open('masses2Sim.tsv','r')
local dipoleFile = io.open('moments2Sim.tsv', 'r')

local numberOfParticles=1

while g:read('*line') ~= nil do --this initially finds how many particles will be tested so that SIMION will stop at the end of the list
      numberOfParticles = numberOfParticles+1
end
g:close()

local currentParticleNum = 1 --keeps track of which simulation is running

local dataFile = 'initParamData.tsv' --This file captures the initial parameters of particles that are trapped

local outputFile = io.open(dataFile, 'w')
outputFile:close()
local outputFile = io.open(dataFile, 'a') --records initial parameters of trapped molecules. Append mode is better because the file is updated while the simulation is running

local initPosx, initPosy, initVx, initVy --stores initial parameters

local quadfrequency
local quadHalfPeriod
local constantMass = 20 --use this if you want to run with just a constant mass

simNum = 0 --keeps track of how many simulations will be ran

----------------------------------------Beginning of program----------------------------------------


simion.workbench_program()

--these variables should be adjustable in the workbench but they are not...

adjustable simTime = 10000 --Limits simulation time to 'x' us. Max Planck sim uses 10 ms or 10000 us
adjustable zSpeed = .001 --Velocity component in z direction. At around .001 the micro/macro-motion is very observable
adjustable frequencyStep = 1000
adjustable frequencyMax = 10000
adjustable startfrequency = 10000--allows the user to control different frequencies to be tested
adjustable simMax = 0 --maximum number of simulations to run if running simulation where voltage is variable
adjustable effectiveVoltage = 5000 --initial voltage. Keeps track of voltage that would be applied if "charge" was changed



if not quadfrequency then --intializes frequency and half period
	quadfrequency = startfrequency
	quadHalfPeriod = (1/(2*quadfrequency))*10^6 --Have to convert the period from seconds to micro seconds. half because voltage changes every half period.
end

--This gives the molecules the velocities and positions imported from mathematica
function segment.initialize()

      --velocities
	ion_vx_mm = velocityFile:read('*number')
	ion_vy_mm = velocityFile:read('*number')
      ion_vz_mm = zSpeed

      --positions
      ion_px_mm = positionFile:read('*number')
      ion_py_mm = positionFile:read('*number')

      --mass
      ion_mass = massFile:read('*number')

      --dipole moment
      ion_charge = dipoleFile:read('*number')

      --stores these values in global variables so they can be accessed at terminate segment
	initPosx = ion_px_mm
	initPosy = ion_py_mm
	initVx = ion_vx_mm
	initVy = ion_vy_mm

end



--This keeps track of the molecule's position and controls the voltage change
function segment.other_actions()

          --**simulates flipping of voltages**
		timeAfterLastSwitch = ion_time_of_flight - halfCycleNum*quadHalfPeriod

		--This switches the particle between the two PA instances; simulating the flipping of voltage.
		if timeAfterLastSwitch / quadHalfPeriod >= 1 and isIn2 == false then
			ion_px_mm = ion_px_mm + 3 --workbench coordinates are in mm. The pa's themself have been scaled by a factor of 100
			isIn2 = true
			halfCycleNum = halfCycleNum+1
		elseif timeAfterLastSwitch / quadHalfPeriod >= 1 and isIn2 == true then
			ion_px_mm = ion_px_mm - 3
			isIn2 = false
			halfCycleNum = halfCycleNum+1
		end

		--the following checks if the molecule reaches a distance 1 mm away from the center of the quadrupole and counts it as dead if so
		if ion_px_mm < 3 then
			if math.sqrt((1.5-ion_px_mm)^2 + (1.5-ion_py_mm)^2) >= 1 then
				ion_splat = 1
			end
		end
		if ion_px_mm > 3 then
			if math.sqrt((4.5-ion_px_mm)^2 + (1.5-ion_py_mm)^2) >= 1 then
				ion_splat = 1
			end
		end

		--this keeps the ion inside the simulation region if you add a z component to velocity
		if ion_pz_mm >= 2.99 then
			ion_pz_mm = -2.99
		end

		--Stops simulation after a certain amount of time.
		if ion_time_of_flight >= simTime then
			ion_splat = 2
			trapped = 1
		end
end



--ensures that the time step smaller than the period of the quadrupole
function segment.tstep_adjust()
	ion_time_step = math.min(ion_time_step, .01*quadHalfPeriod)
end



--grabs initial parameters if molecule is trapped and stops simulation when velocity/positions lists are empty
function segment.terminate()

	if trapped == 1 then
		outputFile:write(quadfrequency..'\t'..ion_mass..'\t'..ion_charge..'\t'..initPosx..'\t'..initVx..'\t'..initPosy..'\t'..initVy..'\n')
		outputFile:flush()--updates file while simulation is running
		trapped = 0
	end

	halfCycleNum = 0 --resets the timeAfterLastSwitch variable for next molecule

	isIn2 = false -- resets starting point

	currentParticleNum = currentParticleNum+1 --moves to next particle

	if currentParticleNum == numberOfParticles then --this keeps SIMION from throwing an error at the end
            currentParticleNum = 1
            velocityFile:close()
            positionFile:close()
            massFile:close()
            dipoleFile:close()

            velocityFile = io.open('Bigvelocities2Sim.tsv','r')
            positionFile = io.open('positions2Sim.tsv','r')
            massFile = io.open('masses2Sim.tsv','r')
            dipoleFile = io.open('moments2Sim.tsv','r')

            quadfrequency = quadfrequency + frequencyStep --updates frequency
		quadHalfPeriod = (1/(2*quadfrequency))*10^6

	end

      sim_rerun_flym = 1

      --specifies when to stop simulation
      if quadfrequency > frequencyMax then
            print('I am done')
            sim_rerun_flym = 0
            velocityFile:close()
            positionFile:close()
            outputFile:close()
      end

end
