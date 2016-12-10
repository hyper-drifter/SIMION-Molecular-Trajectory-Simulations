--This program should be used to simulate the quadrupole created by the Max-Planck institute for quantum optics
--Taylor Grubbs 3/19/16


local isIn2 = false --determines which state the quadrupole is in
local k = 0
local switch
local trapped = 0 --keeps track of whether or not a molecule was trapped

--this reads in velocities and positions from txt files.

local g = io.open('velocities2Sim.tsv','r') --need to specify velocity and position files
local h = io.open('velocities2Sim.tsv','r')
local f = io.open('positions2Sim.tsv','r')
local i=1

while g:read('*line') ~= nil do --this initially finds how many particles will be tested so that SIMION will stop at the end of the list
      i = i+1
end
g:close()

local v = 1 --keeps track of which simulation is running

local dataFile = 'initParamData.tsv' --This file captures the initial parameters of particles that are trapped

local j = io.open(dataFile, 'w') --resets timeAlive.tsv if it already exists in the directory
j:close()

j = io.open(dataFile, 'a')

local initPosx, initPosy, initVx, initVy --stores initial parameters

local quadrupole_frequency
local quadHalfPeriod

local scaleFactor = 1. --Use this number to scale down to micrometer dimensions



----------------------------------------Beginning of program----------------------------------------



simion.workbench_program()

--these variables should be adjustable in the workbench but they are not...

adjustable simTime = 200 --Limits simulation time to <x> us. Max Planck sim uses 10 ms or 10000 us
adjustable zSpeed = .001 --Velocity component in z direction. At around .001 the micro/macro-motion is very observable
adjustable frequencyStart = 100000.-- Hz, Allows the user to control different frequencies to be tested
adjustable frequencyStep = 1000
adjustable frequencyMax = 500000.--Sets upper limit to frequencies
adjustable maxSpeed = .015 --Maximum allowed speed

if not quadrupole_frequency then --gives quadrupole initial starting frequency
      quadrupole_frequency = frequencyStart
      quadHalfPeriod = (1/(2*quadrupole_frequency))*10^6 --Have to convert the period from seconds to micro seconds. half because voltage changes every half period.
end



--This gives the molecules the velocities and positions imported from mathematica
function segment.initialize()

	ion_vx_mm = h:read('*number')
	ion_vy_mm = h:read('*number')
      ion_vz_mm = zSpeed

      ion_px_mm = f:read('*number')/scaleFactor
      ion_py_mm = f:read('*number')/scaleFactor
      ion_pz_mm = 0

      --Throws out any particles that are moving too fast
      if math.abs(ion_vx_mm) >= maxSpeed or math.abs(ion_vy_mm) >= maxSpeed then
            ion_splat = 1
      end

	initPosx = ion_px_mm
	initPosy = ion_py_mm
	initVx = ion_vx_mm
	initVy = ion_vy_mm
end
--This keeps track of the molecule's position and control the voltage switch
function segment.other_actions()

		switch = ion_time_of_flight - k*quadHalfPeriod
		--This switches the particle between the two PA instances; simulating the flipping of voltage.

		if switch / quadHalfPeriod >= 1. then

                  if isIn2 == false then

			      ion_px_mm = ion_px_mm + .05/scaleFactor --workbench coordinates are in mm. The pa's themself have been scaled by a factor of 100

			      isIn2 = true

			      k = k+1

		      elseif isIn2 == true then

			      ion_px_mm = ion_px_mm - .05/scaleFactor

			      isIn2 = false

			      k = k+1
                  end
		end

		--the following checks if the molecule gets too close to an electrode or exits the trap entirely. It then counts the molecule as terminated

            if ion_py_mm <= .0005/scaleFactor or ion_py_mm >= .0495/scaleFactor then
                  ion_splat = 1
            end

		if isIn2 == false then
                  if ion_px_mm <= .0005/scaleFactor or ion_px_mm >= 0.0495/scaleFactor then
				ion_splat = 1
			end
		end

		if isIn2 == true then
                  if ion_px_mm <= .0505/scaleFactor or ion_px_mm >= 0.0995/scaleFactor then
				ion_splat = 1
			end
		end


		--this keeps the ion inside the simulation region if you add a z component to velocity
		if ion_pz_mm >= .049/scaleFactor then
			ion_pz_mm = -.049/scaleFactor
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
--grabs initial paramters if molecule is trapped and stops simulation when velocity/positions lists are empty
function segment.terminate()

      if v == 1 then
            j:write(quadrupole_frequency..'\n') --writes quadrupole frequency to trapped.tsv in order to differentiate data from different frequencies
      end

	if trapped == 1 then
		j:write(initPosx..'	'..initVx..'	'..initPosy..'	'..initVy..'\n')
		trapped = 0
	end


	k = 0 --resets the switch variable for next molecule

	isIn2 = false -- resets starting point

	v = v+1 --moves to next particle


	if v == i then --this keeps SIMION from throwing an error at the end
            v = 1
            quadrupole_frequency = quadrupole_frequency+frequencyStep
            print(quadrupole_frequency)
            quadHalfPeriod = (1/(2*quadrupole_frequency))*10^6
            h:close()
            f:close()
            h = io.open('velocities2Sim.tsv','r')
            f = io.open('positions2Sim.tsv','r')
	end

      sim_rerun_flym = 1

      if quadrupole_frequency > frequencyMax then
            print('I am done')
            sim_rerun_flym = 0
            h:close()
            f:close()
            j:close()
      end
end
