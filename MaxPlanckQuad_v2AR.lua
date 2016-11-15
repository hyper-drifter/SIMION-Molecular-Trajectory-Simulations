--This program should be used to simulate the quadrupole created by the Max-Planck institute for quantum optics
--Taylor Grubbs 3/19/16


local isIn2 = false --determines which state the quadrupole is in
local k = 0
local switch
local trapped = 0 --keeps track of whether or not a molecule was trapped

--this reads in velocities and positions from txt files. Velocites were generated randomly, weighted by a Maxwell-Boltzmann curve at 20mK. Positions are randomly uniform in a circle with radius .4 mm centered around 1.5,1.5 mm which is the center of the quadrupole

local g = io.open('Bigvelocities2Sim.tsv','r') --need to specify velocity and position files
local h = io.open('Bigvelocities2Sim.tsv','r')
local f = io.open('positions2Sim.tsv','r')
local i=1

while g:read('*line') ~= nil do --this initially finds how many particles will be tested so that SIMION will stop at the end of the list
      i = i+1
end
g:close()

local v = 1 --keeps track of which simulation is running

local dataFile = 'initParamData.tsv' --This file captures the initial parameters of particles that are trapped

local j = io.open(dataFile, 'w') --resets initParamData.tsv if it already exists in the directory
j:close()

j = io.open(dataFile, 'a') --records initial parameters of trapped molecules

local initPosx, initPosy, initVx, initVy --stores initial parameters

local quadrupole_frequency
local quadHalfPeriod



----------------------------------------Beginning of program----------------------------------------



simion.workbench_program()

--these variables should be adjustable in the workbench but they are not...

adjustable simTime = 10000 --Limits simulation time to 'x' us. Max Planck sim uses 10 ms or 10000 us
adjustable zSpeed = .001 --Velocity component in z direction. At around .001 the micro/macro-motion is very observable
adjustable frequencyStart = 5000--Allows the user to control different frequencies to be tested
adjustable frequencyStep = 1000
adjustable frequencyMax = 25000 --Sets upper limit to frequencies

if not quadrupole_frequency then --gives quadrupole initial starting frequency
      quadrupole_frequency = frequencyStart
      quadHalfPeriod = (1/(2*quadrupole_frequency))*10^6 --Have to convert the period from seconds to micro seconds. half because voltage changes every half period.
end



--This gives the molecules the velocities and positions imported from mathematica
function segment.initialize()

	ion_vx_mm = h:read('*number')
	ion_vy_mm = h:read('*number')
      ion_vz_mm = zSpeed

      ion_px_mm = f:read('*number')
      ion_py_mm = f:read('*number')


	initPosx = ion_px_mm
	initPosy = ion_py_mm
	initVx = ion_vx_mm
	initVy = ion_vy_mm
end
--This keeps track of the molecule's position and control the voltage switch
function segment.other_actions()

		switch = ion_time_of_flight - k*quadHalfPeriod
		--This switches the particle between the two PA instances; simulating the flipping of voltage.

		if switch / quadHalfPeriod >= 1 and isIn2 == false then

			ion_px_mm = ion_px_mm + 3 --workbench coordinates are in mm. The pa's themself have been scaled by a factor of 100

			isIn2 = true

			k = k+1

		elseif switch / quadHalfPeriod >= 1 and isIn2 == true then

			ion_px_mm = ion_px_mm - 3

			isIn2 = false

			k = k+1
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
            h = io.open('Bigvelocities2Sim.tsv','r')
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
