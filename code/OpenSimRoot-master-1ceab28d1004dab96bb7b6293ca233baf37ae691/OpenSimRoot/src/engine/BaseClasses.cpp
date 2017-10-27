/*
Copyright © 2016, The Pennsylvania State University
All rights reserved.

Copyright © 2016 Forschungszentrum Jülich GmbH
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted under the GNU General Public License v3 and provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

Disclaimer
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You should have received the GNU GENERAL PUBLIC LICENSE v3 with this file in license.txt but can also be found at http://www.gnu.org/licenses/gpl-3.0.en.html

NOTE: The GPL.v3 license requires that all derivative work is distributed under the same license. That means that if you use this source code in any other program, you can only distribute that program with the full source code included and licensed under a GPL license.

 */

#include "BaseClasses.hpp"
#include "../cli/Messages.hpp"
#include "../math/MathLibrary.hpp"
#include "../math/InterpolationLibrary.hpp"
#include "Origin.hpp"
#if _WIN32 || _WIN64
#include <algorithm>
#endif

DerivativeBase::DerivativeBase(SimulaDynamic* const new_pSD) :
		pSD(new_pSD)
{
}

DerivativeBase::~DerivativeBase() {
}

void DerivativeBase::calculate(const Time &t, double &var) {
	msg::error(
			"DerivativeBase::calculate: basefunction called. Object name: "
					+ pSD->getName());
}

void DerivativeBase::calculate(const Time &t, const Coordinate &pos, double &var) {
	calculate(t,var);
}

void DerivativeBase::calculate(const Time &t, int &var) {
	msg::error("DerivativeBase::calculate: basefunction called");
}

void DerivativeBase::calculate(const Time &t, std::string &var) {
	msg::error("DerivativeBase::calculate: basefunction called");
}

void DerivativeBase::calculate(const Time &t, Coordinate &var) {
	msg::error("DerivativeBase::calculate: basefunction called");
}

bool DerivativeBase::postIntegrationCorrection(SimulaVariable::Table &table) {
	return false;
}

bool DerivativeBase::postIntegrationCorrection(SimulaPoint::Table &table) {
	return false;
}

void DerivativeBase::addObject(SimulaBase *newObject){
	msg::error("DericatiBase::addObject called, but not implemented for class "+getName());
}

SimulaBase* DerivativeBase::getNext(const Time & t) {
	msg::error(
			"DerivativeBase: no getNext function implemented for this class with name "
					+ getName());
	return nullptr;
}

/* Instead of locking - we could try to switch to an estimate mode - which would allow for a more iterative method to resolve
 * note that when using such methods user much be sure that the 'estimate' tag is preserved throughout the calculations
 * note that we can still try to detect loop by counting how many estimates are asked for - this number should be small (4 for rungkuta, 1 for heuns etc).
 */

std::string DerivativeBase::getName() const {
	msg::error(
			"DerivativeBase::getName: programming error: forgot to include a getName() function in class used by "
					+ pSD->getName());
	return "";
}


SimplePredictor::SimplePredictor(SimulaVariable::Table & data,
		const Time& deltaT) {
	//last state rate couple
	SimulaVariable::Table::reverse_iterator it(data.rbegin());
	s1 = (it->second);
	t1 = (it->first);

	//last state rate couple
	if (data.size() > 1)
		++it;
	s0 = (it->second);
	t0 = (it->first);

	//check estimate consistency
	if (data.size() > 1 && s0.rate.estimate && !s1.estimate) {
		msg::warning("SimplePredictor:: fixing estimate in s0");
		s0.rate.estimate = false;
		it->second.rate.estimate = false;
	}

	//set predictor
	t2 = t1 + deltaT;
	newTime=t2;
	s2 = s1;
	if (data.size() > 1) {
		linearInterpolation(t0, t1, t2, s0.rate, s1.rate, s2.rate);
		linearInterpolation(t0, t1, t2, s0.state, s1.state, s2.state);
		//second order correction
		s2.state += (s2.rate - s1.rate) / 2;
	}

	//set step
	step = 1;
}
SimplePredictor::SimplePredictor(SimulaPoint::Table & data,
		const Time& deltaT) {
	//last state rate couple
	SimulaPoint::Table::reverse_iterator it(data.rbegin());
	p1 = (it->second);
	t1 = (it->first);

	//second last state rate couple
	if (data.size() > 1)
		++it;
	p0 = (it->second);
	t0 = (it->first);

	//check estimate consistency
	if (data.size() > 1 && p0.rate.estimate && !p1.estimate) {
		msg::warning("SimplePredictor:: fixing estimate in p0",2);
		s0.rate.estimate = false;
		it->second.rate.estimate = false;
	}

	//set predictor
	t2 = t1 + deltaT;
	newTime=t2;
	p2 = p1;
	if (data.size() > 1) {
		//first order
		linearInterpolation(t0, t1, t2, p0.rate, p1.rate, p2.rate, true);
		linearInterpolation(t0, t1, t2, p0.state, p1.state, p2.state);
	}

	//set step
	step = 1;
}

IntegrationBase::IntegrationBase() :
		pSD(nullptr), predictor(nullptr), multiplier(nullptr), rateMultiplier(nullptr) {
}

IntegrationBase::~IntegrationBase() {
}




void IntegrationBase::initiate(SimulaTimeDriven* pSP,
		DerivativeBase & rateCalculator) {
	if (!pSD) {
		pSD = pSP;
		//set defaultSpatialIntegrationLength
		if (defaultSpatialIntegrationLength < 0) {
			//read default precision for spatial integration - this is a parameter in cm
			SimulaBase* o(
					ORIGIN->getChild("simulationControls")->existingChild(
							"integrationParameters"));
			if(o) o = o->existingChild("defaultSpatialIntegrationLength");
			if (o) {
				o->get(defaultSpatialIntegrationLength);
				if (o->getUnit().order != Unit("cm").order)
					msg::error(
							"IntegrationBase: Unit for IntegrationParameters/defaultSpatialIntegrationLength must be of order length i.e. use cm or mm etc");
				//check this is unit is of order length
				//convert to meters
				defaultSpatialIntegrationLength *= o->getUnit().factor;
				//convert to cm
				defaultSpatialIntegrationLength /= Unit("cm").factor;
			}
			if (defaultSpatialIntegrationLength <= 0) {
				defaultSpatialIntegrationLength = 1;
			}
			msg::warning("IntegrationBase:: setting defaultSpatialIntegrationLength to "+std::to_string(defaultSpatialIntegrationLength)+".",1);
		}
	}
	if(pSP->existingChild("multiplier")){
		multiplier = pSP->getChild("multiplier");
	}
	if(pSP->existingChild("rateMultiplier")){
		rateMultiplier = pSP->getChild("rateMultiplier");
	}
}


void IntegrationBase::integrate(SimulaVariable::Table & data,
		DerivativeBase & rateCalculator) {
	msg::error(
			"IntegrationBase:: integrate (SimulaVariable) is not implemented, please provide implementation for "
					+ this->getName()
					+ " or use different integration method for "
					+ pSD->getPath());
}

void IntegrationBase::integrate(SimulaPoint::Table & data,
		DerivativeBase & rateCalculator) {
	msg::error(
			"IntegrationBase:: integrate (SimulaPoint) is not implemented, please provide implementation for "
					+ this->getName()
					+ " or use different integration method for "
					+ pSD->getPath());
}

void IntegrationBase::predict(const Time&t, State &ps) {//, SimulaVariable::Table & data){
		/*	bool dp=false;
		 if(!predictor) {
		 Time lastTime=data.rbegin()->first, newTime, deltaT;
		 timeStep(lastTime, newTime, deltaT);
		 predictor=new SimplePredictor(data,deltaT);
		 dp=true;
		 }
		 */
	//timestep
	Time deltaT = t - predictor->t1;
	//integration
	if (predictor->step == 1) {				//forward
		ps = predictor->s1.state + predictor->s1.rate * deltaT;
	} else {				//backward
		ps = predictor->s1.state + predictor->s2.rate * deltaT;
	}
	//estimate is true
	ps.estimate = true;

	/*	 if(dp){
	 delete predictor;
	 predictor=nullptr;
	 }*/
}

void IntegrationBase::predictRate(const Time&t, State &ps) {//, SimulaVariable::Table & data){
		/*	bool dp=false;
		 if(!predictor) {
		 Time lastTime=data.rbegin()->first, newTime, deltaT;
		 timeStep(lastTime, newTime, deltaT);
		 predictor=new SimplePredictor(data,deltaT);
		 dp=true;
		 }*/

	//interpolate the rate
	if (predictor->step == 1) {				//forward
		if (predictor->t0 == predictor->t1) {
			ps = predictor->s1.rate;
		} else {
			//linear extrapolation based on the change in the previous timestep
			linearInterpolation(predictor->t0, predictor->t1, t,
					predictor->s0.rate, predictor->s1.rate, ps);
			//avoid sign changes
			if (predictor->s0.rate >= 0 && predictor->s1.rate >= 0 && ps < 0)
				ps.value = 0;
			if (predictor->s0.rate <= 0 && predictor->s1.rate <= 0 && ps > 0)
				ps.value = 0;
		}
	} else {
		//linear interpolation based on the estimated rate in s2
		linearInterpolation(predictor->t1, predictor->t2, t, predictor->s1.rate,
				predictor->s2.rate, ps);
	}
	//estimate is true
	ps.estimate = true;

	/*	 if(dp){
	 delete predictor;
	 predictor=nullptr;
	 }*/
}

void IntegrationBase::predict(const Time&t, Coordinate& pp) {//, SimulaPoint::Table & data){
		/*	booTime IntegrationBase::stepLimiter(-1);
		 * l dp=false;
		 if(!predictor) {
		 Time lastTime=data.rbegin()->first, newTime, deltaT;
		 timeStep(lastTime, newTime, deltaT);
		 predictor=new SimplePredictor(data,deltaT);
		 dp=true;
		 predictor->step=1;
		 }*/

	//delta t
	Time deltaT = t - predictor->t1;
	//integration
	if (predictor->step == 1) {				//forward
		pp = predictor->p1 + predictor->p1.rate * deltaT;
	} else {				//backward
		pp = predictor->p1 + predictor->p2.rate * deltaT;
	}
	//estimate is true, for this is an prediction
	pp.estimate = true;

	/*if(dp){
	 delete predictor;
	 predictor=nullptr;
	 }*/
}

void IntegrationBase::predictRate(const Time&t, Coordinate&pp) {//, SimulaPoint::Table & data){
		/*	bool dp=false;
		 if(!predictor) {
		 Time lastTime=data.rbegin()->first, newTime, deltaT;
		 timeStep(lastTime, newTime, deltaT);
		 predictor=new SimplePredictor(data,deltaT);
		 dp=true;
		 }*/

	//interpolate the rate
	if (predictor->step == 1) {				//forward
		if (predictor->t0 == predictor->t1) {
			pp = predictor->p1.rate;			///@todo does this ever happen?
		} else {
			//linear extrapolation based on the change in the previous timestep
			linearInterpolation(predictor->t0, predictor->t1, t,
					predictor->p0.rate, predictor->p1.rate, pp, true);
		}
	} else {
		//linear interpolation based on the estimated rate in s2
		linearInterpolation(predictor->t1,
				predictor->t2 / (Time) predictor->step, t, predictor->p1.rate,
				predictor->p2.rate, pp, true);
	}
	//estimate is true
	pp.estimate = true;

	/*	 if(dp){
	 delete predictor;
	 predictor=nullptr;
	 }*/
}

void IntegrationBase::getTimeStepInfo(Time & st, Time & rt) const {
	if (predictor) {
		///todo this is not safe as the predictor could easily have another implementation.
		st = predictor->t1;
		rt = predictor->newTime;
	} else {
		st = -1;
		rt = -1;
	}
}

void IntegrationBase::timeStep(const Time & lastTime, Time & newTime,
		Time & deltaT, const double & l) {
	//determine timestep deltaT
	/*if(pSD->preferedTimeStep()>pSD->minTimeStep()){
	 //use the prefered time step
	 deltaT=pSD->preferedTimeStep();
	 }else if(l>0 && defaultSpatialIntegrationLength>0){
	 //scale in reference to l
	 deltaT=defaultSpatialIntegrationLength.value/l;
	 }else{
	 //use timestepping module (first time use minimum timestep)
	 timeSteppingModule(deltaT);
	 }*/
	deltaT = pSD->preferedTimeStep();
	///todo a variable timestep causes problems
	//if(deltaT<TIMEERROR) deltaT=pSD->lastTimeStep()*1.05;

	//adjust to spatial resolution
	if (l > 1e-5) {
		//scale in reference to l
		deltaT = defaultSpatialIntegrationLength / l;
	}

	//maxmimum and minimum timesteps
	deltaT = std::max(pSD->minTimeStep(), deltaT);
	deltaT = std::min(pSD->maxTimeStep(), deltaT);

	//newTime
	newTime = lastTime + deltaT;

	//round off timesteps to multiple of min time step, so that processes with same minTimeStep are in sync
	/*	if(pSD->syncTimeStep()){
	 newTime=pSD->minTimeStep()*floor((newTime+TIMEERROR)/pSD->minTimeStep());
	 deltaT=newTime-lastTime;
	 if(deltaT<TIMEERROR) {
	 newTime+=pSD->preferedTimeStep();
	 deltaT+=pSD->preferedTimeStep();
	 }
	 }else{*/
	Time synctime = MINTIMESTEP * floor((newTime + TIMEERROR) / MINTIMESTEP);
	if (synctime - lastTime > TIMEERROR) {
		//snapping the timestep to the default rhytm at which the other processes synchronize.
		newTime = synctime;
		deltaT = newTime - lastTime;
	} else {
		//trying to balance the timestep so that the timestep snapping does not cause large oscillations in the timestep.
		//note this is possible if the process uses timesteps smaller than the default MINTIMESTEP, that is
		//MINTIMESTEP>pSD->getMinTimeStep() and we want to reduce the timestep so that multiple steps will
		//cause synchronization with MINTIMESTEP.
		double steps = ceil(MINTIMESTEP / deltaT);
		deltaT = MINTIMESTEP / steps;
		newTime = lastTime + deltaT;
	}
	//   }

	//do not run ahead of objects that are calling this.
	Time wt = pSD->getWallTime(newTime);
	newTime = wt;
	deltaT = newTime - lastTime;
	//just for debuggin info
	/*if(deltaT+TIMEERROR<MAXTIMESTEP) {
	 msg::warning("Using smaller timestep");
	 }else{
	 msg::warning("Using larger timestep");
	 }*/
}

void IntegrationBase::timeSteppingModule(Time & deltaT) {
	msg::error(
			"IntegrateBase::timeSteppingModule: timeSteppingModule not implemented for this integration method ("
					+ getName()
					+ "). Make sure preferedTimeStep>minTimeStep or use different integration method for "
					+ pSD->getPath());
}

double IntegrationBase::defaultSpatialIntegrationLength(-1);

std::map<std::string, derivativeBaseInstantiationFunction> & BaseClassesMap::getDerivativeBaseClasses() {
	static std::map<std::string, derivativeBaseInstantiationFunction> derivativeBaseClasses;
	return derivativeBaseClasses;
}
std::map<std::string, objectGeneratorInstantiationFunction> & BaseClassesMap::getObjectGeneratorClasses() {
	static std::map<std::string, objectGeneratorInstantiationFunction> objectGeneratorClasses;
	return objectGeneratorClasses;
}
std::map<std::string, integrationInstantiationFunction> & BaseClassesMap::getIntegrationClasses() {
	static std::map<std::string, integrationInstantiationFunction> integrationClasses;
	return integrationClasses;
}
