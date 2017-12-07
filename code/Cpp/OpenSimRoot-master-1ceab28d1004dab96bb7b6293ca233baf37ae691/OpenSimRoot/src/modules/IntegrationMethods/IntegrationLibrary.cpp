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

#include "IntegrationLibrary.hpp"
#include "../../math/InterpolationLibrary.hpp"
#include "../../engine/SimulaVariable.hpp"
#include "../../engine/SimulaPoint.hpp"
#include "../../cli/Messages.hpp"
#include "../../tools/StringExtensions.hpp"
#include "../../math/VectorMath.hpp"
#include <math.h>
#include "../../engine/Origin.hpp"
#include "../../math/MathLibrary.hpp"


//======================Integration functions===================================
//Integration using forward euler
ForwardEuler::ForwardEuler():
	IntegrationBase(),totError(0),totChange(0),updateRate(0),callCount(0){
}
void ForwardEuler::integrate(SimulaVariable::Table & data, DerivativeBase &rateCalculator){
	//name (stored for pSD is destructed before this one and we want to use the name in the destructor)
	name=pSD->getName(); //can not be done during construction.
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	State dummy(predictor->s1.rate);
	if(dummy.estimate) {
		rateCalculator.calculate(predictor->t1,dummy.value);
		if (rateMultiplier){
			double temp;
			rateMultiplier->get(predictor->t1, temp);
			dummy.value *= temp;
		}
		//check estimate & clear locks
		pSD->keepEstimate(dummy.estimate);
		//insert updated predictor
		predictor->s1.rate=dummy;
		data.rbegin()->second.rate=dummy;
	}

	//forward euler
	StateRate newValue;
	newValue.state=predictor->s1.state+dummy*deltaT;
	predictRate(newTime,newValue.rate);
	data.insert(std::pair<Time,StateRate>(newTime,newValue));

	//post integration correction
	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		//use forward euler now
		SimulaVariable::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->s1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	data.rbegin()->second.estimate=newValue.estimate;
	data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			State dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy.value);
			if (rateMultiplier){
				double temp;
				rateMultiplier->get(newTime, temp);
				dummy.value *= temp;
			}
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}

	delete predictor;
	predictor=nullptr;
}

void ForwardEuler::integrate(SimulaPoint::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	State l(0);
	if(defaultSpatialIntegrationLength>0) l=data.rbegin()->second.rate.length();
	timeStep(lastTime, newTime, deltaT, l.value);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	Coordinate dummy(predictor->p1.rate);
	if(dummy.estimate) {
		rateCalculator.calculate(predictor->t1,dummy);
		//check estimate & clear locks
		pSD->keepEstimate(dummy.estimate);
		//insert updated predictor
		predictor->p1.rate=dummy;
		data.rbegin()->second.rate=predictor->p1.rate;
	}

	//forward euler
	MovingCoordinate newValue;
	newValue=predictor->p1+predictor->p1.rate*deltaT;
	predictRate(newTime,newValue.rate);
	data.insert(std::pair<Time,MovingCoordinate>(newTime,newValue));

	//post integration correction
	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		SimulaPoint::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->p1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	data.rbegin()->second.estimate=newValue.estimate;
	data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			Coordinate dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy);
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}

	delete predictor;
	predictor=nullptr;
}

//numerical error estimate
///@todo this needs to be updated and inserted into timeSteppingModule()
//eulers error=1/2 * dt^2 * |y''| = 0.5*dt*dt*(rate(t2)-rate(t1))/dt
/*double err(0);
if(!last.rate.estimate && data.size()>1) err=((last.rate.value-last2.rate.value)*(lastTime-lastTime2)*0.5);
totError+=err;
err=fabs(err);
double ref(fabs(last.value-last2.value));
totChange+=ref;*/
//numerical error estimate for simulaPoint
//double ref(vectorlength(last.state-last2.state).value);
//totChange+=ref;

//error based on length, not direction
//double currentError(((vectorlength(last.rate)-vectorlength(last2.rate))*(lastTime-lastTime2)*0.5).value);
//if(!updateRate && !last.rate.estimate && ref>0.1) totError+=currentError;
//double err(vectorlength(currentError));
//double l(last.rate.length().value);

ForwardEuler::~ForwardEuler(){
	//warn about total errors
	if(callCount>20){//only report for objects that existed for significant amount of time
		if(1e2*totError>totChange)
			msg::warning("ForwardEuler: results may not be accurate for "+name);
		if(1e8*totError<totChange)
			msg::warning("ForwardEuler: you could probably use a lager (minimum) time step for "+name);
	}
}
std::string ForwardEuler::getName()const{
	return "ForwardEuler";
}



BackwardEuler::BackwardEuler():
	IntegrationBase(){}
void BackwardEuler::integrate(SimulaVariable::Table & data, DerivativeBase &rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	State dummy(predictor->s1.rate);
	if(dummy.estimate) {
		rateCalculator.calculate(predictor->t1,dummy.value);
		if (rateMultiplier){
			double temp;
			rateMultiplier->get(predictor->t1, temp);
			dummy.value *= temp;
		}
		//check estimate & clear locks
		pSD->keepEstimate(dummy.estimate);
		//insert updated predictor
		predictor->s1.rate=dummy;
		data.rbegin()->second.rate=dummy;
	}

	//update rate
	rateCalculator.calculate(newTime,dummy.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(newTime, temp);
		dummy.value *= temp;
	}
	//check estimate & clear locks
	pSD->keepEstimate(dummy.estimate);

	//backward euler
	StateRate newValue;
	newValue.state=predictor->s1.state+dummy*deltaT;
	newValue.rate=dummy;
	data.insert(std::pair<Time,StateRate>(newTime,newValue));

	//post integration correction
	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		//use forward euler now
		SimulaVariable::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->s1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	data.rbegin()->second.estimate=newValue.estimate;
	data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			State dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy.value);
			if (rateMultiplier){
				double temp;
				rateMultiplier->get(newTime, temp);
				dummy.value *= temp;
			}
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}

	delete predictor;
	predictor=nullptr;
}
std::string BackwardEuler::getName()const{
	return "BackwardEuler";
}



//heun's method - this is the average of a forward and backward euler.
Heuns::Heuns():IntegrationBase(){}
void Heuns::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	State r1(predictor->s1.rate);
	if(r1.estimate) {
		rateCalculator.calculate(predictor->t1,r1.value);
		if (rateMultiplier){
			double temp;
			rateMultiplier->get(predictor->t1, temp);
			r1.value *= temp;
		}
		//it's still possible to get an estimate here when time steps are not synchronized
		pSD->keepEstimate(r1.estimate);
		//insert updated predictor
		predictor->s1.rate=r1;
		data.rbegin()->second.rate=r1;
	}

	//backward
	State r4;
	rateCalculator.calculate(newTime,r4.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(newTime, temp);
		r4.value *= temp;
	}
	pSD->keepEstimate(r4.estimate);
	predictor->s2.rate=r4;
	predictor->t2=newTime;

	//heuns: average of forward and backward euler
	StateRate newValue(predictor->s1);
	newValue.state+=(r1+r4)*deltaT/2;
	newValue.rate=r4;

	data.insert(std::pair<Time,StateRate>(predictor->t2,newValue));

	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		//use forward euler now
		SimulaVariable::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->s1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	if(newValue.estimate) data.rbegin()->second.estimate=newValue.estimate;
	if(newValue.rate.estimate) data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			State dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy.value);
			if (rateMultiplier){
				double temp;
				rateMultiplier->get(newTime, temp);
				dummy.value *= temp;
			}
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}

	delete predictor;
	predictor=nullptr;
}

void Heuns::integrate(SimulaPoint::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	State l(0);
	if(defaultSpatialIntegrationLength>0) l=data.rbegin()->second.rate.length();
	timeStep(lastTime, newTime, deltaT, l.value);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	MovingCoordinate r1(predictor->p1.rate);
	if(r1.estimate) {
		rateCalculator.calculate(predictor->t1,r1);
		//it's still possible to get an estimate here when time steps are not synchronized
		pSD->keepEstimate(r1.estimate);
		//insert updated predictor
		predictor->p1.rate=r1;
		data.rbegin()->second.rate=r1;
	}

	//backward
	MovingCoordinate r4;
	rateCalculator.calculate(newTime,r4);
	pSD->keepEstimate(r4.estimate);
	predictor->p2.rate=r4;
	predictor->t2=newTime;

	//heuns: average of forward and backward euler
	MovingCoordinate newValue(predictor->p1);
	l=(r4.length()+r1.length())/2;
	newValue.rate=(r1+r4)*deltaT/2;
	newValue.rate.setLength(l);
	newValue.state+=newValue.rate*deltaT;
	newValue.rate=r4;

	data.insert(std::pair<Time,MovingCoordinate>(predictor->t2,newValue));

	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		SimulaPoint::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->p1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	if(newValue.estimate) data.rbegin()->second.estimate=newValue.estimate;
	if(newValue.rate.estimate) data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			Coordinate dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy);
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}

	delete predictor;
	predictor=nullptr;
}
std::string Heuns::getName()const{
	return "Heuns";
}

//HeunsII
///@todo numerical error estimate: HeunsII error comes from the assumption that second derivative is linear ( not constant as in forward euler)
HeunsII::HeunsII():IntegrationBase(){}
void HeunsII::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	msg::error("HeunsII: heunsII is used but is known to cause oscillations. Use a different method for "+pSD->getPath());
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//calculate rate
	Time ht(lastTime+deltaT*0.5);
	State r1;
	rateCalculator.calculate(ht,r1.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(ht, temp);
		r1.value *= temp;
	}
	pSD->keepEstimate(r1.estimate);

	//fix estimated rate
	//update rate
	State r0;
	linearInterpolation(predictor->t0,ht,predictor->t1,predictor->s0.rate,r1,r0);
	r0.estimate=r1.estimate;
	data.rbegin()->second.rate=r0;

	//insert half time
	StateRate newValue(predictor->s1);
	newValue.state+=r1*ht;
	newValue.rate=r1;
	data.insert(std::pair<Time,StateRate>(ht,newValue));

	//integrate
	newValue=(predictor->s1);
	newValue.state+=r1*deltaT;
	newValue.rate=r1;//*2-r0; //this can lead to oscillations and or negative rates, but anything else will cause wrong rate interpolation.
	newValue.rate.estimate=true;

	//update
	data.insert(std::pair<Time,StateRate>(newTime,newValue));

	//post-integration correction
	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		//use forward euler now
		SimulaVariable::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->s1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	data.rbegin()->second.estimate=newValue.estimate;
	data.rbegin()->second.rate.estimate=newValue.estimate;//heuns 11 the rate's estimate does not get updated again

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			State dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy.value);
			if (rateMultiplier){
				double temp;
				rateMultiplier->get(newTime, temp);
				dummy.value *= temp;
			}
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}

	//delete predictor
	delete predictor;
	predictor=nullptr;
}


std::string HeunsII::getName()const{
	return "HeunsII";
}

//Runge-kutta 4
///@todo - it would be better if we used a higher order interpolation method with this.
RungeKutta4::RungeKutta4():IntegrationBase(){}
void RungeKutta4::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	State r1(predictor->s1.rate);
	if(r1.estimate) {
		rateCalculator.calculate(predictor->t1,r1.value);
		if (rateMultiplier){
			double temp;
			rateMultiplier->get(predictor->t1, temp);
			r1.value *= temp;
		}
		//it's still possible to get an estimate here when time steps are not synchronized
		pSD->keepEstimate(r1.estimate);
		//insert updated predictor
		predictor->s1.rate=r1;
		data.rbegin()->second.rate=r1;
	}

	//half time
	Time ht(lastTime+deltaT/2);

	//step 2
	State r2;
	rateCalculator.calculate(ht,r2.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(ht, temp);
		r2.value *= temp;
	}
	pSD->keepEstimate(r2.estimate);
	predictor->s2.rate=r2;
	predictor->step=2;
	predictor->t2=ht;//for extrapolation of predicted rates - but also reduces the timestep of called models

	//step 3
	State r3;
	rateCalculator.calculate(ht,r3.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(ht, temp);
		r3.value *= temp;
	}
	pSD->keepEstimate(r3.estimate);
	predictor->s2.rate=r3;

	//step 4
	State r4;
	rateCalculator.calculate(newTime,r4.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(newTime, temp);
		r4.value *= temp;
	}
	pSD->keepEstimate(r4.estimate);
	predictor->s2.rate=r4;
	predictor->t2=newTime;


	//rungekutta 4
	StateRate newValue(predictor->s1);
	newValue.rate=(r1+((r2+r3)*2)+r4)/6;
	newValue.state+=newValue.rate*deltaT;
	newValue.rate=r4;



	//check stability
	///todo implement correct checks for the integration error.
/*	if((r2+r3)!=0){
		State rd((r2-r3)/(r2+r3));
		if(fabs(rd.value) > 0.05 ){
			msg::warning("RungeKutta4:: large (>5%) difference between forward and backward euler, results may be inaccurate for"+pSD->getName());
		}
	}

	//check stability
	if(!newValue.estimate && newTime>1+pSD->getStartTime() && fabs(r1.value)>1e-4) {
		double rd1(deltaT*(r4-r1).value);
		double rd2(deltaT*(((r2+r3)/2)-r1).value);
		double rdm(maximum(rd1,rd2));
		if( rdm/r1.value > 0.01){ //the rate is changing with at least 1 percent during the timestep
			if(fabs((rd1-rd2)/rdm)> 0.1 ){ //
				msg::warning("RungeKutta: steep second derivative detected for "+pSD->getName()+". Reduce timestep?");
			}
		}
	}
*/

	data.insert(std::pair<Time,StateRate>(newTime,newValue));

	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		//use forward euler now
		SimulaVariable::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->s1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	if(newValue.estimate) data.rbegin()->second.estimate=newValue.estimate;
	if(newValue.rate.estimate) data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate){// && pSD->numberOfDependents()){ ///@todo somehow this seems right but still causes objects to trail behind and thereby cause late updates the branches (for estimates do not create branches) and than an error.
			//allow possible dependent object to clear their estimates
			State dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy.value);
			if (rateMultiplier){
				double temp;
				rateMultiplier->get(newTime, temp);
				dummy.value *= temp;
			}

/*	note that the estimate might not clear when is produced by a forward stepping, in forexample SimulaDerivative<>::getRate
 * 	if(dummy.estimate) {
				pSD->keepEstimate(dummy.estimate);//unlock, shouldn't be necessary?
				rateCalculator.calculate(rt,dummy);
			}
			if(dummy.estimate){
				msg::warning("RungaKutta4: Possible bug detected: estimates did not clear when calling dependent objects");
			}*/
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}


	delete predictor;
	predictor=nullptr;
}

void RungeKutta4::integrate(SimulaPoint::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	State l(0);
	if(defaultSpatialIntegrationLength>0) l=data.rbegin()->second.rate.length();
	timeStep(lastTime, newTime, deltaT, l.value);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	Coordinate r1(predictor->p1.rate);
	if(r1.estimate) {
		rateCalculator.calculate(predictor->t1,r1);
		//it's still possible to get an estimate here when time steps are not synchronized
		pSD->keepEstimate(r1.estimate);
		//insert updated predictor
		predictor->p1.rate=r1;
		data.rbegin()->second.rate=r1;
	}

	//half time
	Time ht(lastTime+deltaT/2);

	//step 2
	Coordinate r2;
	rateCalculator.calculate(ht,r2);
	pSD->keepEstimate(r2.estimate);
	predictor->p2.rate=r2;
	predictor->step=2;
	predictor->t2=ht;//for extrapolation of predicted rates

	//step 3
	Coordinate r3;
	rateCalculator.calculate(ht,r3);
	pSD->keepEstimate(r3.estimate);
	predictor->p2.rate=r3;

	//step 4
	Coordinate r4;
	rateCalculator.calculate(newTime,r4);
	pSD->keepEstimate(r4.estimate);
	predictor->p2.rate=r4;
	predictor->t2=newTime;

	//rungekutta 4
	MovingCoordinate newValue(predictor->p1);
	l=((r1.length()+((r2.length()+r3.length())*2)+r4.length())/6);
	newValue.rate=(r1+((r2+r3)*2)+r4)/6;
	newValue.rate.setLength(l);
	newValue.state+=newValue.rate*deltaT;
	newValue.rate=r4;

	data.insert(std::pair<Time,MovingCoordinate>(predictor->t2,newValue));

	int count(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		SimulaPoint::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->p1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//fix possible estimate loss during postIntegrationCorrection
	if(newValue.estimate) data.rbegin()->second.estimate=newValue.estimate;
	if(newValue.rate.estimate) data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate){// && pSD->numberOfDependents()){ ///@todo somehow this seems right but still causes objects to trail behind and thereby cause late updates the branches (for estimates do not create branches) and than an error.
			//allow possible dependent object to clear their estimates
			Coordinate dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy);
/*note that this not clear when the model moves ahead in time, which happens in SimulaDerivative<>::getRate
 * 			if(dummy.estimate) {
				pSD->keepEstimate(dummy.estimate);//unlock, shouldn't be necessary?
				rateCalculator.calculate(rt,dummy);
			}
			if(dummy.estimate){
				msg::warning("RungaKutta4: Possible bug detected: estimates did not clear when calling dependent objects");
			}*/
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}


	delete predictor;
	predictor=nullptr;
}

std::string RungeKutta4::getName()const{
	return "RungeKutta4";
}

//convergence solver Rates
ConvergenceSolverRates::ConvergenceSolverRates():IntegrationBase(){
}

std::string ConvergenceSolverRates::getName()const{
	return "explicitConvergence";
}

void ConvergenceSolverRates::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//update rate
	State r1(predictor->s1.rate);
	if(r1.estimate) {
		rateCalculator.calculate(predictor->t1,r1.value);
		if (rateMultiplier){
			double temp;
			rateMultiplier->get(predictor->t1, temp);
			r1.value *= temp;
		}
		//it's still possible to get an estimate here when time steps are not synchronized
		pSD->keepEstimate(r1.estimate);
		//insert updated predictor
		predictor->s1.rate=r1;
		data.rbegin()->second.rate=r1;
	}

	//iteration
	int count(0);
	int maxcount(10);
	double d(1);
	StateRate newValue(predictor->s2);
	while( d>0.00001 && count<maxcount){//
		//counter
		++count;
		//rate calculation
		newValue.rate.estimate=false;
		rateCalculator.calculate(newTime,newValue.rate.value);
		if (rateMultiplier){
			double temp;
			rateMultiplier->get(newTime, temp);
			newValue.rate.value *= temp;
		}
		//convergence criteria
		if(newValue.rate.value!=0){
			d=absolute((newValue.rate.value-predictor->s2.rate.value)/newValue.rate.value);
		}else if(predictor->s2.rate.value!=0){
			d=absolute((newValue.rate.value-predictor->s2.rate.value)/predictor->s2.rate.value);
		}else{
			d=0;
		}
		//if(!newValue.rate.estimate) d=0;
		//integration
		pSD->keepEstimate(newValue.rate.estimate);
		//newValue.state=predictor->s1.state+(predictor->s1.rate+newValue.rate)*deltaT/2;
		newValue.state=predictor->s1.state+(newValue.rate)*deltaT;
		//Preparation for next round
		predictor->s2.rate.value=newValue.rate.value;
		//from now on do backward euler for estimates
		predictor->step=2;
	}
	if(count==maxcount) msg::warning("ConvergenceSolver: needing more than "+std::to_string(maxcount)+" iterations, continuing with less precise results for "+pSD->getName());

	//store values
	data.insert(std::pair<Time,StateRate>(newTime,newValue));

	//postIntegrationCorrection
	count=(0);
	while(rateCalculator.postIntegrationCorrection(data)) {
		SimulaVariable::Table::reverse_iterator eit(data.rbegin()),peit(eit);++peit;
		eit->second.state=predictor->s1.state+peit->second.rate*deltaT;
		++count;
		if(count>10){
			msg::warning("postIntegrationCorrection: needing more than 10 iterations, continuing with imprecise results for "+pSD->getName());
			break;
		}
	}

	//fix possible estimate loss during postIntegrationCorrection
	data.rbegin()->second.estimate=newValue.estimate;
	data.rbegin()->second.rate.estimate=newValue.rate.estimate;

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			State dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy.value);
			if (rateMultiplier){
				double temp;
				rateMultiplier->get(newTime, temp);
				dummy.value *= temp;
			}
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}
	bool r(false);
	pSD->keepEstimate(r); //post integration method might add to predictor list - clear it


	//release memory
	delete predictor;
	predictor=nullptr;
}

IntegrationBase * newInstantiationConvergenceSolverRates(){
   return new ConvergenceSolverRates();
}


//-------------Solvers -----------------

BaseSolver::BaseSolver():IntegrationBase(){

}

std::string BaseSolver::getName()const{
	return "baseSolver";
}

void BaseSolver::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	msg::error("BaseSolver:: base class, please implement integrate (SimulaVariable) for the derived class");
}

void BaseSolver::integrate(SimulaPoint::Table & data, DerivativeBase & rateCalculator){
	msg::error("BaseSolver:: base class, please implement integrate (SimulaPoint) for the derived class");
}

void BaseSolver::predict(const Time&t,Coordinate& pp){
	//interpolate
	linearInterpolation(predictor->t1,predictor->t2,t,predictor->p1.state,predictor->p2.state,pp);
	//prediction
	pp.estimate=true;
}
void BaseSolver::predict(const Time&t,State &ps){
	//interpolate
	linearInterpolation(predictor->t1,predictor->t2,t,predictor->s1.state,predictor->s2.state,ps);
	//avoid sign change
	if(predictor->s0>=0 && predictor->s1>=0 && ps<0) ps.value=0;
	if(predictor->s0<=0 && predictor->s1<=0 && ps>0) ps.value=0;
	//prediction
	ps.estimate=true;
}
void BaseSolver::predictRate(const Time&t,Coordinate&pp){
	msg::error("BaseSolver::predictRate: Rates are not calculated by the solvers. Use different method for "+pSD->getPath());
}
void BaseSolver::predictRate(const Time&t,State &ps){
	msg::error("BaseSolver::predictRate: Rates are not calculated by the solvers. Use different method for "+pSD->getPath());
}

//IterativeSolver -- this not an integration method
IterativeSolver::IterativeSolver():BaseSolver(){

}

std::string IterativeSolver::getName()const{
	return "iterativeSolver";
}

void IterativeSolver::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//let rate calculator do the iteration
	pSD->ignoreLocks(true);
	rateCalculator.calculate(newTime,predictor->s2.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(newTime, temp);
		predictor->s2.value *= temp;
	}
	pSD->ignoreLocks(false);

	//clear dependencies
	pSD->keepEstimate(predictor->s2.estimate);

	//rate not calculated
	predictor->s2.rate=0;

	//store values
	data.insert(std::pair<Time,StateRate>(newTime,predictor->s2));

	//release memory
	delete predictor;
	predictor=nullptr;
}
void IterativeSolver::predict(const Time&t,State &ps){
	ps=predictor->s2;
	ps.estimate=true;
}


IntegrationBase * newInstantiationIterativeSolver(){
   return new IterativeSolver();
}

//SingleStepSolver -- this not a really an integration method
SingleStepSolver::SingleStepSolver():BaseSolver(){
}
std::string SingleStepSolver::getName()const{
	return "singleStepSolver";
}
void SingleStepSolver::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//iteration
	StateRate newValue;
	rateCalculator.calculate(newTime,newValue.value);
	if (rateMultiplier){
		double temp;
		rateMultiplier->get(newTime, temp);
		newValue.value *= temp;
	}
	pSD->keepEstimate(newValue.estimate);

	///@todo rate based on finite difference?
	newValue.rate=0;

	//store values
	data.insert(std::pair<Time,StateRate>(newTime,newValue));

	//release memory
	delete predictor;
	predictor=nullptr;
}

IntegrationBase * newInstantiationSingleStepSolver(){
   return new SingleStepSolver();
}

//convergence solver
ConvergenceSolver::ConvergenceSolver():BaseSolver(){

}

std::string ConvergenceSolver::getName()const{
	return "convergenceSolver";
}

void ConvergenceSolver::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){
	//timestepping
	Time lastTime=data.rbegin()->first, newTime, deltaT;
	timeStep(lastTime, newTime, deltaT);

	//predictor
	predictor=new SimplePredictor(data,deltaT);

	//let rate calculator do the iteration
	int count(0);
	int maxcount(10);
	double d(1);
	StateRate newValue(predictor->s2);
	while(d>0.001 && count<maxcount){
		++count;
		rateCalculator.calculate(newTime,newValue.value);
		if (rateMultiplier){
			double temp;
			rateMultiplier->get(newTime, temp);
			newValue.value *= temp;
		}
		if(newValue.value!=0){
			d=absolute((newValue.value-predictor->s2.value)/newValue.value);
		}else{
			d=absolute((newValue.value-predictor->s2.value)/0.0001);
		}
		predictor->s2.value=newValue.value;
		//release locks
		pSD->keepEstimate(newValue.estimate);
	}
	if(count==maxcount) msg::warning("ConvergenceSolver: needing more than "+std::to_string(maxcount)+" iterations, continuing with less precise results for "+pSD->getName());

	///@todo rate based on finite difference?
	newValue.rate=0;

	//store values
	data.insert(std::pair<Time,StateRate>(newTime,newValue));

	//fix dependent objects, so that can clear their estimates
	if(!newValue.estimate && pSD->numberOfDependents()){
			//allow possible dependent object to clear their estimates
			State dummy;
			//Time rt=pSD->getWallTime(newTime);
			rateCalculator.calculate(newTime,dummy.value);
			if (rateMultiplier){
				double temp;
				rateMultiplier->get(newTime, temp);
				dummy.value *= temp;
			}
			pSD->keepEstimate(dummy.estimate);//unlock, somehow locks are still set when the rate calculator is called and this returns false.
	}

	//release memory
	delete predictor;
	predictor=nullptr;
}

IntegrationBase * newInstantiationConvergenceSolver(){
   return new ConvergenceSolver();
}

//==================auto registration of the functions=================


IntegrationBase * newInstantiationDefaultIntegration(){
   return new ForwardEuler();
}
IntegrationBase * newInstantiationBackwardEuler(){
   return new BackwardEuler();
}
IntegrationBase * newInstantiationHeuns(){
   return new Heuns();
}
IntegrationBase * newInstantiationHeunsII(){
   return new HeunsII();
}
IntegrationBase * newInstantiationRungeKutta4(){
   return new RungeKutta4();
}

///@todo resolve this
//extern std::map<std::string, integrationInstantiationFunction> BaseClassesMap::getIntegrationClasses();

 //Register the module
static class AutoRegisterIntegrationFunctions {
public:
   AutoRegisterIntegrationFunctions() {
		//BaseClassesMap::getIntegrationClasses()["default"] = newInstantiationDefaultIntegration;//this produces less accurate results and is often not really faster
		//BaseClassesMap::getIntegrationClasses()["default"] = newInstantiationHeunsII;//this does not work at all
		//BaseClassesMap::getIntegrationClasses()["default"] = newInstantiationHeuns;//this works with near identical results and is a bit faster
		BaseClassesMap::getIntegrationClasses()["default"] = newInstantiationRungeKutta4;
		//BaseClassesMap::getIntegrationClasses()["default"] = newInstantiationConvergenceSolverRates;//not implemented for simulaPoint.
		BaseClassesMap::getIntegrationClasses()["Euler"] = newInstantiationDefaultIntegration;
		BaseClassesMap::getIntegrationClasses()["ForwardEuler"] = newInstantiationDefaultIntegration;
		BaseClassesMap::getIntegrationClasses()["BackwardEuler"] = newInstantiationBackwardEuler;
		BaseClassesMap::getIntegrationClasses()["Heuns"] = newInstantiationHeuns;
		BaseClassesMap::getIntegrationClasses()["HeunsII"] = newInstantiationHeunsII;
		BaseClassesMap::getIntegrationClasses()["RungeKutta4"] = newInstantiationRungeKutta4;
		BaseClassesMap::getIntegrationClasses()["explicitConvergence"] = newInstantiationConvergenceSolverRates;
		BaseClassesMap::getIntegrationClasses()["iterativeSolver"] = newInstantiationIterativeSolver;
		BaseClassesMap::getIntegrationClasses()["convergenceSolver"] = newInstantiationConvergenceSolver;
		BaseClassesMap::getIntegrationClasses()["singleStepSolver"] = newInstantiationSingleStepSolver;
   }
}p44608510843540385;
