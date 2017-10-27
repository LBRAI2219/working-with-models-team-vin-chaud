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

#include "SimulaPoint.hpp"
#include "BaseClasses.hpp"
#include "../cli/Messages.hpp"
#include "../math/VectorMath.hpp"
#include "../math/InterpolationLibrary.hpp"
#include "../tools/StringExtensions.hpp"
#if _WIN32 || _WIN64
#include <algorithm>
#endif

SimulaPoint::SimulaPoint(const std::string &newName, SimulaBase* newAttributeOf, const std::string &rateF, const std::string &integrationF, const Time &newstartTime, const Time &newEndTime ):
	SimulaTimeDriven(newName, newAttributeOf, rateF, UNITCOORDINATES, newstartTime, newEndTime)
{
	//set integration function
	setIntegrationFunction(integrationF);
	//tell others
	//todo, maybe broadcasting should always be done, simply when a simulabase get's registered?
	SimulaBase::broadcast(this);
}

SimulaPoint::SimulaPoint(SimulaBase* newAttributeOf, const Time &newstartTime, const SimulaPoint* copyThis):
	SimulaTimeDriven(newAttributeOf, newstartTime, copyThis)
{
	//set initial conditions
	table[newstartTime]=copyThis->table.begin()->second;
	//generate a new movement calculator (don't just copy pointers)
	setIntegrationFunction(copyThis->integrationF->getName());
	//tell others
	//todo, maybe broadcasting should always be done, simply when a simulabase get's registered?
	SimulaBase::broadcast(this);
}

SimulaPoint::~SimulaPoint(){
}

std::string SimulaPoint::getType()const{
	 return "SimulaPoint";

}

//TODO overload this so that typecast from a static point is possible. The user will have to do the return cast though
///@todo getRate - and get could be merged? if the return type was StateRate?
//units should be included in this call - make sure the rate units are divided by time.

void SimulaPoint::getRate(const Time &t, Coordinate &p){
	//is this a call back?
	if(callBack){
		if((t-table.rbegin()->first)<1e-5){
			//interpolate
			p=interpolate<Colum1,Colum2>(t,table).rate;
			if(p.estimate) addToPredictorList(this);
		}else{
			if(callBack==-1){//comes from getRate with only 1 entry
				msg::error("SimulaPoint::getRate: eternal loop detected.");
			}else{
				//use predictor
				addToPredictorList(this);
				integrationF->predictRate(t,p);
			}
		}
	}else{
		//forward simulation
		integrate_(t,table);

		if(table.size()==1 && table.rbegin()->second.rate.estimate){//we never did a rate calculation so do it now
			//make sure we initiated so that initial values are set (in case time=startTime)
			///@todo this is not save when postIntegration method is called
			bool pcb=callBack;
			callingOrder.push_back(this);
			callBack=-1;
			getRateCalculator()->calculate(t,p);
			callBack=pcb;
			keepEstimate(p.estimate);
			callingOrder.pop_back();
			//we store this so in case it is not an estimate integration can keep it
			table.rbegin()->second.rate=p;
		}else{
			if(t<table.rbegin()->first+TIMEERROR){
				//interpolate
				p=interpolate<Colum1,Colum2>(t,table).rate;
				if(p.estimate){
					p=interpolate<Colum1,Colum2>(t,table).rate;
					Table::reverse_iterator pit(table.rbegin());++pit;
					if(t>pit->first && !pit->second.rate.estimate && !table.rbegin()->second.estimate)p.estimate=false;
				}
			}else{
				//integrate_ refused to move forward in time beyond the newtime of the first called integration routine - simple make a prediction
				//use predictor
				p=table.rbegin()->second.rate;
				p.estimate=true;
				addToPredictorList(this);
			}

		}

	}
	if(multiplier_) {
		msg::error("SimulaPoint: not clear way to use the multiplier for "+getPath());
		double m;
		multiplier_->get(t,m);
		p*=m;
	}

	if(p.estimate){
		Time rt= getWallTime(t);
		if(rt>0 && rt-MAXTIMESTEP>t+TIMEERROR && !callBack){
			msg::warning("SimulaPoint: returning estimates even though time is still young");
		}
	}
}

void SimulaPoint::get(const Time &t, Coordinate &p){
	//is this a call back?
	if(callBack){
		if((t-table.rbegin()->first)<1e-5){
			p=interpolate<Colum1,Colum2>(t,table);
			if(p.estimate) 	addToPredictorList(this);
		}else{
			if(callBack==-1){//comes from getRate with only 1 entry
				p=table.rbegin()->second;
			}else{
				//use predictor
				addToPredictorList(this);
				integrationF->predict(t,p);
			}
		}
	}else{
		//forward simulation
		integrate_(t,table);

		if(t<table.rbegin()->first+TIMEERROR){
			p=interpolate<Colum1,Colum2>(t,table);
		}else{
			//integrate_ refused to move forward in time beyond the newtime of the first called integration routine - simple make a prediction
			//use predictor
			typename Table::reverse_iterator it(table.rbegin());
			if(table.size()>2){
				typename Table::reverse_iterator pit(it);++pit;
				linearInterpolation(pit->first,it->first,t,pit->second.state,it->second.state,p);
			}else{
				p=it->second;
			}
			p.estimate=true;
			addToPredictorList(this);
		}
	}
	if(multiplier_) {
		msg::error("SimulaPoint: not clear way to use the multiplier for "+getPath());
		double m;
		multiplier_->get(t,m);
		p*=m;
	}
	if(p.estimate){
		Time rt= getWallTime(t);
		if(rt>0 && rt-MAXTIMESTEP>t+TIMEERROR && !callBack){
			msg::warning("SimulaPoint: returning estimates even though time is still young");
		}
	}
}

void SimulaPoint::getAbsolute(const Time &t, MovingCoordinate &p){
	//current relative position
	get(t,p);
	//if referencPoint is NULL-pointer, than this point is the origin.
	if (parent_!=nullptr){
		//current absolute position & direction of the referencePoint
		Coordinate abs;
		parent_->getAbsolute(t,abs);
		//current absolute position & direction
		p+=abs;
	}
}
void SimulaPoint::getAbsolute(const Time &t, Coordinate &p){
	//current relative position
	get(t,p);
	//if referencPoint is NULL-pointer, than this point is the origin.
	if (parent_!=nullptr){
		//current absolute position & direction of the referencePoint
		Coordinate abs;
		parent_->getAbsolute(t,abs);
		//current absolute position & direction
		p+=abs;
	}
}

SimulaBase* SimulaPoint::getReference(){
	return this;
}


void SimulaPoint::getBase(const Time &t, Coordinate &p){
	//if referencPoint is NULL-pointer, than this point is the origin.
	p*=0;
	if (parent_!=nullptr){
		//current absolute position & direction of the referencePoint
		Coordinate abs;
		parent_->getAbsolute(t,abs);
		//current absolute position & direction
		p+=abs;
	}
}



//gets time based on when position was closest to p - assuming p and this have same referencepoint
void SimulaPoint::getTime(const Coordinate &p, Time &r, Time tmin, Time tmax){
	//never collect garbage
	collectGarbage_=false;
	//TODO look critical at this, it is to heavy for a simple question?
	State d(10E30);
	//search for the nearest line segment
	if(table.empty()) msg::error("SimulaPoint::getTime: Table is empty");
	if (table.size()==1) {
		r=table.begin()->first;
		return;
	}
	Table::const_iterator itMatch=table.end();
	double rl,frl(0);
	Coordinate p1,p2(table.begin()->second);
	Time t1(table.begin()->first),t2(table.begin()->first);
	for (Table::const_iterator it(table.begin()++) ; it != table.end() ; ++it){
		p1=p2;
		p2=it->second;
		t1=t2;
		t2=it->first;
		//incase point did not move here
		if(p1==p2) continue;
		//this is a distance between a point and a line.
		State nd=distance(p1,p2,  p, rl);
		if(nd<d && rl<=1 && rl>=0){
			d=nd;
			frl=rl;
			itMatch=it;
			if(d<1e-4) break;
		}
	}
	//check if r was found
	if (itMatch==table.end()) {
		if(d>1){
			msg::error("SimulaPoint::getTime: Projection of point lays wide (>10cm) outside the current simulated path");
		}else{
			msg::warning("SimulaPoint::getTime: Projection of point lays outside the current simulated path");
		}
	}

	//interpolate time
	r= t1+((t2-t1)*frl);

}

//todo, this does not forward in time, is also not needed for current use, but does require care of the programmer.
void SimulaPoint::getTangent(const Coordinate &p, Coordinate &r, Coordinate &projectionPoint){
	//never collect garbage
	collectGarbage_=false;
	//TODO look critical at this, it is to heavy for a simple question?
	State d(10E30);
	//search for the nearest line segment
	if(table.empty()) msg::error("SimulaPoint::getTime: Table is empty");
	if (table.size()==1) {
		r*=0;
		return;
	}
	Table::const_iterator itMatch=table.end();
	double rl,frl;
	Coordinate p1,p2(table.begin()->second);
	Coordinate pt1,pt2(table.begin()->second);
	for (Table::const_iterator it(table.begin()++) ; it != table.end() ; ++it){
		pt1=pt2;
		pt2=it->second;
		//incase point did not move here
		if(pt1==pt2) continue;
		//this is a distance between a point and a line.
		State nd=distance(pt1,pt2,  p, rl);
		if(nd<d && rl<=1 && rl>=0){
			p1=pt1;
			p2=pt2;
			d=nd;
			frl=rl;
			itMatch=it;
			if(d<1e-4) break;
		}
	}
	if(itMatch==table.end()){
		//never found one, just determine nearest point
		Coordinate pt1((table.begin()++)->second),pt2(table.begin()->second);
		for (Table::const_iterator it(table.begin()++) ; it != table.end() ; ++it){
			pt2=it->second;
			State nd=vectorlength(p,pt2);
			if(nd<d){
				p1=pt1;
				p2=pt2;
				d=nd;
			}
			pt1=pt2;
		}
		rl=1;
	}

	//check if r was found
	if (itMatch==table.end()) {
		if(d>10){
			msg::error("SimulaPoint::getTangent: Projection of point lays wide (>10cm) outside the current simulated path");
		}else{
			msg::warning("SimulaPoint::getTangent: Projection of point lays outside the current simulated path");
		}
	}

	//interpolate time
	r= p2-p1;
	projectionPoint=p1+r*frl;
}


Time SimulaPoint::lastTimeStep(){
	if(table.size()>1){
		Table::reverse_iterator it(table.rbegin()), pit(it);++pit;
		return it->first-pit->first;
	}else{
		return preferedT;
	}
}

SimulaPoint::Table const * SimulaPoint::getTable()const{
	return &table;
}

void SimulaPoint::setInitialValue(const Colum2 &pos){
	if(pos.estimate) msg::error("SimulaPoint::setInitialPosition: refusing to use estimate for initial position");
	//check if point already moved
	if(table.size()>1)msg::warning("SimulaPoint::setInitialPosition: Changing the initial position of '"+getName()+"' which has already moved to other positions.");
	//change position
	table[startTime]=pos;
}


void SimulaPoint::collectGarbage(const Time & t){
	if(collectGarbage_){
		//check time
		Time s=std::min(table.rbegin()->first,t);
		//make sure we keep at least 3 timesteps
		s=std::min(s,table.rbegin()->first-(5*maxT));
		if(s <= startTime) return;
		//find iterator that points to
	    table.erase(table.begin(),table.lower_bound(s));
	}
}


std::ostream &operator<<(std::ostream &os, const SimulaPoint &obj){
	return os;
}

SimulaBase* SimulaPoint::createAcopy(SimulaBase * attributeOf, const Time & startTime)const{
      return new SimulaPoint(attributeOf,startTime,this);
}

void SimulaPoint::getXMLtag(Tag &tag){
	//name
	tag.name="SimulaPoint";
	tag.closed=true;
	//attributes
	tag.attributes.clear();
	tag.attributes["name"]=getName();
	tag.attributes["prettyName"]=getPrettyName();
	tag.attributes["unit"]=getUnit().name;
	if(getStartTime()>0) tag.attributes["startTime"]=std::to_string(getStartTime());
	if(getEndTime()>0) tag.attributes["endTime"]=std::to_string(getEndTime());
	tag.attributes["function"]=getRateFunctionName();
	if(integrationF) tag.attributes["integrationFunction"]=integrationF->getName();
	if(!collectGarbage_) tag.attributes["garbageCollectionOff"]="true";
	std::string on(getObjectGeneratorName());
	if(!on.empty()) tag.attributes["objectGenerator"]=on;
	//data
	if(!table.empty()){
		std::string d;
		for(auto &dp:table){
			d+=convertToString<Colum1>(dp.first);
			d+="\t";
			d+=convertToString<Colum2>(dp.second);
			d+="\n";
		}
		tag.attributes["_DATA_"]=d;
		tag.closed=false;
	}
	//children
	if(getNumberOfChildren()) tag.closed=false;

}
