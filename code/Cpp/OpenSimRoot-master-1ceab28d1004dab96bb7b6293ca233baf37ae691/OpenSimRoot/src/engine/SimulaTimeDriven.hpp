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

#ifndef SIMULATIMEDRIVEN_HPP_
#define SIMULATIMEDRIVEN_HPP_

#include "SimulaDynamic2.hpp"
#include <set>
class IntegrationBase;
///@todo, now hard coded, finally this should be part of a function
// ///@todo this schould not be larger than max timestep, ever
#define PREFEREDTIMESTEP 0.20000000 //by default zero - it will be ignored - if set - it will be used
///@todo: the smaller this the less warnings about min timestep but the more warnings about error from forward euler and the faster execution time - something is wrong
#define MINTIMESTEP 0.1000000000
#define MAXTIMESTEP 0.20000000000
#define LASTTIMESTEP 0.00000000


class SimulaTimeDriven:public SimulaDynamic{
protected:
	typedef std::map<Time, StateRate > SimulaVariableTable;
	typedef std::map<Time, MovingCoordinate > SimulaPointTable;
	typedef std::set<SimulaBase* > UniquePointerList;
	Time minT;
	Time maxT;
	Time preferedT;
	//bool syncTimeStep_;
	std::size_t predictorPosition;
	void integrate_(const Time & t, SimulaPointTable & data);
	void integrate_(const Time & t, SimulaVariableTable & data);
	IntegrationBase * integrationF;
	int callBack;
	SimulaBase* dependLock_;
	virtual void releaseDependents();
	bool forwardOnly; //allows this object to clear estimates and go only forward in time
	static bool ignoreLocks_;
private:
	UniquePointerList dependentList;
public:
	static  PointerList callingOrder;
	virtual Time &minTimeStep();
	virtual Time &maxTimeStep();
	virtual Time &preferedTimeStep();
	virtual Time lastTimeStep();
	virtual void getTimeStepInfo(Time &st,Time &rt)const;
	static Time getWallTime(const Time &t);
	static void getGlobalTimeStepInfo(Time &st,Time &rt);
	//virtual bool &syncTimeStep();///todo redundant now we do timestep snapping.
	virtual void keepEstimate(bool &);
	void setIntegrationFunction(const std::string &nameIF);
	std::string getIntegrationFunction()const;
	virtual void avoidPredictorCorrectedLoops();

	//timestep function
	SimulaTimeDriven(SimulaBase* const newAttributeOf, const std::string &functionName, const Unit &newUnit, const Time &newstartTime);
	SimulaTimeDriven(const std::string &newName, SimulaBase* newAttributeOf, const std::string &functionName, const Unit &newUnit, const Time &newStartTime);
	SimulaTimeDriven(const std::string &newName, SimulaBase* newAttributeOf, const std::string &functionName, const Unit &newUnit, const Time &newStartTime, const Time &newEndTime);
	SimulaTimeDriven(SimulaBase* newAttributeOf, const Time &newstartTime, const SimulaTimeDriven* copyThis);
	virtual ~SimulaTimeDriven();
	virtual bool setDependent(SimulaBase* me);
	virtual void unLockDependency();
	void ignoreLocks(const bool &);
	std::size_t numberOfDependents()const;

	virtual SimulaBase* createAcopy(SimulaBase * attributeOf, const Time & startTime)const;
	virtual void getXMLtag(Tag &tag);

};



#endif
