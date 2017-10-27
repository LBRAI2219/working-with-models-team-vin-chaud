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
#include "ObjectGenerator.hpp"
#include "SimulaBase.hpp"
#include "../cli/Messages.hpp"
#if _WIN32 || _WIN64
#include <algorithm>
#endif

ObjectGenerator::ObjectGenerator(SimulaBase* const npSB):
	lock(-1),
	pSB(npSB),
	//set last time for updatefunction to -1 (function has not been called)
	name("notSet"),
	run(true)
{}
Time ObjectGenerator::wallTime_(-2);
ObjectGenerator::~ObjectGenerator(){}

void ObjectGenerator::initialize(const Time &t){
	//nothing todo
}
void ObjectGenerator::generate(const Time &t){
	//error message
	msg::error("ObjectGenerator: base class called - did you forget to include a generator method in your module?");
}
void ObjectGenerator::getWallTime(Time &t){
  if(wallTime_>-1 && t>wallTime_){
	  t=wallTime_;
  }
}
void ObjectGenerator::startUpdating(){
  wallTime_=-1;
}
void ObjectGenerator::generateObjects(const Time &t){
	//check for loops
	if(wallTime_<-1.5) return; //do not run during reading of input files.
	if(lock==1) {
		//Time t0,t1;
		//SimulaTimeDriven::getGlobalTimeStepInfo(t0,t1);
		//TODO this is normally where the a predictor should be used -
		//but we can not predict objects other than no change
		//I don't think we need this check, since if the object generator than later inserts objects with wrong time, the database will complain.
		//if(lastTime_+MINTIMESTEP<t0){
		//	msg::error("ObjectGenerator::generateObjects: Update calls update and no way to signal estimates.");
		//}// //otherwise implicit step calls this again but we can not handle that.
        //todo signaling predictions here us correct, but the model tends to start lingering behind. Note, that if you want this, you should remove the reentry lock in Database::updateAll
		//pSB->addToPredictorList(pSB);
		return;
	}

	//check if we should destruct.
	if(!run) {
		//self destruct by deregistering
		pSB->setFunctionPointer();
		return;
	}

	//initialize if not already done so
	Time oldWallTime(wallTime_);
	std::size_t ps(pSB->getPredictorSize());
	if(lock==-1) {
		lock=1;
		//SimulaTimeDriven::callingOrder.push_back(pSB);
		initialize(pSB->getStartTime());
		//SimulaTimeDriven::callingOrder.pop_back();
		wallTime_=pSB->getStartTime();
		pSB->clearPredictions(ps);//this is necessary, as initialization called earlier could result in a prediction
		if(t>pSB->getStartTime()+TIMEERROR){
			//SimulaTimeDriven::callingOrder.push_back(pSB);
			generate(pSB->getStartTime());
			//SimulaTimeDriven::callingOrder.pop_back();
			pSB->clearPredictions(ps);//this is necessary, as initialization called earlier could result in a prediction
		}
	}

	lock=1;
	//run the generator //will only run if time is positive.
	if(oldWallTime>0){
		wallTime_=std::min(wallTime_,t);
	}else{
		wallTime_=t;
	}
	//SimulaTimeDriven::callingOrder.push_back(pSB);
	generate(wallTime_);
	//SimulaTimeDriven::callingOrder.pop_back();
	lock=0;
	wallTime_=oldWallTime;
	pSB->clearPredictions(ps);//this is necessary, as initialization called earlier could result in a prediction

}





