/*
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
#include "TillerFormation.hpp"
#include "../../cli/Messages.hpp"
#include "../../engine/Origin.hpp"
#include "../../engine/SimulaConstant.hpp"
#include "../PlantType.hpp"

//create tillers
TillerFormation::TillerFormation(SimulaBase* const pSB) :
	ObjectGenerator(pSB), timingOfTillers(nullptr), crownDiameter(nullptr), stress(nullptr),
	gp(nullptr), timingOfLastTiller(pSB->getStartTime()),
	tillernumber(0), newPosition(0, 0, 0) ,pot(false){
}
TillerFormation::~TillerFormation() {
}
void TillerFormation::initialize(const Time &t) {
	auto found=pSB->getName().find("potential");
	if(found!=std::string::npos)
		pot=true;

	//top of plant
	SimulaBase *top(pSB);
	PLANTTOP(top);
	//plant type
	std::string plantType;
	top->getChild("plantType")->get(plantType);
	//pointer to root parameters of the parent
	SimulaBase *shootPar(ORIGIN->getChild("rootTypeParameters")->getChild(plantType)->getChild("shoot"));
	if(pot)	{
		timingOfTillers = shootPar->existingChild("potentialTimeDelayBetweenTillers");
		if(!timingOfTillers) msg::warning("TillerFormation: potentialTimeDelayBetweenTillers not found, using timeDelayBetweenTillers");
	}
	if(!timingOfTillers) timingOfTillers = shootPar->getChild("timeDelayBetweenTillers");

	double delay;
	timingOfTillers->get(0, delay);
	timingOfLastTiller += delay;

	//newPosition, relative to seed position.
	pSB->getAbsolute(pSB->getStartTime(), newPosition);
	newPosition.y*=-1;//tillers form at the surface a
	newPosition.x=0;//tillers form at the surface a
	newPosition.z=0;//tillers form at the surface a

	//hypocotyl
	gp=dynamic_cast<SimulaPoint*>	( top->getChild("plantPosition")->getChild("hypocotyl")->getChild("growthpoint"));

	//todo positions should be read from input files
	xs.push_back(1);zs.push_back(1);
	xs.push_back(-1);zs.push_back(1.2);
	xs.push_back(0);zs.push_back(0.8);
	xs.push_back(0.1);zs.push_back(-0.9);
	xs.push_back(-0.5);zs.push_back(0.5);
	xs.push_back(0.5);zs.push_back(-0.5);
	xs.push_back(-0.5);zs.push_back(-0.5);
	xs.push_back(0.5);zs.push_back(0.5);
	xs.push_back(0.2);zs.push_back(-0.8);
	xs.push_back(0.2);zs.push_back(0.8);
	xs.push_back(-0.2);zs.push_back(-0.8);
	xs.push_back(-0.2);zs.push_back(0.8);
	xs.push_back(0.8);zs.push_back(0.2);
	xs.push_back(0.8);zs.push_back(-0.2);
	xs.push_back(-0.8);zs.push_back(0.2);
	xs.push_back(-0.8);zs.push_back(-0.2);

	//todo, this is a simple hack, as long as we do not have the true geometry of tillers, and secondary order of tillers.
	crownDiameter = shootPar->getChild("crownDiameter","cm");

	//stress factor
	if(!pot) stress=pSB->getParent(3)->existingChild("stressFactor:impactOn:tillerFormation");
}

void TillerFormation::generate(const Time &t) {
	if (t > timingOfLastTiller-TIMEERROR) {
		double rs;
		crownDiameter->get(t - pSB->getStartTime(),rs);
		rs/=2;

		if(gp){//should only happen when hypocotyl has grown sufficiently
			Coordinate directionParent,projectionPoint;
			gp->getTangent(newPosition, directionParent,projectionPoint);
			//note that both the tiller positions and the hypocotyl here are relative to the same base, namely the seed position
			newPosition.x=projectionPoint.x;
			newPosition.z=projectionPoint.z;
			gp=nullptr;
		}

		while (t > timingOfLastTiller-TIMEERROR) {
			int pos= tillernumber%xs.size();
			Coordinate shift;
			shift.x=xs[pos];
			shift.z=zs[pos];
			shift.normalize();
			shift*=rs;
			Coordinate newPositionRot(newPosition-shift);
			//name
			tillernumber += 1;
			if (tillernumber > 100)
				msg::error("TillerFormation: more than 100 tillers generated. Check your tillering rules");
			std::string nameNewTiller("tiller" + std::to_string(tillernumber));

			//create the tiller
			//create new branch by copying it from the template
			SimulaConstant<Coordinate>* pNewBranch = new SimulaConstant<Coordinate>(nameNewTiller, pSB, newPositionRot, timingOfLastTiller);
			if(!pot) pNewBranch->setFunctionPointer("tillerDevelopment");

			//check for next time
			double delay;
			timingOfTillers->get(timingOfLastTiller - pSB->getStartTime(), delay);
			if(stress){ //todo this is numerically and physiologically probably not ideal during recovery from stress.
				double s;
				stress->get(timingOfLastTiller,s);
				if(s<1e-3) s=1e-3;
				delay/=s;
			}
			timingOfLastTiller += delay;

		}}

}

//the maker function for instantiation of the class
ObjectGenerator * newInstantiationTillerFormation(SimulaBase* const pSB) {
	return new TillerFormation(pSB);
}

//have tillers form nodal roots etc.
TillerDevelopment::TillerDevelopment(SimulaBase* const pSB) :
		ObjectGenerator(pSB), timingOfNodalRoots(nullptr), hypo(nullptr), timingOfLastNodalRoot(pSB->getStartTime()), nodalRootNumber(0), newPosition(0, 0, 0) {
}
TillerDevelopment::~TillerDevelopment() {
}
void TillerDevelopment::initialize(const Time &t) {
	//top of plant
	SimulaBase *top(pSB);
	PLANTTOP(top);
	//plant type
	top->getChild("plantType")->get(plantType);
	//hypocotyl
	hypo = top->getChild("plantPosition")->getChild("hypocotyl");
	//todo, now for only one root class
	branchRootType="nodalrootsOfTillers";

	//pointer to root parameters of the parent
	SimulaBase *rootPar(
			ORIGIN->getChild("rootTypeParameters")->getChild(plantType)->getChild("hypocotyl")->getChild("branchList")->getChild(
					branchRootType));
	timingOfNodalRoots = rootPar->getChild("branchingDelay");

	//tiller template
	SimulaBase *tillerTemplate(	ORIGIN->existingChild("tillerTemplate"));
	if(tillerTemplate){
		pSB->copyAttributes(pSB->getStartTime(), tillerTemplate);
	}


}

void TillerDevelopment::generate(const Time &t) {
	//generate roots on the tillers

	auto container = hypo->getChild("branches")->getChild(branchRootType);
	//newPosition
	pSB->getAbsolute(pSB->getStartTime(),newPosition);//position of the tiller
	Coordinate refPosition;
	container->getAbsolute(pSB->getStartTime(),refPosition);//sowing position (base of hypocotyl)
	newPosition-=refPosition;

	while (t > timingOfLastNodalRoot) {
		double delay;
		timingOfNodalRoots->get(timingOfLastNodalRoot - pSB->getStartTime(), delay);
		timingOfLastNodalRoot += delay;

		//name
		nodalRootNumber += 1;
		if (nodalRootNumber > 100)
			msg::error("TillerDevelopment: more than 100 nodalroots generated. Check your tillering rules");
		std::string nameNewNodalRoot(pSB->getName() + ".nodal" + std::to_string(nodalRootNumber));

		//create the root
		SimulaConstant<Coordinate>* pNewBranch = new SimulaConstant<Coordinate>(nameNewNodalRoot, container, newPosition, timingOfLastNodalRoot);
		pNewBranch->setFunctionPointer("generateRoot");

		//insert rootType attribute of new branch
		new SimulaConstant<std::string>("rootType", pNewBranch, branchRootType, UnitRegistry::noUnit(), timingOfLastNodalRoot);
	}
	pSB->updateRecursively(t);

}

//the maker function for instantiation of the class
ObjectGenerator * newInstantiationTillerDevelopment(SimulaBase* const pSB) {
	return new TillerDevelopment(pSB);
}



//--------------------------------------------------------------------------------------------------------

NumberOfTillers::NumberOfTillers(SimulaDynamic* pSD) :
		DerivativeBase(pSD) {
	auto found=pSD->getName().find("potential");
	if(found!=std::string::npos){
		ref=pSD->getSibling("potentialTillers");
	}else{
		ref=pSD->getSibling("tillers");
	}

}
std::string NumberOfTillers::getName() const {
	return "numberOfTillers";
}

void NumberOfTillers::calculate(const Time &t, int &ci) {
	SimulaBase::List list;
	ref->getAllChildren(list,t);
	ci=(int)list.size();
}
void NumberOfTillers::calculate(const Time &t, double &ci) {
	SimulaBase::List list;
	ref->getAllChildren(list,t);//todo implement a get size
	ci=(double)list.size();
}
DerivativeBase * newInstantiationNumberOfTillers(SimulaDynamic* const pSD) {
	return new NumberOfTillers(pSD);
}





NumberOfRoots::NumberOfRoots(SimulaDynamic* pSD) :
		DerivativeBase(pSD) {
	std::string rootClass=pSD->getName().substr(8,std::string::npos);
	ref=pSD->getSibling("plantPosition")->getChild("hypocotyl")->getChild("branches")->getChild(rootClass);
}
std::string NumberOfRoots::getName() const {
	return "numberOfRoots";
}

void NumberOfRoots::calculate(const Time &t, int &ci) {
	SimulaBase::List list;
	ref->getAllChildren(list,t);
	ci=(int)list.size();
}
void NumberOfRoots::calculate(const Time &t, double &ci) {
	SimulaBase::List list;
	ref->getAllChildren(list,t);
	ci=(double)list.size();
}
DerivativeBase * newInstantiationNumberOfRoots(SimulaDynamic* const pSD) {
	return new NumberOfRoots(pSD);
}


class AutoRegisterTillerFormationInstantiationFunctions {
public:
	AutoRegisterTillerFormationInstantiationFunctions() {
		// register the maker with the factory
		BaseClassesMap::getObjectGeneratorClasses()["tillerFormation"] = newInstantiationTillerFormation;
		BaseClassesMap::getObjectGeneratorClasses()["tillerDevelopment"] = newInstantiationTillerDevelopment;
		BaseClassesMap::getDerivativeBaseClasses()["numberOfTillers"] = newInstantiationNumberOfTillers;
		BaseClassesMap::getDerivativeBaseClasses()["numberOfRoots"] = newInstantiationNumberOfRoots;
	}
};

static AutoRegisterTillerFormationInstantiationFunctions p679kik438382;




