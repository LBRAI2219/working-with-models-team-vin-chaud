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
#include "NutrientUptake.hpp"
#include "../../cli/Messages.hpp"
#include "../../engine/Origin.hpp"
#include "../PlantType.hpp"
#include <math.h>		// std::pow
#if _WIN32 || _WIN64
#include <algorithm>
#endif

/* Cash-Karp method
	this means two methods of orders 5 and 4

 Butcher Tableau:
 	0    |
	1/5  |	1/5
	3/10 |	3/40 	    9/40
	3/5  |  3/10 	   -9/10 	      6/5
    1    |-11/54 	    5/2 	    -70/27 	    35/27
    7/8  |1631/55296 	175/512 	575/13824 	44275/110592 	253/4096
 -------------------------------------
		 | 37/378 	    0 	250/621 	    125/594 	          0 	512/1771
		 | 2825/27648   0 	18575/48384 	13525/55296 	277/14336 	1/4
 */



Barber_cushman_1981_nutrient_uptake_explicit::Barber_cushman_1981_nutrient_uptake_explicit(SimulaDynamic* pSD) :TotalBaseLabeled(pSD){
}

void Barber_cushman_1981_nutrient_uptake_explicit::calculate(const Time &t, double &fluxDensity) {
}

//***************************************************************************/


std::string Barber_cushman_1981_nutrient_uptake_explicit::getName()const{
	return "barber_cushman_1981_nutrient_uptake_explicit";
}

DerivativeBase * newInstantiationBarber_cushman_1981_nutrient_uptake_explicit(SimulaDynamic* const pSD) {
	return new Barber_cushman_1981_nutrient_uptake_explicit(pSD);
}



//Register the module
static class AutoRegisterNutrientClassesInstantiationFunctionszj56 {
public:
	AutoRegisterNutrientClassesInstantiationFunctionszj56() {
		BaseClassesMap::getDerivativeBaseClasses()["barber_cushman_1981_nutrient_uptake_explicit"]	= newInstantiationBarber_cushman_1981_nutrient_uptake_explicit;
	}

}p7bg876h868h876h;


double BarberCushmanSolver::theta = -1;

const double BarberCushmanSolver::A11 = 0.2;    // 1/5;
const double BarberCushmanSolver::A21 = 0.075;  // 3/40;
const double BarberCushmanSolver::A22 = 0.225;  // 9/40;
const double BarberCushmanSolver::A31 = 0.3;    // 3/10;
const double BarberCushmanSolver::A32 = -0.9;   // -9/10;
const double BarberCushmanSolver::A33 = 1.2;    // 6/5 ;
const double BarberCushmanSolver::A41 = -11./54.;
const double BarberCushmanSolver::A42 = 2.5;    // 5/2;
const double BarberCushmanSolver::A43 = -70./27.;
const double BarberCushmanSolver::A44 = 35./27.;
const double BarberCushmanSolver::A51 = 1631./55296.;
const double BarberCushmanSolver::A52 = 175./512.;
const double BarberCushmanSolver::A53 = 575./13824.;
const double BarberCushmanSolver::A54 = 44275./110592.;
const double BarberCushmanSolver::A55 = 253./4096.;

/* order 5 */
const double BarberCushmanSolver::B11 = 37./378.;

const double BarberCushmanSolver::B13 = 250./621.;
const double BarberCushmanSolver::B14 = 125./594. ;

const double BarberCushmanSolver::B16 = 512./1771.;

/* order 4 */
const double BarberCushmanSolver::B21 = 2825./27648.;

const double BarberCushmanSolver::B23 = 18575./48384.;
const double BarberCushmanSolver::B24 = 13525./55296.;
const double BarberCushmanSolver::B25 = 277./14336.;
const double BarberCushmanSolver::B26 = 0.25;

const double BarberCushmanSolver::C1 = 0.2;
const double BarberCushmanSolver::C2 = 0.3;
const double BarberCushmanSolver::C3 = 0.6;
const double BarberCushmanSolver::C4 = 1;
const double BarberCushmanSolver::C5 = 0.875;

bool BarberCushmanSolver::MIXED = true;
bool BarberCushmanSolver::automatic_startstep = true;



BarberCushmanSolver::BarberCushmanSolver():
IntegrationBase(),dr(0),b(0),Ci(0),r0(0),r1(0),circ(0),rootHairDensityConversionFactor(0),
Cmin(0),res(0),increaseTimeStep(0),rt(-1),st(-1),
heuns(false),
sKm(nullptr),
sImax(nullptr),
sDe(nullptr),
sv0(nullptr),
competition(nullptr),
slh(nullptr),
sNh(nullptr),
srh0(nullptr),
rootDiameter(nullptr),
spatialRootDensity(nullptr),
rootLength(nullptr),
depletionZone(nullptr),
depth(0,0,0),
lt(-1),
k(0),
msgcount(0),
firststep(true),
h(0),
h_new(0),
os()
{}

void BarberCushmanSolver::doDelayedConstruction(){
	//initialization (t=0)
	Time t=pSD->getStartTime();
	//plant type
	std::string plantType;
	PLANTTYPE(plantType,pSD);
	//root type
	std::string rootType;
	pSD->getParent(4)->getChild("rootType")->get(rootType);
	//get nutrient
	std::string nutrient(pSD->getParent()->getName());
	//position
	pSD->getAbsolute(pSD->getStartTime(), depth);

	//plant parameters
	SimulaBase *param(GETROOTPARAMETERS(plantType,rootType));
	double dummy;
	//TODO assuming diameter to be constant
	rootDiameter=pSD->getParent(2)->getChild("rootDiameter", "cm");
	rootDiameter->get(pSD->getStartTime(),r0);r0/=2.;
	if(r0<=0)msg::error("BarberCushmanExplicit: root diameter is <= 0");

	if(pSD->getUnit()=="umol"){
		rootLength=pSD->getParent(2)->getChild("rootSegmentLength", "cm");
	}else{
		rootLength=nullptr;
		if(pSD->getUnit()!="umol/cm") msg::error("BarberCushmanExplicit: "+pSD->getPath()+" should have unit umol or umol/cm");
	}

	param=param->getChild(nutrient);
	sImax=param->getChild("Imax", "umol/cm2/day");
	sKm=param->getChild("Km", "umol/ml");
	param->getChild("Cmin", "umol/ml")->get(Cmin);
	competition=param->existingChild("competition");
	if(competition){
		bool cmp;
		competition->get(cmp);
		if(!cmp) competition=nullptr;
	}else{
		competition=pSD; //default is on
	}
	sv0=pSD->getParent(2)->existingChild("rootWaterUptake","cm3/cm");
	//soil parameters
	SimulaBase* paramsoil=ORIGIN->getChild("environment")->getChild("soil")->getChild(nutrient);
	sDe=paramsoil->getChild("diffusionCoefficient", "cm2/day");
	paramsoil->getChild("bufferPower")->get(t,depth,b);
	paramsoil->getChild("r1-r0", "cm")->get(t,depth,r1);r1+=r0;

	//volumetric water content
	//note barber-cushman can not simulate variable water content.
	if(theta<0){
		theta=1; //default is one, which is not a good default, but does make the code backward compatible.
		SimulaBase* p(paramsoil->existingSibling("water"));
		if(p){
			p=p->existingChild("volumetricWaterContentInBarberCushman");
			if(p){
				p->get(t,depth,theta);
			}else{
				msg::warning("BarberCushmanExplicit: ../soil/water/volumetricWaterContentInBarberCushman not given, using 1 for backward compatibility reasons.");
			}
		}
	}

	paramsoil->getChild("concentration","umol/ml")->get(t,depth,Ci);
	Ci*=theta;//correction of Ci for volumentric water content.
	if (Ci < 0)msg::error("BarberCushmanExplicit: initial concentration < 0");

	//root competition
	spatialRootDensity=pSD->getParent(2)->existingChild("spatialRootDensity","cm/cm3");
	if(spatialRootDensity) {
		spatialRootDensity->get(pSD->getStartTime(),dummy);
		if(dummy>0){
			//change r1 so it is the mid-distance between roots
			dummy=1/sqrt(dummy);
			if(dummy<r1){
				r1=std::min(r1,dummy);
				//change the Ci concentration
				//The idea is that this is a new root growing in a already some what depleted soil volume
				//Add up the uptake of the other root segments and recalculate Ci
				Coordinate center;
				pSD->getAbsolute(pSD->getStartTime(),center);

				SimulaBase::Positions list,list2;
				double radius(1);
				pSD->getPositionsWithinRadius(pSD->getStartTime(),center,radius,list);
				//eliminate nodes from same root
				SimulaBase *root1(pSD->getParent(3)), *root2;
				for(auto & it:list){
					root2 = it.second->getParent(2);
					//if(root1==root2) list.erase(it--);//note this must be post decrement of iterator, it will decrement the iterator and than return old copy to erase.
					if(root1!=root2) list2.insert(it);
				}
				//total uptake in umol
				double uptake(0), u;
				SimulaBase *p;
				for(auto& it:list2){
					p=it.second->getChild(nutrient)->getChild("rootSegmentNutrientUptake");
					if(p!=pSD){//do not include selve.
						p->get(pSD->getStartTime(),u);
						uptake+=u;
					};
				}
				//adjust Ci by subtracting uptake.
				double volume=radius*radius*radius*(4/3)*PI; // volume of the sphere in which uptake took place in cm3
				double CiOri(Ci);
				Ci-= ((uptake/volume)/b);
				if(Ci<0.1*CiOri) msg::warning("BarberCushmanExplicit: adjusting Ci to less than 10% of original Ci");
				if(Ci<0) Ci=0;
			}
		}
	}

	//exudates factor
	SimulaBase* p=param->existingChild("exudatesFactor");
	if(p){//found factor for exudates
		p->get(t,depth,dummy);
		Ci*=dummy;
		b/=dummy;
	}

	//estimate dr automatically
	k=31;//from experience 31 is a bit small, 61 is better, but takes much memory and time
	dr=(r1-r0)/(double)k;
	/*if(nutrient!="nitrate"){
		if(dr>0.9*r0) dr=0.9*r0;
		if(dr<0.3*r0){
			dr=0.3*r0;
		}
	}else{
		if(dr>0.9*r0) msg::warning("BarberCushman: running nitrate, but this model does not do well on nitrate. Using too large dr to avoid a crash");
	}*/

	//determine k
	//k=(int)ceil((r1-r0)/dr);
	//if(k>100) msg::warning("Barber-Cushman: using large arrays: k>100");
	r1=r0+(double)k*dr;


	//circumference
	circ=r0*2*PI;

	//set local time
	lt=pSD->getStartTime();

	//set initial concentrations
	X_old.resize(k);
	for (std::size_t i(0); i < k; ++i) X_old[i] = Ci;

	//roothair parameter
	slh = pSD->getParent(2)->existingChild("rootHairLength","cm");
	if(slh){
		srh0 = pSD->getParent(2)->getChild("rootHairDiameter","cm");
		sNh = pSD->getParent(2)->getChild("rootHairDensity");
	}else{
		sNh=nullptr;
	}
	if(sNh){
		if(sNh->getUnit()=="#/cm2"){
			rootHairDensityConversionFactor=r0*2*PI;
		}else if(sNh->getUnit()=="#/cm"){
			rootHairDensityConversionFactor=1;
		}else{
			msg::error("BarberCushmanExplicit: unknown unit for roothair density. Should be #/cm or #/cm2.");
		}
	}

	//initial flux density at time 0
	double Imax,Km;
	sImax->get(lt,Imax);
	sKm->get(lt,Km);
	Km*=theta; //correction for volumetric water content
	Cmin*=theta; //correction for volumetric water content
	b/=theta; //it's probably correct to adjust b as well, but it would confuse the definition of b. A fixed b is not very good anyway, and in the future may need to change.
	res=circ*(Imax*Cmin/(Km+Cmin));

	//depletion zone
	SimulaBase* dz(pSD->existingSibling("radiusDepletionZone"));
	if(dz){
		depletionZone= dynamic_cast < SimulaTable<Time> * > (dz);
	}else{
		depletionZone=nullptr;
	}
	if(depletionZone) {
		depletionZone->set(pSD->getStartTime(),r0);
	}else{
		msg::warning("Depletionzone not found");
	}

	//try  to choose 'smartly' the right timestep.
	p=paramsoil->existingChild("increaseTimeStep");
	if(p) {
		p->get(t,depth,increaseTimeStep);
		//guess initial timestep
		pSD->preferedTimeStep()=std::max(std::min(100*dr*pSD->preferedTimeStep(),pSD->maxTimeStep()),pSD->minTimeStep());
	}else{
		increaseTimeStep=1.05;
		msg::warning("BarberCushmanExplicit:: increaseTimeStep not set, defaulting to 1.05 and not adjusting timestep for root radius",3);
	}

	msgcount=1;

	//method
	SimulaTimeDriven * std=dynamic_cast < SimulaTimeDriven * > (pSD);
	std::string m=std->getIntegrationFunction();
	if(m!="HeunsII"){
		heuns=false;
	}

	//debugging files (additional output)
	bool writeFiles=false;
	p=param->existingChild("writeFiles");
	if(p) p->get(writeFiles);

	if(writeFiles){
		std::string filename("barber_cushman_explicit.tab");
		os.open( filename.c_str() );
		if ( !os.is_open() ) msg::warning("generateTable: Failed to open "+filename);
		if ( os.is_open() ){
			os<<std::endl<<0.0;
			for(auto Cn_i: X_old){
				os<<"\t"<<Cn_i;
			}
		}
	}

} // end constructor


std::string BarberCushmanSolver::getName()const{
	return "BarberCushmanSolver";
}


void BarberCushmanSolver::predict(const Time&,State &){
	msg::error("BarberCushmanExplicit: Predicting requested but no code for a prediction implemented");
}
void BarberCushmanSolver::predictRate(const Time&,State &){
	msg::error("BarberCushmanExplicit: Predicting requested but no code for a prediction implemented");
}
void BarberCushmanSolver::getTimeStepInfo(Time & st1, Time & rt1)const{
	st1=st;
	rt1=rt;
}

BarberCushmanSolver::~BarberCushmanSolver(){
	if ( os.is_open() ) os.close();
}


void BarberCushmanSolver::integrate(SimulaVariable::Table & data, DerivativeBase & rateCalculator){

	if (firststep){
		doDelayedConstruction();
		//firststep is set false in automatic timestep, but what if that is not done....?
	}

	//time
	Time t=data.rbegin()->first;

	//adjust k if competition increased
	if(spatialRootDensity) {
		double dummy;
		spatialRootDensity->get(pSD->getStartTime(),dummy);
		if(dummy>1){
			dummy=1/sqrt(dummy);
			r1=std::min(r1,dummy);
			std::size_t ok=k;
			k=(int)ceil((r1-r0)/dr);
			if(k<ok){
				r1=r0+(double)k*dr;
			}else if(k>ok){
				k=ok;
			}else{ //do nothing
				;
			}
		}
	}

	//current root radius
	double rdz;
	Time bt(t-0.2);if(bt<pSD->getStartTime()) bt=pSD->getStartTime();//todo see if this hack is still needed.
	rootDiameter->get(bt,rdz);rdz/=2;

	//check the soil concentration
	double fluxDensity;
	bool exit(false);
	if(Ci<=Cmin) {
		fluxDensity=0;
		exit=true;
	}else if(rdz>2*r0){
		//secondary growth destroyed cortex, assume uptake is zero
		fluxDensity=0;
		exit=true;
	}else if(k<4){//severe competition
		msg::warning("Barber-Cushman: severe competition: k<4, you may want to use smaller delta r for "+pSD->getParent()->getName());
		fluxDensity=0;
		exit=true;
	}

	if(exit){
		firststep = false;
		StateRate & oldValue=data.rbegin()->second;
		StateRate newValue;
		Time newTime(t+pSD->maxTimeStep());
		newValue.state.value=oldValue.state.value+h*(oldValue.rate.value+fluxDensity)/2; // trapezoidal sum rule
		pSD->keepEstimate(newValue.state.estimate); //true if the values depend on predictors from objects earlier in the calling sequence
		newValue.rate.value=fluxDensity;
		pSD->keepEstimate(newValue.rate.estimate); //true if the values depend on predictors from objects earlier in the calling sequence
		data.insert(std::pair<Time,StateRate>(newTime,newValue));
		if (depletionZone) depletionZone->set(newTime,rdz);
		return; //dont' run the model at all, ever
	}

	//root radius and the distance of all finite element nodes from the center of the root node
	//can go into constructor
	//the only part does goes here is that k is changed, so the array might have to be shortened.
	std::valarray<double> r(0.,k);
	r[0] = r0;
	for (unsigned int i(1); i < k; ++i) {
		r[i] = r[i - 1] + dr;
	}

	//Get the different parameters
	///@todo this may be numerically not completely valid - barber assumes Imax, Km and E V0 De constant over time.
	//Now it is not, but we assume them constant during the timestep
	double Imax,Km,v0(0.00864),De;
	sImax->get(t-pSD->getStartTime(),Imax);
	sKm->get(t-pSD->getStartTime(),Km);Km*=theta;//todo!
	if(sv0) {
		sv0->getRate(t,v0);// cm3 water/cm root/day
		//todo check this: divide v0 with the assumed constant volumetric water content so amount of water taken up is expressed in cm3 soil volume. This because the concentration is expressed in cm3 soil volume
		//We need to do this to make cm3 equal not refer to amount of water but an amount of water per cubic cm3 of soil.
		v0/=theta; //rocksprings 2009 about 30% volumetric water content.
		v0/=r[0];//cm3/cm/day / cm -> cm/day
	}
	sDe->get(t,depth,De);

	//*************************************
	//constants

	//roothairs
	//TODO Imax and Km etc could be different for roothairs
	double Imaxh(Imax);

	//double Kmh(Km), Cminh(Cmin);
	std::valarray<double> vol(0.,k),Ah(0.,k),S4(0.,k),Ch(0.,k),Ih(0.,k);
	if(slh){ //if there are root hairs
		std::valarray<double> rh1(0.,k);
		double Nh;
		sNh->get(t,Nh);
		if(Nh){
			double rh0,lh;
			slh->get(t,lh);
			srh0->get(t,rh0);rh0/=2;
			sNh->get(t,Nh);Nh*=rootHairDensityConversionFactor;
			for (unsigned int i(0); i < k; ++i){
				rh1[i] = sqrt(r[i]*PI/(2*Nh)); // half distance
				//vol[i] = (PI*square(r[i]+0.5*dr))-(PI*square(r[i]-0.5*dr));
				vol[i]=(PI*square(r[i]+dr))-(PI*square(r[i]));
				Ah[i]  = (Nh*std::min(dr,std::max(0.,(lh+r0)-r[i]))*2.*PI*rh0)/vol[i];
				S4[i]  = ((Imaxh*rh0)/(De*b))*log(rh1[i]/(1.6487212707*rh0)); // exp(0.5) = 1.6487212707
			}
		}
	}

	// RK solving procedure ...

	std::valarray<double> Pe(0.,k);
	double vmax(0.);
	for(std::size_t i=0; i < k; ++i) {
		double v = fabs( 1/r[i]*(De+v0*r0/b) ); // theoretically v0 can be negative
		Pe[i] = v*dr/De ; // <= 2
		if(Pe[i] > 2) msg::warning("BarberCushmanExplicit: Peclet criteria not fulfilled. Upwind used in compartment.",2);
		vmax = vmax > v ? vmax : v;
	}

	std::valarray<double> X(Ci,k);

	std::valarray<double> w1(0.,k);
	std::valarray<double> w2(0.,k);
	std::valarray<double> w3(0.,k);
	std::valarray<double> w4(0.,k);
	const double dr2 = dr*dr;
	for(size_t i=0; i < k; ++i){
		double advection_term =  1.0/r[i]*(De+v0*r0/b) /dr;
		if(advection_term > 0){ //  upwind in root direction
			if(MIXED && Pe[i]<=1.8){
				// central
				w1[i] = advection_term*(-0.5)       +De/dr2;
				w2[i] =                           -2.*De/dr2;
				w3[i] = advection_term*(0.5)        +De/dr2;
				//w4[i] = 0;
			}
			else if(MIXED==false || (MIXED && Pe[i]>2)){
				// upwind
				w1[i] =                                De/dr2;
				w2[i] = advection_term*(-1.5)       -2.*De/dr2;
				w3[i] = advection_term*2.0            +De/dr2;
				w4[i] = advection_term*(-0.5);
			}
			else if(MIXED && Pe[i]>1.8){ // linear combination
				advection_term = 0.5*advection_term;
				w1[i] = advection_term*(-0.5)     + De/dr2;
				w2[i] = advection_term*(-1.5)    -2.*De/dr2;
				w3[i] = advection_term*(2.5)      + De/dr2;
				w4[i] = advection_term*(-0.5);
			}
		}
		else{ // downwind
			// use central differences as workaround !
			w1[i] = advection_term*(-0.5)       +De/dr2;
			w2[i] =                           -2.*De/dr2;
			w3[i] = advection_term*(0.5)        +De/dr2;
			//w4[i] = 0.;
			msg::warning("BarberCushmanExplicit: downwind! but central differences are used for advection.",2);
		}
	}

	// initialization of RK coeff's
	std::valarray<double> k1(0.,k);
	std::valarray<double> k2(0.,k);
	std::valarray<double> k3(0.,k);
	std::valarray<double> k4(0.,k);
	std::valarray<double> k5(0.,k);
	std::valarray<double> k6(0.,k);

	std::valarray<double> Xo4(0.,k);
	std::valarray<double> Xo5(0.,k);

	const double TOL =  1e-4; // 1e-5; // relative error tolerance of RK45

	const double h_max =  pSD->maxTimeStep();
	const double h_min =  pSD->minTimeStep();


	// first time step size guess (ref.: Hairer I)
	double h=pSD->preferedTimeStep();//timestep (dt)
	if (firststep && automatic_startstep) {

		// FUNCTION EVALUATION AT THE INITIAL POINT
		std::valarray<double> IhRK = calcIhRK(X, S4, Imaxh, Ah, Km, Cmin, k);

		std::valarray<double> f1 = rk_step(X, w1, w2, w3, w4, De, b, Cmin, r0, v0, Km, Imax, dr, IhRK, k, r);

		double d0 = L2Norm(X);
		double d1 = L2Norm(f1);
		double h0;
		if (d0 < 1e-5 || d1 < 1e-5) {
			h0 = 1e-6;
		} else {
			h0 =  d0 / d1; // norm(X)/norm(f1);
			if(h0>pSD->maxTimeStep()) h0=pSD->maxTimeStep();
			h0*=0.01;
		}

		// PERFORM EXPICIT EULER STEP y1= y0+ h0*f1(x0,y0)
		X = X + h0 * f1;

		// COMPUTE f(x0+h0, y1) , where y1 := X
		IhRK = calcIhRK(X, S4, Imaxh, Ah, Km, Cmin, k);

		std::valarray<double> f2 = rk_step(X, w1, w2, w3, w4, De, b, Cmin, r0, v0, Km, Imax, dr, IhRK, k, r);

		double d2 = L2Norm(f2 - f1) / h0; // as an estimate of the second derivative of the solution

		double h1;
		if (std::max(d1, d2) < 1e-15) {
			h1 = std::max(1e-6, h0 * 1e-3);
		} else {
			h1 = std::pow(0.01 / std::max(d1, d2), 0.2); // q+1=5 is order of local error
		}
		//todo don't want h to be larger than hmax
		h = std::min(100 * h0, h1);
		//h = 0.01;
		//h = std::max(h_min, h);
		pSD->preferedTimeStep()=h;
	}
	firststep = false;
	// else step size h is user input (or "small" default value.)


	//make sure the timestep is synchronized.
	//pSD->preferedTimeStep()=0.01;
	double newTime = t + h;
	st=t;
	rt=newTime;
	double h_ori=h;
	timeStep(t, newTime, h); // synchroniziation

	bool rejection=true;
	double facmax = 2.0; // usually between 1.5 and 5
	double facmin = 0.2; // not more than [1/facmin]-times reduction

	double safety = 0.9; // or (0.25)^(0.25); // given safety factor

	std::valarray<double> IhRK(0.,k);

	// ********* loop for re-calculating am h-step *************************
	while(rejection){
		//X_0 = X_old;

		// step 1
		// X = X_old;
		IhRK = calcIhRK( X_old, S4, Imaxh, Ah, Km, Cmin, k);

		// F( t , X_0 )
		k1 = rk_step(X_old, w1,w2,w3,w4,De,b,Cmin,r0,v0,Km,Imax ,dr, IhRK, k, r);

		// step 2
		X = X_old + h*A11 *k1;

		IhRK = calcIhRK( X, S4, Imaxh, Ah, Km, Cmin, k);

		// F( t+h*C1 , X_0+h*A11 *k1 )
		k2 = rk_step(X, w1,w2,w3,w4,De,b,Cmin,r0,v0,Km,Imax ,dr, IhRK, k, r);

		// step 3
		X = X_old + h*(A21*k1 + A22*k2);

		IhRK = calcIhRK( X, S4, Imaxh, Ah, Km, Cmin, k);

		// F( t+ C2*h , X_0 + h*(A21*k1 + A22*k2) )
		k3 = rk_step(X, w1,w2,w3,w4,De,b,Cmin,r0,v0,Km,Imax ,dr, IhRK, k, r);

		//  step 4
		X = X_old + h*(A31*k1 +A32*k2 + A33*k3);

		IhRK = calcIhRK(X, S4, Imaxh, Ah, Km, Cmin, k);

		// F( t+C3*h , X_0 + h*(A31*k1 +A32*k2 + A33*k3) )
		k4 = rk_step(X, w1,w2,w3,w4,De,b,Cmin,r0,v0,Km,Imax ,dr, IhRK, k, r);

		//  step 5
		X = X_old + h*(A41*k1 +A42*k2 +A43*k3 + A44*k4);

		IhRK = calcIhRK( X, S4, Imaxh, Ah, Km, Cmin, k);

		// F( t+C4*h , X_0+ h*(A41*k1 +A42*k2 +A43*k3 + A44*k4) )
		k5 = rk_step(X, w1,w2,w3,w4,De,b,Cmin,r0,v0,Km,Imax ,dr, IhRK, k, r);
		// k5 is needed for the lower order solution and is computed with a
		// full h-step, see C4:=1.

		//  step 6
		X = X_old + h*(A51*k1 + A52*k2 + A53*k3 + A54*k4 + A55*k5);

		IhRK = calcIhRK( X, S4, Imaxh, Ah, Km, Cmin, k);

		// F( t+C5*h , X_0+ h*(A51*k1 + A52*k2 + A53*k3 + A54*k4 + A55*k5) )
		k6 = rk_step(X, w1,w2,w3,w4,De,b,Cmin,r0,v0,Km,Imax ,dr, IhRK, k, r);

		Xo5 = X_old + h*(B11*k1+ B13*k3 + B14*k4          + B16*k6); // 5th
		Xo4 = X_old + h*(B21*k1+ B23*k3 + B24*k4 + B25*k5 + B26*k6); // 4th order solution

		double local_error = L2Norm(Xo5 - Xo4)/L2Norm(Xo5);


		// step size control algorithm:
		// ----------------------------
		if (local_error <= TOL){    // continue with larger step
			if(h+TIMEERROR>h_ori){
				//this should only happen when we did not use a small h because of synchronization of the timestep.
				double h0 = 1 + safety* (std::pow((TOL/local_error),0.2)-1);
				h_new = h*std::min(h0, facmax); // prevents the code from too large steps

				if( h_new > h_max){
					h_new = h_max;
				}
			}//else leave the timestep for what it is, for local_error is based on a small timestep to snap to the synchronization time.
			X_old = Xo5;//update x_old to new solution
			rejection = false;
			//std::cout<<std::endl<<this<<"h "<<std::scientific<<h<<" hnew "<<h_new<<" local error "<<local_error<<" Csurf "<<X_old[1]<<"K"<<k;
		}else if(h<h_min+TIMEERROR){
			//do not change h
			msg::warning("BarberCushmanExplicit: using miniminum timestep and not meeting tolerance criteria. Results might be instable");
			rejection = false;
		}else{  // if(local_error >= err_max)    //  repeat with smaller step

			double h0 =  safety* std::pow((TOL/local_error),0.25);//how (100%) much to decrease the step size

			h_new =  h*std::max(h0, facmin); // prevents the code from too small steps
			if( h_new < h_min){
				h_new = h_min;
			}
			h=h_new;
			newTime=t + h; // set new time
			rejection = true;
			facmax=1;
			msg::warning("BarberCushmanExplicit: retrying with a smaller timestep",3);
		}
	} // **** end while loop ***********************

	if(fabs( vmax*(h_new/dr) ) > 1. ){
		h_new = dr/vmax; //meet the CFL criteria
		if(h_new<pSD->minTimeStep()) msg::warning("BarberCushmanExplicit: min timestep does not allow CFL criteria to be met for at least one compartment. A smaller min dt is needed or a larger dr."); // <=1, CFL condition, same as courant number
	}
	pSD->preferedTimeStep() = h_new; // set dt for global stepping precedure


	//calculate flux into the root using second order (r0 is currently circumference - see above)
	double Ca(X_old[0]-Cmin);
	fluxDensity=circ*(Imax*Ca/(Km+Ca));
	if(Ah[0]){
		//include uptake by hairs
		for(unsigned int i(0); i<k && Ah[i]>0 ; ++i){
			double d1=calcCh(X_old[i],Km,Cmin,S4[i]);
			double d3=Ah[i]*vol[i];
			d1=calcIh(d1,Imaxh,Km,Cmin,d3);
			fluxDensity+=(d1);
		}
	}

	if ( os.is_open() ){
		os<<std::endl<<t;
		for(auto Cn_i: X_old){
			os<<"\t"<<Cn_i;
		}
	}

	//convert to per cm to per root segment length
	if(rootLength){
		double l;
		rootLength->get(t,l);
		fluxDensity*=l;
	}

	//do some checks
	if(std::isnan(fluxDensity))msg::error("BarberCushmanExplicit: Results is NaN");
	if(fluxDensity<0) {
		msg::warning("BarberCushmanExplicit: getting negative values for uptake, returning 0");
		fluxDensity=0;
	}

	//check for increasing concentrations. This happens when dr is too large
	if(X_old[0]>1.1*Ci)
		msg::warning("BarberCushmanExplicit: concentrations at root surface are increasing. massflow>uptake? or dr too large?");

	//find edge of depletion zone
	if(depletionZone){
		std::size_t edgeDepletionZone(0);
		for (int i(k-1); i >= -1; --i) {
			if (X_old[i]< 0.95*Ci){
				edgeDepletionZone=i;
				break;
			}
		}
		if(edgeDepletionZone>=k-1){ //edge is not found and competition is effective
			msg::warning("BarberCushman: competition effective.");
		}
		rdz+=r[edgeDepletionZone]-r[0];
		depletionZone->set(lt,rdz);
	}

	StateRate & oldValue=data.rbegin()->second;
	StateRate newValue;
	newValue.state.value=oldValue.state.value+h*(oldValue.rate.value+fluxDensity)/2; // trapezoidal sum rule
	pSD->keepEstimate(newValue.state.estimate); //true if the values depend on predictors from objects earlier in the calling sequence
	newValue.rate.value=fluxDensity;
	pSD->keepEstimate(newValue.rate.estimate); //true if the values depend on predictors from objects earlier in the calling sequence

	data.insert(std::pair<Time,StateRate>(newTime,newValue));
	st=-1;
	rt=-1;


}


void BarberCushmanSolver::timeStep(const Time & lastTime, Time & newTime,
		Time & deltaT, const double & l) {

	deltaT = pSD->preferedTimeStep();

	//adjust to spatial resolution
	if (l > 0) {
		//scale in reference to l
		deltaT = defaultSpatialIntegrationLength / l;
	}

	//maxmimum and minimum timesteps
	deltaT = std::max(pSD->minTimeStep(), deltaT);
	deltaT = std::min(pSD->maxTimeStep(), deltaT);

	//newTime
	newTime = lastTime + deltaT;

	Time synctime = MINTIMESTEP * floor((newTime + TIMEERROR) / MINTIMESTEP);
	if (synctime - lastTime > TIMEERROR) {
		//snapping the timestep to the default rhytm at which the other processes synchronize.
		newTime = synctime;
		deltaT = newTime - lastTime;
	}
	// TODO
	//do not run ahead of objects that are calling this.
	//Time wt = pSD->getWallTime(newTime);
	//newTime = wt;
	//deltaT = newTime - lastTime;
}



std::valarray<double> BarberCushmanSolver::rk_step(const std::valarray<double> &X,const std::valarray<double> &w1,const std::valarray<double> &w2,const std::valarray<double> &w3,const std::valarray<double> &w4,const double &De, const double &b, const double &Cmin, const double &r0, const double &v0, const double &Km,const double &Imax ,const double &dr, const std::valarray<double> &IhRK,const std::size_t &k,const std::valarray<double> &r)
{
	std::valarray<double> k_(0.,k);

	double Cea = std::max(0.0, X[0]-Cmin) ;
	double C_0 = X[1] -2*dr/De/b *( Imax*(Cea)/(Km+Cea) -v0*X[0]  );

	k_[0] = w1[0]*C_0+w2[0]*X[0]+w3[0]*X[1]+w4[0]*X[2] - IhRK[0]/b;
	for( std::size_t i=1; i < k-2; ++i) { // because k-1 is the end of the vectors
		k_[i] = w1[i]*X[i-1]+w2[i]*X[i]+w3[i]*X[i+1]+w4[i]*X[i+2] - IhRK[i]/b;
	}
	double C_outer1 = X[k-2] - (2*dr)/(De*b)*r0*v0/r[k-1]*X[k-1];
	double C_outer2 = X[k-1] - (2*dr)/(De*b)*r0*v0/r[k-1]*C_outer1;  // 2bd GHOST POINT --> workaround for last inner upwind point (central differences for ghost)

	k_[k-2] =   w1[k-2]*X[k-3] + w2[k-2]*X[k-2]+w3[k-2]*X[k-1]  +w4[k-2]*C_outer1  - IhRK[k-2]/b;
	k_[k-1] =   w1[k-1]*X[k-2] + w2[k-1]*X[k-1]+w3[k-1]*C_outer1+w4[k-1]*C_outer2  - IhRK[k-1]/b;
	return k_;
}

//***************************************************************************/
std::valarray<double> BarberCushmanSolver::calcIhRK(const std::valarray<double> &X, const std::valarray<double> &S4, const double &Imaxh, const std::valarray<double> &Ah, const double &Km, const double  &Cmin, const std::size_t &k){

	std::valarray<double> Ih(0.,k);

	for(std::size_t i = 0; i < k && Ah[i]>1e-8; ++i){
		double p = X[i]-Km+Cmin-S4[i];
		double z = (p*p)/4 + X[i]*(Km-Cmin) + S4[i]*Cmin ;
		if(z < 0){
			msg::error("calcClrV sqrt of negative value.");
		}
		double Clr   = 0.5*(p)+ sqrt( z ) - Cmin;;
		Ih[i] = Imaxh*Clr/(Km+Clr)*Ah[i];
	}
	return Ih;
}

//***************************************************************************/
double BarberCushmanSolver::L2Norm(const std::valarray<double> &a) {
	double sum;
	std::size_t n = a.size();
	if (n == 0){
		return 0.;
	}
	sum = a[0] * a[0];

	for (std::size_t i = 1; i < n; ++i){
		sum += a[i] * a[i];
	}
	sum = sqrt(sum);
	return sum;
}
//***************************************************************************/

double BarberCushmanSolver::calcCh(const double &C,const double &Km,const double &Cmin,const double &S4){
	return ( ( (   (C-Km+Cmin-S4)+sqrt(  square(-C+Km-Cmin+S4)-4*C*Cmin+4*S4*Cmin+4*C*Km) )/2) - Cmin);
}

double BarberCushmanSolver::calcIh(const double &Ch,const double &Imaxh,const double &Km,const double &Cmin,const double &Ah){
	return  ((Ah*Imaxh*Ch)/(Km+Ch));
}




IntegrationBase * newInstantiationBarberCushmanSolver(){
	return new BarberCushmanSolver();
}

//Register the module
static class AutoRegisterIntegrationFunctionsbcs {
public:
	AutoRegisterIntegrationFunctionsbcs() {
		BaseClassesMap::getIntegrationClasses()["BarberCushmanSolver"] = newInstantiationBarberCushmanSolver;
	}
}p44608510843654654046454385;

