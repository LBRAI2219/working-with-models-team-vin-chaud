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


#include "Solute.hpp"

#include "../../math/BiCGSTAB.hpp"
//#include <iomanip> /* setprecision */

#include "Mapping.hpp"


Solute::Solute(const Watflow & water, const std::string & nameNutrient, const int meshrefinement) :
		d95Solute(nullptr), totSoluteInColumn(nullptr), totSoluteInColumnDissolved(nullptr),
		totSoluteInColumnAbsorbed(nullptr), totSoluteChange(nullptr), massBalanceError_(nullptr),
		totSink_(nullptr), sumTopBoundaryFlux_(nullptr), sumBottomBoundaryFlux_(nullptr),
		totalMineralization_(nullptr), sPrec_(nullptr), sinterception_(nullptr), scPrec_(nullptr),
		water_(water), factor(meshrefinement),
		soluteMesh(factor), mineral(nullptr), // for example mineralization on solute Mesh because it interacts with cSink
		NumNP(soluteMesh.getNumNP()), NumEl(soluteMesh.getNumEl()),
		NumBP(soluteMesh.getNumBP()), x(soluteMesh.getCordX()), y(soluteMesh.getCordY()), z(soluteMesh.getCordZ()),
		soluteMatrix(NumNP), PeCr(5.0), /* TODO what is a good criterion? */
		PeCrMx(0.), Peclet(0.), Courant(0.), dtMaxC(1.e+30), Conc(0., NumNP),  // this is set and used in this class
		cSink(0., NumNP), KXB(soluteMesh.getKXB()), // NumBP size
		Dispxx(0., NumNP), Dispyy(0., NumNP), Dispzz(0., NumNP), Dispxy(0., NumNP), Dispxz(0., NumNP), Dispyz(0., NumNP),
		ListNE(soluteMesh.getListNE()), bulkDensity(0., NumNP), satDiffusionCoeff(0., NumNP),
		longitudinalDispersivity(0., NumNP), transverseDispersivity(0., NumNP), adsorptionCoeff(0., NumNP), nutrient(nameNutrient),
		//these are currently not used at all.
		/*mu_w(0., NumNP),     First-order rate constant for dissolved phase
		 mu_s(0., NumNP),     First-order rate constant for solid phase
		 gamma_w(0., NumNP),  Zero-order rate constant for dissolved phase
		 gamma_s(0., NumNP),   Zero-order rate constant for solid phase      */
		seedReserves(0), ConVolI(0.), matchmode(1), Nsub(soluteMesh.getNsubel()), fine(factor > 1 ? true : false) {

	// need matrix structure
	soluteMesh.easy_init_matrix(soluteMatrix);

	//mineralization model
	if (nutrient == "nitrate") {
		mineral = new Mineralization(soluteMesh);
	}

	double Y;

	//set d95 solute container
	SimulaBase* p = ORIGIN->getChild("soil")->existingChild(nutrient);
	if (p) {
		SimulaBase *d = p->existingChild("D90", "cm");
		if (d)
			d95Solute = dynamic_cast<SimulaTable<Time>*>(d);

		d = p->existingChild("totalSoluteInColumn", "umol");
		if (d)
			totSoluteInColumn = dynamic_cast<SimulaTable<Time>*>(d);
		d = p->existingChild("totalDissolvedSoluteInColumn", "umol");
		if (d)
			totSoluteInColumnDissolved = dynamic_cast<SimulaTable<Time>*>(d);
		d = p->existingChild("totalAbsorbedSoluteInColumn", "umol");
		if (d)
			totSoluteInColumnAbsorbed = dynamic_cast<SimulaTable<Time>*>(d);

		//todo come up with better names for these
		d = p->existingChild("totalSoluteChange", "umol");
		if (d)
			totSoluteChange = dynamic_cast<SimulaTable<Time>*>(d);

		d = p->existingChild("topBoundaryFluxRate", "umol/day");
		if (d)
			sumTopBoundaryFlux_ = dynamic_cast<SimulaTable<Time>*>(d);
		d = p->existingChild("bottomBoundaryFluxRate", "umol/day");
		if (d)
			sumBottomBoundaryFlux_ = dynamic_cast<SimulaTable<Time>*>(d);

		d = p->existingChild("totalMineralizationRate", "umol/day");
		if (d)
			totalMineralization_ = dynamic_cast<SimulaTable<Time>*>(d);

		d = p->existingChild("massBalanceError", "umol");
		if (d)
			massBalanceError_ = dynamic_cast<SimulaTable<Time>*>(d);

		d = p->existingChild("totalSinkRate", "umol/day");
		if (d)
			totSink_ = dynamic_cast<SimulaTable<Time>*>(d);

	}

	totalNutrientUptake_ = ORIGIN->existingPath("/plants/nitrate/plantNutrientUptake", "umol"); //"umol"
	//todo for backwards compatibility, simply insert this, instead of not calculating the balance
	if (totalNutrientUptake_)
		totalNutrientUptake_->get(0, seedReserves);	//todo these are a problem, when plants are not planted at time 0
	sPrec_ = ORIGIN->existingPath("/environment/atmosphere/precipitation", "cm/day");
	sinterception_ = ORIGIN->existingPath("/atmosphere/interception", "cm/day");
	scPrec_ = ORIGIN->existingPath("/atmosphere/concentrationInPrecipitation", "umol/ml");

	p = ORIGIN->getPath("environment/soil/" + nutrient + "/concentration");
	for (tIndex i = 0; i != NumNP; ++i) {
		p->get(0, Coordinate(x[i], z[i], y[i]), Y);
		Conc[i] = Y;
	}

	p = ORIGIN->getPath("environment/soil/bulkDensity");
	for (tIndex i = 0; i != NumNP; ++i) {
		p->get(0, Coordinate(x[i], z[i], y[i]), Y);
		bulkDensity[i] = Y;
	}
	// Ionic or molecular diffusion coefficient in free water
	p = ORIGIN->getPath("environment/soil/" + nutrient + "/saturatedDiffusionCoefficient");
	for (tIndex i = 0; i != NumNP; ++i) {
		p->get(0, Coordinate(x[i], z[i], y[i]), Y);
		satDiffusionCoeff[i] = Y*(double)soluteMesh.getpotmask()[i];
	}

	p = ORIGIN->getPath("environment/soil/" + nutrient + "/longitudinalDispersivity");
	for (tIndex i = 0; i != NumNP; ++i) {
		p->get(0, Coordinate(x[i], z[i], y[i]), Y);
		longitudinalDispersivity[i] = Y;
	}

	p = ORIGIN->getPath("environment/soil/" + nutrient + "/transverseDispersivity");
	for (tIndex i = 0; i != NumNP; ++i) {
		p->get(0, Coordinate(x[i], z[i], y[i]), Y);
		transverseDispersivity[i] = Y;
	}

	p = ORIGIN->getPath("environment/soil/" + nutrient + "/adsorptionCoefficient");
	for (tIndex i = 0; i != NumNP; ++i) {
		p->get(0, Coordinate(x[i], z[i], y[i]), Y);
		adsorptionCoeff[i] = Y;
	}

	//try to determine initial timestep based on peclet courant numbers
	double dt = 0.001; //just an initial guess to work with

	if (fine) {
		tVector hNew(0., NumNP), Con(0., NumNP), thOld(0., NumNP);

		Mapping::interpolate(water_.getMesh(), water_.gethNew(), hNew, factor);
		Mapping::interpolate(water_.getMesh(), water_.getCon(), Con, factor);
		Mapping::interpolate(water_.getMesh(), water_.getthOld(), thOld, factor);

		tVector Vx(0., NumNP), Vy(0., NumNP), Vz(0., NumNP);
		velocities(hNew, Con, Vx, Vy, Vz);
		dispersion_coeff(dt, thOld, Vx, Vy, Vz);
		peclet_courant(dt, thOld, Vx, Vy, Vz);
	} else {
		tVector Vx(0., NumNP), Vy(0., NumNP), Vz(0., NumNP);
		velocities(water_.gethNew(), water_.getCon(), Vx, Vy, Vz); // set Vx, Vy, Vz
		dispersion_coeff(dt, water_.getthOld(), Vx, Vy, Vz); // set needed in Peclet_courant, and upstream_weighing_factors
		peclet_courant(dt, water_.getthOld(), Vx, Vy, Vz); //used for initial timestep
	}

	// find out what method to use for matching roots to fem nodes
	p = ORIGIN->existingPath("environment/dimensions/rootMatchingMode");
	if (p) {
		std::string mode;
		p->get(mode);
		if (mode == "postma and lynch, 2011") {
			matchmode = 1;
		} else if (mode == "ignore root placement") {
			matchmode = 2;
			msg::warning("Solute: using rootMatchingMode: 'ignore root placement' as set in environment/dimensions/rootMatchingMode");
		} else {
			msg::error("Solute: unknown string for environment/dimensions/rootMatchingMode ");
		}
	}

	//report initial values for solid in the column
	this->reportTotSolidColumn(0);
	tVector B(0., NumNP);
	this->cbound_condition(0., B);
	if (totSink_)
		totSink_->set(0, 0); //normally should be zero, will be reset as soon as we do a timestep,
}

//time t at the start of the timestep and the size of the timestep
void Solute::calculate(const double & t, const double & dt) {

	// if no interpolation, these are empty
	tVector hNew;
	tVector Con;
	tVector thNew;
	tVector thOld;
	tVector ConOld;
	tVector hOld;
	if (fine) {
		hNew.resize(NumNP);
		Con.resize(NumNP);
		thNew.resize(NumNP);
		Mapping::interpolate(water_.getMesh(), water_.gethNew(), hNew, factor);
		Mapping::interpolate(water_.getMesh(), water_.getCon(), Con, factor);
		Mapping::interpolate(water_.getMesh(), water_.getthNew(), thNew, factor);
	}

	this->setSoluteSink(t, t + dt, cSink);

	// (Re-)Initialisation
	tVector B(0., NumNP);
	soluteMatrix.resetToZero(); // set soluteMatrix to zero with the right entries

	tVector DS(0., NumNP);
	tVector Hhdt(0), Conhdt(0), Thhdt(0); //dummy vectors for half timesteps

	// NLEVEL is 1 or if epsi<1 it is 2
	// epsi should be 0.5 or 1 for 0=explicit method, 0.5 heuns implicit, 1 fully implicit method
	double incr = 1. / (double) NLEVEL;
	double eps = 0.;
	for (tIndex Level = 1; Level <= NLEVEL; ++Level) { // this for loop, loops only twice or once. NLEVEL\in{1,2}
		// intermediate time
		eps += incr;
		// set pointers to either the values at beginning of the timestep, or the end of the timestep.
		const tVector * hhead = nullptr;
		const tVector * hconductivity = nullptr;
		const tVector * theta = nullptr;
		if (Level == NLEVEL) {
			if (fine) {
				hhead = &hNew;
				hconductivity = &Con;
				theta = &thNew;
			} else {
				hhead = &water_.gethNew();
				hconductivity = &water_.getCon();
				theta = &water_.getthNew();
			}
		} else {
			double alf = 1. - eps;
			if (fine) {
				thOld.resize(NumNP);
				ConOld.resize(NumNP);
				hOld.resize(NumNP);
				Mapping::interpolate(water_.getMesh(), water_.getthOld(), thOld, factor);
				Mapping::interpolate(water_.getMesh(), water_.getthNew(), thNew, factor);
				Thhdt.resize(NumNP);
				Thhdt = thOld * eps + thNew * alf;
				Hhdt.resize(NumNP);
				Mapping::interpolate(water_.getMesh(), water_.gethOld(), hOld, factor);
				Mapping::interpolate(water_.getMesh(), water_.gethNew(), hNew, factor);
				Hhdt = hOld * eps + hNew * alf;
				Conhdt.resize(NumNP);
				Mapping::interpolate(water_.getMesh(), water_.getConOld(), ConOld, factor);
				Mapping::interpolate(water_.getMesh(), water_.getCon(), Con, factor);
				Conhdt = ConOld * eps + Con * alf;
				hhead = &Hhdt;
				hconductivity = &Conhdt;
				theta = &Thhdt;
			} else {
				//average values
				Thhdt.resize(NumNP);
				Thhdt = water_.getthOld() * eps + water_.getthNew() * alf;
				Hhdt.resize(NumNP);
				Hhdt = water_.gethOld() * eps + water_.gethNew() * alf;
				Conhdt.resize(NumNP);
				Conhdt = (water_.getConOld() * eps + water_.getCon() * alf);
				hhead = &Hhdt;
				hconductivity = &Conhdt;
				theta = &Thhdt;
			}
		}

		// parameter g and f
		// gc is in mass/volume/time
		//currently gamma_w and gamma_s are set to 0
		//tVector Gc = gamma_w * (*theta) + bulkDensity * gamma_s - cSink;
		const tVector Gc(cSink * -1.);

		// fc is in 1/time (sink is in per time) - i think this is vol water / element vol / time
		//const tVector &Sink(water_.getwSink()); // get Sink from water
		//currently mu_w and mu_s are 0
		//tVector Fc = mu_w * (*theta) + bulkDensity * mu_s * adsorptionCoeff + Sink;

		tVector wSink;
		if (fine) {
			wSink.resize(NumNP);
			Mapping::interpolate(water_.getMesh(), water_.getwSink(), wSink, factor);
		}
		const tVector & Fc = (fine) ? wSink : (water_.getwSink());

		// adsorption= bulkDensity* adsorptionCoeff
		// used to have (ThOld[i]*alf+ThNew[i]*epsi)
		const tVector Ac = -(*theta) - bulkDensity * adsorptionCoeff;

		DS = 0;

		// summing up of root extraction terms,
		// note this is purely to compute a balance and is further not used in the routine
		// double tmpr=CumChR;

		/* todo: this reduces simply to rootch=sum(cSink*vol*dt) which is the total nutrient uptake. There is no point in this check
		 * unless gamma_s, gamma_w, mu_s or mu_w are not 0.
		 */

		/* tIndex count = 0; tIndex it_subelem = 0;
		 *
		 for (list_const_iter it = soluteMesh.getSubelementlist().begin(); it != soluteMesh.getSubelementlist().end(); ++it) { // Loop over sub-elements
		 const tIndex i = (*it)[0], j = (*it)[1], k = (*it)[2], l = (*it)[3];

		 // volume * time for integrating over space and time
		 const double VE = dt * incr * soluteMesh.getSubelemetVolumes()[it_subelem]; //volumeTetraHedron(NumNP,i,j,k,l,x,y,z)
		 // uptake as determined by simroot
		 // root ch is in mass - csink is mass/vol/time : vol * vol * time * mass/vol/time / vol
		 const double RootCh ( VE * (cSink[i] * global_weight[count]
		 + cSink[j] * global_weight[count+1]
		 + cSink[k] * global_weight[count+2]
		 + cSink[l] * global_weight[count+3]));
		 // Gc is in mass/volume/time -- CumCh0 is in mass
		 CumCh0 += -VE * (Gc[i] * global_weight[count] + Gc[j] * global_weight[count+1]
		 + Gc[k] * global_weight[count+2] + Gc[l] * global_weight[count+3]) - RootCh;
		 // fc is in vol/vol/time -- cumch1 is in mass
		 // we are substracting the sink term which we previously added.
		 CumCh1 += -VE * (( (Fc[i] - Sink[i]) * Conc[i] * global_weight[count])
		 + ( (Fc[j] - Sink[j]) * Conc[j] * global_weight[count+1])
		 + ( (Fc[k] - Sink[k]) * Conc[k] * global_weight[count+2])
		 + ( (Fc[l] - Sink[l]) * Conc[l] * global_weight[count+3]) );

		 CumChR += RootCh;
		 count += 4; // increase counter to get the next four (nodal) weights
		 ++it_subelem;
		 }*/

		//calculate velocity vx,vy,vz
		tVector Vx(0., NumNP), Vy(0., NumNP), Vz(0., NumNP);
		velocities(*hhead, *hconductivity, Vx, Vy, Vz);
		//dispersion & diffusion is calculated from the velocity
		dispersion_coeff(dt, (*theta), Vx, Vy, Vz);
		//computes the maximum local Peclet and Courant numbers and the maximum permissible time step
		peclet_courant(dt, (*theta), Vx, Vy, Vz);
		//weighing factors for all sides of all elements.
		//effects matrix WeTab
#ifdef lUpW
		tVector WeTab0(Nsub * NumEl), WeTab1(Nsub * NumEl), WeTab2(Nsub * NumEl), WeTab3(Nsub * NumEl), WeTab4(Nsub * NumEl), WeTab5(
				Nsub * NumEl);
		upstream_weighing_factors(0, Vx, Vy, Vz, WeTab0);
		upstream_weighing_factors(0, Vx, Vy, Vz, WeTab1);
		upstream_weighing_factors(0, Vx, Vy, Vz, WeTab2);
		upstream_weighing_factors(0, Vx, Vy, Vz, WeTab3);
		upstream_weighing_factors(0, Vx, Vy, Vz, WeTab4);
		upstream_weighing_factors(0, Vx, Vy, Vz, WeTab5);
#else
		// no artificial dispersion and no upward scaling
		const tVector DPom(dt / 6. / ((*theta) + bulkDensity * adsorptionCoeff));			//elementwise
		Dispxx += Vx * Vx * DPom;
		Dispyy += Vy * Vy * DPom;
		Dispzz += Vz * Vz * DPom;
		Dispxy += Vx * Vy * DPom;
		Dispxz += Vx * Vz * DPom;
		Dispyz += Vy * Vz * DPom;
#endif

		tVector F(0., NumNP);

		//       Loop on elements
		// in our case these are cm cubes

		tIndex numSEl = 0;
		unsigned int count = 0;
		for (list_const_iter it = soluteMesh.getSubelementlist().begin(); it != soluteMesh.getSubelementlist().end(); ++it) { // Loop over sub-elements
			//index values of the nodes of this tetrahedal ellement
			const tIndex i((*it)[0]), j((*it)[1]), k((*it)[2]), l((*it)[3]);

			// The volume of a tetrahedron is equal to the determinant formed by
			// writing the coordinates of the vertices as columns and
			// then appending a row of ones along the bottom.

			const double VE(soluteMesh.getSubelemetVolumes()[numSEl]);
			const double Det(VE * 6.);

			// Calculate Velocities - these are simply differences in hhead in x,y and z direction divided
			// by the determinant which is either the volume of the cube or twice the volume of the cube

			const tVector & Bi(soluteMesh.getBi()), &Ci(soluteMesh.getCi()), &Di(soluteMesh.getDi());
			const double Vxx((Bi[count] * (*hhead)[i] + Bi[count + 1] * (*hhead)[j] + Bi[count + 2] * (*hhead)[k] + Bi[count + 3] * (*hhead)[l]) / Det);
			const double Vyy((Ci[count] * (*hhead)[i] + Ci[count + 1] * (*hhead)[j] + Ci[count + 2] * (*hhead)[k] + Ci[count + 3] * (*hhead)[l]) / Det);
			const double Vzz(((Di[count] * (*hhead)[i] + Di[count + 1] * (*hhead)[j] + Di[count + 2] * (*hhead)[k] + Di[count + 3] * (*hhead)[l]) / Det) + 1); // Kazz=1, det=6*Ve, see eq 4.25 manual, assuming no further anisotropy

			// save velocity based on nodal values
			std::array<double, 4> VxE, VyE, VzE;
			VxE[0] = -(*hconductivity)[i] * Vxx;
			VyE[0] = -(*hconductivity)[i] * Vyy;
			VzE[0] = -(*hconductivity)[i] * Vzz;
			VxE[1] = -(*hconductivity)[j] * Vxx;
			VyE[1] = -(*hconductivity)[j] * Vyy;
			VzE[1] = -(*hconductivity)[j] * Vzz;
			VxE[2] = -(*hconductivity)[k] * Vxx;
			VyE[2] = -(*hconductivity)[k] * Vyy;
			VzE[2] = -(*hconductivity)[k] * Vzz;
			VxE[3] = -(*hconductivity)[l] * Vxx;
			VyE[3] = -(*hconductivity)[l] * Vyy;
			VzE[3] = -(*hconductivity)[l] * Vzz;

			// calculate velocity based on elemental values
			// the vol weighing here is probably at par with the linear assumption above, but probably not worse than 1/4 weighing
			//average conductivity in the element
			const double ConE(0.25 * ((*hconductivity)[i] + (*hconductivity)[j] + (*hconductivity)[k] + (*hconductivity)[l]));

			const double VxEE = -ConE * Vxx;
			const double VyEE = -ConE * Vyy;
			const double VzEE = -ConE * Vzz;

			// Diffusion & dispersion?
			const double FMul(VE / 5.);  // this is a weight, for GcE and F , see Manual Equation 5.8
			const double SMul1 = -1. / VE / 36.;
			const double SMul2 = VE / 30.;
			const double GcE = 0.25 * (Gc[i] + Gc[j] + Gc[k] + Gc[l]);
			const double FcE = 0.25 * (Fc[i] + Fc[j] + Fc[k] + Fc[l]);
			const double AcE = Level == NLEVEL ? (Ac[i] + Ac[j] + Ac[k] + Ac[l]) * 0.25 : 0.;
			const double EC1 = 0.25 * (Dispxx[i] + Dispxx[j] + Dispxx[k] + Dispxx[l]);
			const double EC2 = 0.25 * (Dispyy[i] + Dispyy[j] + Dispyy[k] + Dispyy[l]);
			const double EC3 = 0.25 * (Dispzz[i] + Dispzz[j] + Dispzz[k] + Dispzz[l]);
			const double EC4 = 0.25 * (Dispxy[i] + Dispxy[j] + Dispxy[k] + Dispxy[l]);
			const double EC5 = 0.25 * (Dispxz[i] + Dispxz[j] + Dispxz[k] + Dispxz[l]);
			const double EC6 = 0.25 * (Dispyz[i] + Dispyz[j] + Dispyz[k] + Dispyz[l]);

#ifdef lUpW
			//upwind
			//todo this could at the cost of memory be expanded by storing the A arrays.
			std::array<double, 4> Wx = {0., 0., 0., 0.}, Wy = {0., 0., 0., 0.}, Wz = {0., 0., 0., 0.};
			const double W12 = WeTab0[numSEl];
			const double W13 = WeTab1[numSEl];
			const double W14 = WeTab2[numSEl];
			const double W23 = WeTab3[numSEl];
			const double W24 = WeTab4[numSEl];
			const double W34 = WeTab5[numSEl];
			const double A11 = -2. * W12 + 2. * W14 + 2. * W13;
			const double A12 = -2. * W12 + W14 + W13;
			const double A13 = -W12 + W14 + 2. * W13;
			const double A14 = -W12 + 2. * W14 + W13;
			const double A21 = -W23 + 2. * W12 + W24;
			const double A22 = -2. * W23 + 2. * W12 + 2. * W24;
			const double A23 = -2. * W23 + W12 + W24;
			const double A24 = -W23 + W12 + 2. * W24;
			const double A31 = -W34 + W23 - 2. * W13;
			const double A32 = -W34 + 2. * W23 - W13;
			const double A33 = -2. * W34 + 2. * W23 - 2. * W13;
			const double A34 = -2. * W34 + W23 - W13;
			const double A41 = -2. * W14 + W34 - W24;
			const double A42 = -W14 + W34 - 2. * W24;
			const double A43 = -W14 + 2. * W34 - W24;
			const double A44 = -2. * W14 + 2. * W34 - 2. * W24;
			Wx[0] = VxE[0] * A11 + VxE[1] * A12 + VxE[2] * A13 + VxE[3] * A14;
			Wx[1] = VxE[0] * A21 + VxE[1] * A22 + VxE[2] * A23 + VxE[3] * A24;
			Wx[2] = VxE[0] * A31 + VxE[1] * A32 + VxE[2] * A33 + VxE[3] * A34;
			Wx[3] = VxE[0] * A41 + VxE[1] * A42 + VxE[2] * A43 + VxE[3] * A44;
			Wy[0] = VyE[0] * A11 + VyE[1] * A12 + VyE[2] * A13 + VyE[3] * A14;
			Wy[1] = VyE[0] * A21 + VyE[1] * A22 + VyE[2] * A23 + VyE[3] * A24;
			Wy[2] = VyE[0] * A31 + VyE[1] * A32 + VyE[2] * A33 + VyE[3] * A34;
			Wy[3] = VyE[0] * A41 + VyE[1] * A42 + VyE[2] * A43 + VyE[3] * A44;
			Wz[0] = VzE[0] * A11 + VzE[1] * A12 + VzE[2] * A13 + VzE[3] * A14;
			Wz[1] = VzE[0] * A21 + VzE[1] * A22 + VzE[2] * A23 + VzE[3] * A24;
			Wz[2] = VzE[0] * A31 + VzE[1] * A32 + VzE[2] * A33 + VzE[3] * A34;
			Wz[3] = VzE[0] * A41 + VzE[1] * A42 + VzE[2] * A43 + VzE[3] * A44;
#endif

			// building the matrix looping over all combinations
			for (tIndex j1 = 0; j1 < 4; ++j1) { // here we simply looping through the node points of the element

				const tIndex i1 = (*it)[j1];
				F[i1] += FMul * (GcE + Gc[i1] / 4.);
				if (Level == NLEVEL) {
					DS[i1] += FMul * (AcE + Ac[i1] / 4.);
				}
				for (tIndex j2 = 0; j2 < 4; ++j2) { //here we simply looping through the node points of the element
					const tIndex i2 = (*it)[j2];
					double m_ij = SMul1
							* (EC1 * Bi[count + j1] * Bi[count + j2] + EC2 * Ci[count + j1] * Ci[count + j2] + EC3 * Di[count + j1] * Di[count + j2]
									+ EC4 * (Bi[count + j1] * Ci[count + j2] + Ci[count + j1] * Bi[count + j2])
									+ EC5 * (Bi[count + j1] * Di[count + j2] + Di[count + j1] * Bi[count + j2])
									+ EC6 * (Ci[count + j1] * Di[count + j2] + Di[count + j1] * Ci[count + j2]));
					m_ij -= (Bi[count + j2] / 30. * (VxEE + VxE[j1] / 4.) + Ci[count + j2] / 30. * (VyEE + VyE[j1] / 4.)
							+ Di[count + j2] / 30. * (VzEE + VzE[j1] / 4.));
#ifdef lUpW
					m_ij -= (Bi[j2] / 240. * Wx[j1] + Ci[j2] / 240. * Wy[j1] + Di[j2] / 240. * Wz[j1]);
#endif

					if (i1 == i2) {
						m_ij += SMul2 * 2. * (FcE + (Fc[i1] + Fc[i2]) / 4.);
					} else {
						m_ij += SMul2 * 1. * (FcE + (Fc[i1] + Fc[i2]) / 4.);
					}

					if (Level != NLEVEL) {
						B[i1] -= incr * m_ij * Conc[i2];
					} else {
						soluteMatrix.addValueUnsafely(i1, i2, 1. / (double) NLEVEL * m_ij);
					}

				}
			}
			count += 4;
			++numSEl;
		} // tetrahedra loop (sub elements)

		if (Level != NLEVEL) {
			B -= incr * F;
		} else {
			B += DS / dt * Conc - (1. / (double) NLEVEL) * F;
		}
	} // end of level loop, which is the loop around the integration method (estimate, and redo in this case)

	for (tIndex i = 0; i < NumNP; ++i) {
		soluteMatrix.addValueUnsafely(i, i, DS[i] / dt);
	}

	// Boundary condition
	cbound_condition(t, B);

	// Debug only: Look at first matrix
//	    std::ofstream smatr;
//	    smatr << std::scientific;
//	    smatr.open("soluteMatrix.txt");
//	    soluteMatrix.print_sparse(smatr);
//	    smatr.close();
//	    exit(0);

	BiCGSTAB::solve(soluteMatrix, B, Conc, 1.e-8, NumNP, NO_PRECOND); // conc is the solution, B is rhs

	// the SolInf calculations for output
	double t1 = t + dt;
	this->reportTotSolidColumn(t1);
	this->setrootsurfacevalues(t1);		//new values are for the end of the timestep

	return;
}

void Solute::peclet_courant(const double & dt, const tVector & theta_, const tVector & Vx, const tVector & Vy, const tVector & Vz) {

	double xmax, xmin, ymax, ymin, zmax, zmin, delX, delY, delZ;
	double DxE, DyE, DzE, VxE, VyE, VzE;
	double PecX, PecY, PecZ, VxMax, VyMax, VzMax;
	double R1, R2, R3, R4, RMin, dt1, dt2, dt3, CourX, CourY, CourZ, Cour1, Cour2, Cour3;

	// Loop over sub elements
	for (list_const_iter it = soluteMesh.getSubelementlist().begin(); it != soluteMesh.getSubelementlist().end(); ++it) { // Loop over sub-elements
		//index values of the nodes of this tetrahedral element
		const tIndex i((*it)[0]), j((*it)[1]), k((*it)[2]), l((*it)[3]);

		xmax = std::max( { x[i], x[j], x[k], x[l] });
		xmin = std::min( { x[i], x[j], x[k], x[l] });
		ymax = std::max( { y[i], y[j], y[k], y[l] });
		ymin = std::min( { y[i], y[j], y[k], y[l] });
		zmax = std::max( { z[i], z[j], z[k], z[l] });
		zmin = std::min( { z[i], z[j], z[k], z[l] });
		delX = xmax - xmin;
		delY = ymax - ymin;
		delZ = zmax - zmin;
		DxE = std::max(0.25 * (Dispxx[i] + Dispxx[j] + Dispxx[k] + Dispxx[l]), 1.e-30);
		DyE = std::max(0.25 * (Dispyy[i] + Dispyy[j] + Dispyy[k] + Dispyy[l]), 1.e-30);
		DzE = std::max(0.25 * (Dispzz[i] + Dispzz[j] + Dispzz[k] + Dispzz[l]), 1.e-30);
		VxE = std::fabs(0.25 * (Vx[i] + Vx[j] + Vx[k] + Vx[l]));
		VyE = std::fabs(0.25 * (Vy[i] + Vy[j] + Vy[k] + Vy[l]));
		VzE = std::fabs(0.25 * (Vz[i] + Vz[j] + Vz[k] + Vz[l]));
		VxMax = std::max( { std::fabs(Vx[i]) / theta_[i], std::fabs(Vx[j]) / theta_[j], std::fabs(Vx[k]) / theta_[k], std::fabs(Vx[l]) / theta_[l] });
		VyMax = std::max( { std::fabs(Vy[i]) / theta_[i], std::fabs(Vy[j]) / theta_[j], std::fabs(Vy[k]) / theta_[k], std::fabs(Vy[l]) / theta_[l] });
		VzMax = std::max( { std::fabs(Vz[i]) / theta_[i], std::fabs(Vz[j]) / theta_[j], std::fabs(Vz[k]) / theta_[k], std::fabs(Vz[l]) / theta_[l] });

		PecX = VxE * delX / DxE;
		PecY = VyE * delY / DyE;
		PecZ = VzE * delZ / DzE;
		Peclet = std::max( { Peclet, PecX, PecY, PecZ });
		Peclet = std::min(Peclet, 99999.);

		R1 = 1. + bulkDensity[i] * adsorptionCoeff[i] / theta_[i];
		R2 = 1. + bulkDensity[j] * adsorptionCoeff[j] / theta_[j];
		R3 = 1. + bulkDensity[k] * adsorptionCoeff[k] / theta_[k];
		R4 = 1. + bulkDensity[l] * adsorptionCoeff[l] / theta_[l];
		RMin = std::min( { R1, R2, R3, R4 });

		dt1 = 1.e+30;
		dt2 = 1.e+30;
		dt3 = 1.e+30;
		CourX = VxMax * dt / delX / RMin;
		CourY = VyMax * dt / delY / RMin;
		CourZ = VzMax * dt / delZ / RMin;
		Courant = std::max( { Courant, CourX, CourY, CourZ });

#if defined(lUpW) || defined(lArtD)
		Cour1 = 1.0;
		Cour2 = 1.0;
		Cour3 = 1.0;
#else
		Cour1 = std::min(1., PeCr / std::max(0.5, PecX));
		Cour2 = std::min(1., PeCr / std::max(0.5, PecY));
		Cour3 = std::min(1., PeCr / std::max(0.5, PecZ));
#endif
		if (VxMax > 0.) {
			dt1 = Cour1 * delX * RMin / VxMax;
		}
		if (VyMax > 0.) {
			dt2 = Cour2 * delY * RMin / VyMax;
		}
		if (VzMax > 0.) {
			dt3 = Cour3 * delZ * RMin / VzMax;
		}
		dtMaxC = std::min( { dtMaxC, dt1, dt2, dt3 });

		PeCrMx = std::max( { PeCrMx, PecX * CourX, PecY * CourY, PecZ * CourZ });

	}
	return;
}

#ifdef lUpW

void Solute::upstream_weighing_factors(const int j, const tVector & Vx, const tVector & Vy, const tVector & Vz, tVector &WeTab) const {
	tIndex M1(0), M2(0), M3(0), M4(0), NUS(Nsub);
	std::array<tIndex, 4> List;
	double VV, Vxx, Vyy, Vzz, A, B, C, Aleng, cosx, cosy, cosz, cosDel, Delta;
	double Val, Vel, Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, DAL, Disp, aa;

	tIndex numSEl = 0; //  number of subelements thetrahedrals
	// Loop over sub elements
	for (list_const_iter it = soluteMesh.getSubelementlist().begin(); it != soluteMesh.getSubelementlist().end(); ++it) { // Loop over sub-elements
		M1 = (*it)[0];
		M2 = (*it)[1];
		M3 = (*it)[2];
		M4 = (*it)[3];

		List[0] = M1;
		List[1] = M2;
		List[2] = M3;
		List[3] = M4;

		WeTab[numSEl] = 0.;

		switch (j) {
			case 0: // TODO nothing to do in this case?
			M1 = List[0];
			M2 = List[1];
			break;
			case 1:
			M1 = List[2];
			M2 = List[0];
			break;
			case 2:
			M1 = List[3];
			M2 = List[0];
			break;
			case 3:
			M1 = List[1];
			M2 = List[2];
			break;
			case 4:
			M1 = List[3];
			M2 = List[1];
			break;
			case 5:
			M1 = List[2];
			M2 = List[3];
			break;
		}

		A = x[M2] - x[M1];
		B = y[M2] - y[M1];
		C = z[M2] - z[M1];
		Aleng = std::sqrt(A * A + B * B + C * C);
		cosx = 1. / Aleng * A;
		cosy = 1. / Aleng * B;
		cosz = 1. / Aleng * C;
		Vxx = (Vx[M1] + Vx[M2]) * 0.5;
		Vyy = (Vy[M1] + Vy[M2]) * 0.5;
		Vzz = (Vz[M1] + Vz[M2]) * 0.5;
		VV = std::sqrt(Vxx * Vxx + Vyy * Vyy + Vzz * Vzz);

		if (std::fabs(VV) < 1.e-30) {
			continue;
		}

		cosDel = 1. / VV / Aleng * (A * Vxx + B * Vyy + C * Vzz);
		Delta = std::acos(std::fabs(cosDel));

		if (NUS != Nsub && Delta > 3.1415926536) {
			continue;
		}

		Val = Vxx * cosx + Vyy * cosy + Vzz * cosz;
		Dxx = (Dispxx[M1] + Dispxx[M2]) * 0.5;
		Dyy = (Dispyy[M1] + Dispyy[M2]) * 0.5;
		Dzz = (Dispzz[M1] + Dispzz[M2]) * 0.5;
		Dxy = (Dispxy[M1] + Dispxy[M2]) * 0.5;
		Dxz = (Dispxz[M1] + Dispxz[M2]) * 0.5;
		Dyz = (Dispyz[M1] + Dispyz[M2]) * 0.5;
		DAL = std::fabs(Dxx * cosx * cosx + Dyy * cosy * cosy + Dzz * cosz * cosz + 2. * Dxy * cosx * cosy + 2. * Dxz * cosx * cosz + 2. * Dyz * cosy * cosz);
		Vel = Val * Aleng;
		Disp = 2.0 * DAL;
		aa = 11.;
		if (Disp > 0.) {
			aa = std::fabs(Vel / Disp);
		}
		if (Disp < 1.e-30 || std::fabs(Vel) < 0.001 * VV || std::fabs(aa) > 10.) {
			if (std::fabs(Vel) < 0.001 * VV) {
				WeTab[numSEl] = 0.0;
			}
			if (Vel > 0.001 * VV) {
				WeTab[numSEl] = 1.0;
			}
			if (Vel < -0.001 * VV) {
				WeTab[numSEl] = -1.0;
			}
		} else {
			WeTab[numSEl] = 1.0 / std::tanh(Vel / Disp) - Disp / Vel;
		}

		++numSEl;
	}

	return;
}
#endif

void Solute::dispersion_coeff(const double & dt, const tVector & theta_, const tVector & Vx, const tVector & Vy, const tVector & Vz) {
//	    tVector Tau(NumNP);
	double DispL, DispT, Vabs;

//		for (tIndex i = 0; i < NumNP; ++i) {
//			// tortuosity factor according to Millington and Quirk [1961]
//			Tau[i] = std::pow(theta_[i], (7. / 3.)) / (thSat[i] * thSat[i]);
//			// diffusion = theta * dif in water * Darcian fluid flux density * tau
//			Dif[i] = theta_[i] * satDiffusionCoeff[i] * Tau[i];
//		}

	tVector thSat;
	if (fine) {
		thSat.resize(NumNP);
		Mapping::interpolate(water_.getMesh(), water_.getthSat(), thSat, factor);
	}
	const tVector Dif =
			fine ? (theta_ * satDiffusionCoeff * std::pow(theta_, (7. / 3.)) / std::pow(thSat, 2.)) : (theta_ * satDiffusionCoeff * std::pow(theta_, (7. / 3.))
							/ std::pow(water_.getthSat(), 2.)); // component wise

	// todo this can now be pulled apart into array operations?
	for (tIndex i = 0; i < NumNP; ++i) {
		// dispersion
		DispL = longitudinalDispersivity[i];
		DispT = transverseDispersivity[i];
		Vabs = std::sqrt(Vx[i] * Vx[i] + Vy[i] * Vy[i] + Vz[i] * Vz[i]);
#ifdef lArtD
		if (Vabs > 0.) {
			// adjust dispersion so the Peclet number criterium is met
			DispL = std::max(DispL, Vabs * dt / (theta_[i] + bulkDensity[i] * adsorptionCoeff[i]) / PeCr - Dif[i] / Vabs);
		}
#endif
		Dispxx[i] = Dif[i];
		Dispyy[i] = Dif[i];
		Dispzz[i] = Dif[i];
		Dispxy[i] = 0.;
		Dispxz[i] = 0.;
		Dispyz[i] = 0.;
		if (Vabs > 0.) {
			// velocity is greater than 0
			Dispxx[i] = DispL * Vx[i] * Vx[i] / Vabs + DispT * Vy[i] * Vy[i] / Vabs + DispT * Vz[i] * Vz[i] / Vabs + Dif[i];
			Dispyy[i] = DispL * Vy[i] * Vy[i] / Vabs + DispT * Vx[i] * Vx[i] / Vabs + DispT * Vz[i] * Vz[i] / Vabs + Dif[i];
			Dispzz[i] = DispL * Vz[i] * Vz[i] / Vabs + DispT * Vx[i] * Vx[i] / Vabs + DispT * Vy[i] * Vy[i] / Vabs + Dif[i];
			Dispxy[i] = (DispL - DispT) * Vx[i] * Vy[i] / Vabs;
			Dispxz[i] = (DispL - DispT) * Vx[i] * Vz[i] / Vabs;
			Dispyz[i] = (DispL - DispT) * Vy[i] * Vz[i] / Vabs;
		}
	}
	return;
}

void Solute::velocities(const tVector & hNew_, const tVector & Con_, tVector & Vx, tVector & Vy, tVector & Vz) const {
	double Vxx, Vyy, Vzz, Det;
	std::valarray<tIndex> List(4);
	std::valarray<double> B(4), C(4), D(4);

	// see equation 4.25 in the user manual, computation of nodal fluxes
	// note that we threw out Ka which contains anisotropy of Conductivity values and was
	// present in the original code as conaxx,
	// etc. this is probably only correct for a regular spaced grid with equal distances in all directions

	//Vx = 0.; Vy = 0.; Vz = 0.; // Valarrays assumed to be 0 initialized

	// Loop over sub elements
	tIndex count = 0, sub = 0;
	for (list_const_iter it = soluteMesh.getSubelementlist().begin(); it != soluteMesh.getSubelementlist().end(); ++it) { // Loop over sub-elements
		//index values of the nodes of this tetrahedal ellement
		const tIndex i((*it)[0]), j((*it)[1]), k((*it)[2]), l((*it)[3]);

//	        B[0]=-(y[k]-y[j])*(z[l]-z[j])+(y[l]-y[j])*(z[k]-z[j]);
//	        B[1]=+(y[l]-y[k])*(z[i]-z[k])-(y[i]-y[k])*(z[l]-z[k]);
//	        B[2]=-(y[i]-y[l])*(z[j]-z[l])+(y[j]-y[l])*(z[i]-z[l]);
//	        B[3]=+(y[j]-y[i])*(z[k]-z[i])-(y[k]-y[i])*(z[j]-z[i]);
//	        C[0]=+(x[k]-x[j])*(z[l]-z[j])-(x[l]-x[j])*(z[k]-z[j]);
//	        C[1]=-(x[l]-x[k])*(z[i]-z[k])+(x[i]-x[k])*(z[l]-z[k]);
//	        C[2]=+(x[i]-x[l])*(z[j]-z[l])-(x[j]-x[l])*(z[i]-z[l]);
//	        C[3]=-(x[j]-x[i])*(z[k]-z[i])+(x[k]-x[i])*(z[j]-z[i]);
//	        D[0]=-(x[k]-x[j])*(y[l]-y[j])+(x[l]-x[j])*(y[k]-y[j]);
//	        D[1]=+(x[l]-x[k])*(y[i]-y[k])-(x[i]-x[k])*(y[l]-y[k]);
//	        D[2]=-(x[i]-x[l])*(y[j]-y[l])+(x[j]-x[l])*(y[i]-y[l]);
//	        D[3]=+(x[j]-x[i])*(y[k]-y[i])-(x[k]-x[i])*(y[j]-y[i]);
//	        Det=(x[l]-x[i])*B[3]+(y[l]-y[i])*C[3]+(z[l]-z[i])*D[3];

		Det = soluteMesh.getSubelemetVolumes()[sub] * 6.;

		// tension differences over distance per surface area (and so per Ve)
		const tVector & Bi(soluteMesh.getBi()), &Ci(soluteMesh.getCi()), &Di(soluteMesh.getDi());
		Vxx = (Bi[count] * hNew_[i] + Bi[count + 1] * hNew_[j] + Bi[count + 2] * hNew_[k] + Bi[count + 3] * hNew_[l]) / Det;
		Vyy = (Ci[count] * hNew_[i] + Ci[count + 1] * hNew_[j] + Ci[count + 2] * hNew_[k] + Ci[count + 3] * hNew_[l]) / Det;
		Vzz = ((Di[count] * hNew_[i] + Di[count + 1] * hNew_[j] + Di[count + 2] * hNew_[k] + Di[count + 3] * hNew_[l]) / Det) + 1; // Kazz=1, det=6*Ve, see eq 4.25 manual, assuming no further anisotropy

		Vx[i] -= Vxx;
		Vx[j] -= Vxx;
		Vx[k] -= Vxx;
		Vx[l] -= Vxx;

		Vy[i] -= Vyy;
		Vy[j] -= Vyy;
		Vy[k] -= Vyy;
		Vy[l] -= Vyy;

		Vz[i] -= Vzz;
		Vz[j] -= Vzz;
		Vz[k] -= Vzz;
		Vz[l] -= Vzz;

		count += 4;
		++sub;
	}

	// listNE is number of sub elements the node is a corner of
	// con = Nodal values of the hydraulic conductivity at the new time level
	// basically differences in hydraulic head time conductivity averaged over the number of elements per node and taking above some anisotropy in the conductivity into account
	Vx *= Con_ / ListNE;
	Vy *= Con_ / ListNE;
	Vz *= Con_ / ListNE;

	return;
}

/**
 *  No Dirichlet BC can be set in this code: \n
 *  This is because third-type (von Neumann) conditions, in general,
 *	are physically more realistic and preserve solute mass in the simulated system \n
 *	(e.g., van Genuchten and Parker [ 1984]; Leij et al. [1991]).
 */
void Solute::cbound_condition(const double & t, tVector & B) {

	double alf = 1. - 1. / (double) NLEVEL; // = 0.5, BECAUSE CRANK-NICOLSON

	tVector Qwater;
	if (fine) {
		Qwater.resize(NumNP);
		Mapping::interpolate(water_.getMesh(), water_.getQ(), Qwater, factor);
		Qwater /= pow(factor, 2);
	}
	const tVector & Q = fine ? Qwater : water_.getQ(); //nodal water fluxes in ml/day

	//topboundary
	double cumTopBoundaryFlux(0);

	//reset

	// concentration in precipitation
	double prec;
	if (sPrec_) {
		sPrec_->get(t, prec);
		prec = std::fabs(prec);
		if (sinterception_) {
			double ic;
			sinterception_->get(t, ic);
			prec -= std::fabs(ic);
		}

	} else {
		prec = 0.;
	}
	double cPrec;
	if (scPrec_) {
		scPrec_->get(t, cPrec);
		cPrec = std::fabs(cPrec) * prec;
	} else {
		cPrec = 0.;
	}

	auto & width(soluteMesh.getWidth());

	//add precipitation
	if (cPrec) {
		msg::warning("Solute: fertigating, check this code first to see if signs are right, balance is correct etc");
		for (auto i : soluteMesh.getTopBoundaryNodes()) {
			B[i] += width[i] * cPrec;
			cumTopBoundaryFlux += -width[i] * cPrec;
		}
	}
	for (auto i : soluteMesh.getTopBoundaryNodes()) {
		// evaporation does not cause a solute flux, so correct
		B[i] += Q[i] * (alf * Conc[i]);
		soluteMatrix.addValueUnsafely(i, i, -1. / (double) NLEVEL * Q[i]);
	}
	if (sumTopBoundaryFlux_)
		sumTopBoundaryFlux_->set(t, cumTopBoundaryFlux);

	double cumBottomBoundaryFlux(0);
	for (auto i : soluteMesh.getBottomBoundaryNodes()) {
		//if (Q[i] > 0.) { // if fluid is flowing in
		//	Cauchy boundary condition both neuman and dirichlet are in effect
		B[i] += -Q[i] * (Conc[i] - alf * Conc[i]);
		soluteMatrix.addValueUnsafely(i, i, -1. / (double) NLEVEL * Q[i]);
		cumBottomBoundaryFlux += Q[i] * Conc[i];
	}
	if (sumBottomBoundaryFlux_)
		sumBottomBoundaryFlux_->set(t, cumBottomBoundaryFlux);

	return;
}

void Solute::interpolate(const tVector & source, tVector & target, const int factor) const {

	if (2 != factor) {
		msg::error("interpolate: factor != 2");
	}

	const tIndex & nxw(water_.getMesh().getnx());
	const tIndex & nyw(water_.getMesh().getny());

	tIndex NX = (2 * nxw - 1);
	tIndex NY = (2 * nyw - 1);
	tIndex numnp_2x = NX * NY * (water_.getMesh().getnz() * 2 - 1);
	// tIndex numnp_4x = (2*NX-1)*(2*NY-1)*(2*(water_.getMesh().getnz()*2-1)-1);

	// exception if source and target dimensions do not fit together
	if (NumNP != numnp_2x /*|| NumNP != numnp_4x */) {
		msg::error("interpolate: Fine Mesh is not a two potency multiple of coarse mesh");
	}

	tIndex i = 0;
	tIndex counter = 0;
	tIndex row = 1;
	tIndex layer = 1;
	while (counter < numnp_2x /*&& i < sourceNumNP*/) {
		if (i < row * nxw - 1) {
			target[counter] = source[i];
			target[counter + 1] = 0.5 * (source[i] + source[i + 1]);
			counter += 2;
		} else {
			target[counter] = source[i];
			counter += NX + 1;
			++row;
			if (i == layer * nxw * nyw - 1) {
				counter += (NX * (NY - 1));
				++layer;
			}
		}
		++i;
	}

	i = 0;
	counter = NX;
	row = 1;
	tIndex counter2 = 0;
	while (counter < numnp_2x) {
		if (i < row * nxw - 1) {
			target[counter] = 0.5 * (source[i] + source[i + nxw]);
			target[counter + 1] = 0.25 * (source[i] + source[i + 1] + source[i + nxw] + source[i + nxw + 1]);
			counter += 2;
		} else {
			target[counter] = 0.5 * (source[i] + source[i + nxw]);
			counter += NX + 1;
			++row;
			++counter2;
			if (counter2 == nyw - 1) {
				counter += (NX * NY + NX);
				i += nxw;
				++row; // then i is at the end of the layer
				counter2 = 0;
			}
		}
		++i;
	}

	i = 0;
	counter = NX * NY;
	row = 1;
	layer = 1;
	while (counter < numnp_2x) {
		if (i < row * nxw - 1) {
			target[counter] = 0.5 * (source[i] + source[i + nxw * nyw]);
			target[counter + 1] = 0.25 * (source[i] + source[i + 1] + source[i + nxw * nyw] + source[i + 1 + nxw * nyw]);
			counter += 2;
		} else {
			target[counter] = 0.5 * (source[i] + source[i + nxw * nyw]);
			counter += NX + 1;
			++row;
			if (i == layer * nxw * nyw - 1) {
				counter += (NX * (NY - 1));
				++layer;
			}
		}
		++i;
	}

	i = 0;
	counter = NX * NY + NX;
	row = 1;
	counter2 = 0;
	while (counter < numnp_2x) {
		if (i < row * nxw - 1) {
			target[counter] = 0.25 * (source[i + nxw] + source[i] + source[i + nyw * nxw] + source[i + nyw * nxw + nxw]);
			target[counter + 1] = 0.125
					* (source[i + nxw] + source[i] + source[i + nyw * nxw] + source[i + nyw * nxw + nxw] + source[i + nxw + 1] + source[i + 1]
							+ source[i + nyw * nxw + 1] + source[i + nyw * nxw + nxw + 1]);
			counter += 2;
		} else {
			target[counter] = 0.25 * (source[i + nxw] + source[i] + source[i + nyw * nxw] + source[i + nyw * nxw + nxw]);
			counter += NX + 1;
			++row;
			++counter2;
			if (counter2 == nyw - 1) {
				counter += (NX * NY + NX);
				i += nxw;
				++row; // then i is at the end of the layer
				counter2 = 0;
			}
		}
		++i;
	}

	return;

}

// todo right now this method thinks that meshes are the same, here only soluteMesh is used because it has water values anyway.
void Solute::reportTotSolidColumn(const Time & t) {

	double ConVolSolid(0.); 	/// Amount of solute in the entire flow domain [M](ConVol in Table 9.6).
	double ConVolSolute(0.); 	/// Amount of solute in the entire flow domain [M](ConVol in Table 9.6).

	tVector thNew_water;
	if (fine) {
		thNew_water.resize(NumNP);
		Mapping::interpolate(water_.getMesh(), water_.getthNew(), thNew_water, factor);
	}
	const tVector &thNew = fine ? thNew_water : water_.getthNew();

	tIndex it_subelem = 0;
	for (auto& it : soluteMesh.getSubelementlist()) { // Loop over sub-elements
		//index values of the nodes of this tetrahedral element
		const tIndex i(it[0]), j(it[1]), k(it[2]), l(it[3]);
		const double VE = soluteMesh.getSubelemetVolumes()[it_subelem];

		ConVolSolute += VE / 4. * (thNew[i] * Conc[i] + thNew[j] * Conc[j] + thNew[k] * Conc[k] + thNew[l] * Conc[l]);
		ConVolSolid += VE / 4.
				* (bulkDensity[i] * adsorptionCoeff[i] * Conc[i] + bulkDensity[j] * adsorptionCoeff[j] * Conc[j] + bulkDensity[k] * adsorptionCoeff[k] * Conc[k]
						+ bulkDensity[l] * adsorptionCoeff[l] * Conc[l]);

		++it_subelem;
	}
	double ConVol = ConVolSolid + ConVolSolute; /// Amount of solute in the entire flow domain [M](ConVol in Table 9.6).

	//jouke: write metrics to simroot tables
	if (t < 1e-10) {
		ConVolI = ConVol; // set this here once to keep initial results
	}
	if (totSoluteInColumn)
		totSoluteInColumn->set(t, ConVol);
	if (totSoluteInColumnDissolved)
		totSoluteInColumnDissolved->set(t, ConVolSolute);
	if (totSoluteInColumnAbsorbed)
		totSoluteInColumnAbsorbed->set(t, ConVolSolid);
	if (totSoluteChange)
		totSoluteChange->set(t, ConVol - ConVolI);

	//d95 for solute in column
	//todo with a better interface, this could go into a separate class.
	if (d95Solute) {
		double tresh = 0.10 * ConVol;
		double tot(0.), d90;
		for (auto it = soluteMesh.getSubelementlist().rbegin(); it != soluteMesh.getSubelementlist().rend(); ++it) { // Loop over sub-elements
			//index values of the nodes of this tetrahedal ellement
			it_subelem = 0;
			const tIndex i((*it)[0]), j((*it)[1]), k((*it)[2]), l((*it)[3]);
			const double VE = soluteMesh.getSubelemetVolumes()[it_subelem];
			// weights are 1/4th!
			tot += VE / 4. * (thNew[i] * Conc[i] + thNew[j] * Conc[j] + thNew[k] * Conc[k] + thNew[l] * Conc[l]);
			tot += VE / 4.
					* (bulkDensity[i] * adsorptionCoeff[i] * Conc[i] + bulkDensity[j] * adsorptionCoeff[j] * Conc[j]
							+ bulkDensity[k] * adsorptionCoeff[k] * Conc[k] + bulkDensity[l] * adsorptionCoeff[l] * Conc[l]);

			if (tot > tresh) {
				d90 = std::min(std::min(z[i], z[j]), std::min(z[k], z[l]));
				break;
			}
			++it_subelem;
		}
		d95Solute->set(t, d90); // set the value in the table
	}

	return;
}

void Solute::setSoluteSink(const double & t0, const double & t1, tVector & sink) {
	if (t1 <= t0)
		return;
	sink = 0;
	double cSinkCheck(0);
	for (unsigned int i(0); i < femIndexes_.size(); ++i) {
		if (!nutrientUptake_[i]->evaluateTime(t0))
			continue;
		const std::vector<int> &ind = femIndexes_[i];
		if (ind[0] == -1) {
			// ROOTS OUTSIDE GRID no nutrients and no water
			continue;
		} else {
			double nuptake, nu1;
			nutrientUptake_[i]->get(t1, nuptake); // calls MichaelisMenten::MichaelisMenten(SimulaDynamic* pSD) the first time
			nutrientUptake_[i]->get(t0, nu1);
			nuptake -= nu1; //uptake by true differencing
			nuptake /= (t1 - t0); //uptake per day
			if (nuptake < 1e-30)
				continue;
			cSinkCheck += nuptake;

			std::vector<int> const * pind;
			std::vector<int> indexes;
			std::vector<double> weights;
			std::vector<double> const * pweights;
			double sum;
			//double * psum;
			if (ind[0] == -2) {
				//growthpoint, needs matching
				Coordinate pos;
				nutrientUptake_[i]->getParent(2)->get(t0, pos);
				const std::vector<double> &base = femWeights_[i];
				pos.x += base[0];
				pos.y += base[1];
				pos.z += base[2];
				soluteMesh.matchNodes(pos, indexes, weights, sum); //todo: somewhat ugly, depends on order ofthings in several other places
				pind = &indexes;
				pweights = &weights;
				//psum=&sum;
			} else {
				pind = &ind;
				pweights = &femWeights_[i];
				//psum=&femSumWeights_[i];
			}
			if (pind->size() > 1) {
				double sumweights(0);
				for (unsigned int j = 0; j < pind->size(); ++j) {
					int femi = pind->at(j);
					//if(Conc[femi]>0){
					sumweights += Conc[femi] * (*pweights)[j];
					//}
				}
				//std::cout<<std::endl<< nutrientUptake_[i]->getParent(3)->getName()<<" "<<nutrientUptake_[i]->getParent(5)->getName();
				//std::cout.precision(5);
				for (unsigned int j = 0; j < pind->size(); ++j) {
					int femi = pind->at(j);
					//if(Conc[femi]>0){
					sink[femi] += nuptake * Conc[femi] * (*pweights)[j] / sumweights;
					//}
					//std::cout<<std::endl<<nuptake<<" "<<femi<<" "<<conc[femi]<<" "<<Conc[femi]*(*pweights)[j]/sumweights;
				}
			} else {
				//push the first nodal values
				int femi = pind->at(0);
				sink[femi] += nuptake;
			}

		}

	}
	if (std::abs(cSinkCheck) > 1.) {
		double check = 100 * std::abs((cSinkCheck - sink.sum()) / cSinkCheck);
		if (check > 0.01) {
			msg::error("Solute: sink differs at time=" + std::to_string(t0) + " from uptake by " + std::to_string(check) + " percent.");
		}
	}

	if (totSink_)
		totSink_->set(t0, -1*(sink.sum()) );

	sink /= soluteMesh.getvol();

	//----------------mineralization--------------------------------
	if (mineral) {
		const double reduct = 1.; // this is a reduction coefficient for when net immobilization is limited by N - not used now
		tVector rndiss(mineral->getNumNodesWithMineralization());

		mineral->mineralisationYang1994(t1 - t0, reduct, water_.gethOld(), rndiss);
		// weights not used
		//todo check if this is true, negative values should be OK
		if (rndiss.max() > 0.) {
			msg::error("swms: yang`s model calculated net immobilization instead of mineralisation, but swms can not handle that.");
		}

		//mineralization added to the sink term.
		//only top numnl nodes have mineralization.
		for (unsigned int i(0); i < rndiss.size(); ++i) {
			sink[i] += rndiss[i] / 14.; // cSink += rndiss / 14.; // rndiss is in ug N/cm3 convert to uMol/cm3
		}

		//reporting of mineralizationrate
		double totminrate = (rndiss*soluteMesh.getvol()).sum() / 14.;
		if (totalMineralization_)
			totalMineralization_->set(t0, -1*totminrate);
		if (totminrate > 1e-5)
			msg::warning("Mineralization effective, but code wasn't thoroughly tested yet. Check balance etc.");
	}

}

void Solute::updateLists(const Time &t, const std::deque<SimulaBase *> & rootNodeQueue_) {

	for (auto & it : rootNodeQueue_) {
		SimulaBase* p = it->getChild(nutrient);
		SimulaBase* u = p->getChild("rootSegmentNutrientUptake");
		nutrientUptake_.push_back(u);

		u = p->existingChild("nutrientConcentrationAtTheRootSurface");
		SimulaTable<Time> *tbl;
		if (u) {
			tbl = dynamic_cast<SimulaTable<Time> *>(u);
			if (!tbl)
				msg::error("Solute::updateLists: concentrationAtRootSurface must be a time table");
			nutrientConcentration_.push_back(tbl);
		} else {
			nutrientConcentration_.push_back(nullptr);
			msg::warning("Solute::updateLists: nutrientConcentrationAtTheRootSurface not found");
		}

		//u=p->getChild("Cmin");

		if (it->getName() == "growthpoint") {
			Coordinate pos;
			it->getBase(t, pos);
			std::vector<double> w(3);
			w[0] = pos.x;
			w[1] = pos.y;
			w[2] = pos.z;
			femIndexes_.push_back(std::vector<int>(1, -2));
			femWeights_.push_back(w);
			//femSumWeights_.push_back(-1);
		} else {
			std::vector<int> indexes(27);
			std::vector<double> weights(27);
			double sum; // TODO unused?
			Coordinate pos;
			it->getAbsolute(t, pos);
			soluteMesh.matchNodes(pos, indexes, weights, sum);
			femIndexes_.push_back(indexes);
			femWeights_.push_back(weights);
			//femSumWeights_.push_back(sum);
		}

	}

}

void Solute::setrootsurfacevalues(const Time & t) {
	//update list

	double concsur;
	for (unsigned int i(0); i < femIndexes_.size(); ++i) {
		if (!nutrientUptake_[i]->evaluateTime(t))
			continue;
		std::vector<int> &ind = femIndexes_[i];
		if (ind[0] == -1) {
			//!ROOTS OUTSIDE GRID no nutrients
			concsur = 0;
		} else {
			std::vector<int> * pind;
			std::vector<int> indexes;
			std::vector<double> weights;
			std::vector<double> * pweights;
			double sum;
			//double * psum;
			if (ind[0] == -2) {

				//growthpoint, needs matching
				Coordinate pos;
				nutrientUptake_[i]->getParent(2)->get(t, pos);
				std::vector<double> &base = femWeights_[i];
				pos.x += base[0];
				pos.y += base[1];
				pos.z += base[2];
				soluteMesh.matchNodes(pos, indexes, weights, sum);
				pind = &indexes;
				pweights = &weights;
				//psum=&sum;
			} else {
				pind = &ind;
				pweights = &femWeights_[i];
				//psum=&Simunek::femSumWeights_[i];
			}
			if (pind->size() > 1) {
				//determine weights based on distances
				concsur = 0;
				double tmpsc = 0;
				for (unsigned int j = 0; j < pind->size(); ++j) {
					int femi = pind->at(j);
					double w(pweights->at(j));
					//if(Conc[femi]>0){//ignore nodes with negative concentrations
					tmpsc += w;
					concsur += w * Conc[femi];
					//}
					//std::cout<<std::endl<<pweights->at(j)<<" * "<<Conc[femi]<<" x "<<x[femi]<<" y "<<y[femi]<<" z "<<z[femi];
				}
				//std::cout<<std::endl<<"sum:"<<concsur<<" / "<<*psum<<" / "<<tmps;
				if (tmpsc > 0)
					concsur /= tmpsc;
			} else {
				//push the first nodal values
				int femi = pind->at(0);
				concsur = Conc[femi];
			}
			//sometimes at high sinks (due too high rld) model migth have temporal small negative numbers for concentration
			concsur=std::max(concsur, 0.);

		}
//		std::cout.precision(5);
		if (nutrientConcentration_[i])
			nutrientConcentration_[i]->set(t, concsur);
	}

}

SoluteMassBalanceTest::SoluteMassBalanceTest(SimulaDynamic* pSD) :
		DerivativeBase(pSD), seedReserves(0) {
	std::string nutrient(pSD->getParent()->getName());
	SimulaBase *p = ORIGIN->getChild("soil")->getChild(nutrient);
	SimulaBase *d;

	//todo come up with better names for these
	totSoluteChange = p->getChild("totalSoluteChange", "umol");
	sumTopBoundaryFlux_ = p->getChild("topBoundaryFlux", "umol");
	sumBottomBoundaryFlux_ = p->getChild("bottomBoundaryFlux", "umol");
	totalMineralization_ = p->existingChild("totalMineralization", "umol");
	totSink_ = p->getChild("totalSink", "umol");
	totalNutrientUptake_ = ORIGIN->existingPath("/plants/" + pSD->getParent()->getName() + "/plantNutrientUptake", "umol");
	//todo for backwards compatibility, simply insert this, instead of not calculating the balance
	if (totalNutrientUptake_)
		totalNutrientUptake_->get(0, seedReserves);					//todo these are a problem, when plants are not planted at time 0

	d = p->existingChild("massBalanceError");
	if (d)
		massBalanceError_ = dynamic_cast<SimulaTable<Time>*>(d);

	 p->getChild("totalSoluteInColumn", "umol")->get(0.,ref);
}

void SoluteMassBalanceTest::calculate(const Time & t, double &error) {
	double change(0);
	totSoluteChange->get(t, change);

	double tboundary(0);
	sumTopBoundaryFlux_->get(t, tboundary);

	double bboundary(0);
	sumBottomBoundaryFlux_->get(t, bboundary);

	double mineralization(0);
	if (totalMineralization_)
		totalMineralization_->get(t, mineralization);

	double sink(0);
	totSink_->get(t, sink);

	double uptake(0);
	if (totalNutrientUptake_) {
		totalNutrientUptake_->get(t, uptake);
//note that sink is negative when roots take up but uptake is positive
		if (std::fabs(sink) > 10 && std::fabs((sink + uptake - seedReserves)/sink) > 0.05)
			//std::cout<<std::endl<< uptake << " " << sink <<" "<<seedReserves;
			msg::warning(
					"SoluteMassBalanceTest: uptake may not equal sink for " + pSD->getParent()->getName()
							+ ". Message could be caused by seed reserves added to simulation after t=0, or roots outside the grid etc, but please check");
	} else {
		msg::warning("SoluteMassBalanceTest: not including nutrient uptake as the table is missing");
	}

	//todo sink contains mineralization
	double sumFluxes(tboundary + bboundary + mineralization + sink);
	error = (change - sumFluxes) / ref;
	if (std::fabs(error) > 0.01)
		msg::warning("SoluteMassBalanceTest: mass balance is of by more than 1% of initial "+std::to_string(ref)+" umol in column.");

	if (massBalanceError_)
		massBalanceError_->set(t, change - sumFluxes);
}

std::string SoluteMassBalanceTest::getName() const {
	return "soluteMassBalanceTest";
}

DerivativeBase * newInstantiationSoluteMassBalanceTest(SimulaDynamic* const pSD) {
	return new SoluteMassBalanceTest(pSD);
}

//Register the module
class AutoRegisterSoluteFunctions {
public:
	AutoRegisterSoluteFunctions() {
		BaseClassesMap::getDerivativeBaseClasses()["soluteMassBalanceTest"] = newInstantiationSoluteMassBalanceTest;
	}
};

// our one instance of the proxy
static AutoRegisterSoluteFunctions l8r9h38hr9h9h9;

