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

/* See headerfile for explanation and change log*/
/*
 * The library contains mathematical functions that are used in SimRoot such as interpolation and integration.
 *
 */


#include "MathLibrary.hpp"

#include "../cli/Messages.hpp"

#if _WIN32 || _WIN64
#include <algorithm>
#endif

int quadraticFormula(const double & a,const double & b,const double & c,double & sol1,double & sol2){
	double t1(2*a);
	double t2(-b/t1);
	double t4(square(b)-(4*a*c));
	if(a==0) {
		//this is not a quadratic formula but a linear one
		if(b==0) {
			//msg::error("quadraticFormula: no solution exists");
			return(2);
		}
		sol1= -c/b;
		sol2=sol1;

	}
	if(t4<0){
		//msg::error("quadraticFormula: no solution exists");
		return(1);
	}
	double t3(sqrt(t4)/t1);
	sol1=t2+t3;
	sol2=t2-t3;
	return(0);
}


/*tridiagonal matrix solver of the form A*x = d. Where A is
 *
 * bl c 0 0 0 0 0 0
 * a b c 0 0 0 0 0
 * 0 a b c 0 0 0 0
 * 0 0 a b c 0 0 0
 * 0 0 0 a b c 0 0
 * 0 0 0 0 a b c 0
 * 0 0 0 0 0 a b c
 * 0 0 0 0 0 0 a br
 *
 * and dim A, length x and length b are all equal
 */
// ///@todo note that this could be optimized, e could be saved between timesteps since it does not depend on q or d
void tridiagonalSolver(const double & bl , const double &a, const double &b,const double &c, const double &br, const std::vector<double> &d, std::vector<double> & x){
	//dimensions of the matrix (-1 because in c first element is 0)
	const std::size_t n(d.size()-1);
	if (n<3) msg::error("tridiagonalSolver: matrix should be at least 3x3");
	//declaration of two vectors
	std::vector<double> f(d.size(),0),e(d.size(),0);
	//Gauss elimination on the matrix, saving information needed for the backsolve.
	//initial conditions
	e[0]=-a/bl; f[0]=d[0]/bl;
	for(std::size_t i(1) ; i<n ; i++){
		double den = b+c*e[i-1];
		e[i] = -a/den;
		f[i] = (d[i]-c*f[i-1])/den;
	}
	//steps of "back substitution".
	x[n] = (d[n]-c*f[n-1])/(br+c*e[n-1]);
	for(int i((int)n-1) ; i>=0 ; i--)	x[i] = e[i] * x[i+1] + f[i];
}

//another tridiagonal matix solver for which a and c are vectors
void tridiagonalSolver(std::vector<double> & adiag, const std::vector<double> & aleft, std::vector<double> & arite, std::vector<double> &x){
	//dimensions of the matrix (-1 because in c first element is 0)
	const std::size_t n(x.size()-1);
	if (n<2) msg::error("tridiagonalSolver: matrix should be at least 3x3");
	if(adiag.size() != aleft.size() ||  adiag.size() != arite.size() || adiag.size() != x.size()) msg::error("tridiagonalSolver: vectors of unequal length");
	//Gauss elimination on the matrix, saving information needed for the backsolve.
	arite[0] = arite[0] / adiag[0];
	for (std::size_t i(1); i < n; i++ ){
		adiag[i] = adiag[i] - aleft[i] * arite[i-1];
		arite[i] = arite[i] / adiag[i];
	}
	adiag[n] = adiag[n] - aleft[n] * arite[n-1];
	// Carry out the same elimination steps on F that were done to the matrix.
	x[0] = x[0] / adiag[0];
	for (std::size_t i(1); i <= n ; i++ ) x[i] = ( x[i] - aleft[i] * x[i-1] ) / adiag[i];
	// And now carry out the steps of "back substitution".
	for(int i((int)n-1) ; i >= 0 ; i--)	x[i] = x[i] - arite[i] * x[i+1];
}

void tridiagonalSolver(std::valarray<double> & adiag, const std::valarray<double> & aleft, std::valarray<double> & arite, std::valarray<double> &x){
	//dimensions of the matrix (-1 because in c first element is 0)
	const std::size_t n(x.size()-1);
	//if (n<2) msg::error("tridiagonalSolver: matrix should be at least 3x3");
	//if(adiag.size() != aleft.size() ||  adiag.size() != arite.size() || adiag.size() != x.size()) msg::error("tridiagonalSolver: vectors of unequal length");
	//Gauss elimination on the matrix, saving information needed for the backsolve.
	arite[0] = arite[0] / adiag[0];
	for (std::size_t i(1); i < n; i++ ){
		adiag[i] = adiag[i] - aleft[i] * arite[i-1];
		arite[i] = arite[i] / adiag[i];
	}
	adiag[n] = adiag[n] - aleft[n] * arite[n-1];
	// Carry out the same elimination steps on F that were done to the matrix.
	x[0] = x[0] / adiag[0];
	for (std::size_t i(1); i <= n ; i++ ) x[i] = ( x[i] - aleft[i] * x[i-1] ) / adiag[i];
	// And now carry out the steps of "back substitution".
	for(int i((int)n-1) ; i >= 0 ; i--)	x[i] = x[i] - arite[i] * x[i+1];
}

/*tridiagonal matrix solver of the form A*y = f. Where A is
 *
 * a c 0 0 0 0 0 0
 * b a c 0 0 0 0 0
 * 0 b a c 0 0 0 0
 * 0 0 b a c 0 0 0
 * 0 0 0 b a c 0 0
 * 0 0 0 0 b a c 0
 * 0 0 0 0 0 b a c
 * 0 0 0 0 0 0 b a
 *
 * and dim A, length x and length b are all equal
 */

extern"C" {
	void isnotanumber(double *number,bool *check){
		*check = std::isnan(*number);
	}
}
double elementMultiplication(std::vector<double> &d){
	double r(1);
	for(std::vector<double>::iterator it(d.begin()); it!=d.end() ;++it){
		r*=(*it);
	}
	return r;
}
double vectorsMinimum(std::vector<double> &d){
	double r(1e10);//todo replace with macro containing max for doubles
	for(std::vector<double>::iterator it(d.begin()); it!=d.end() ;++it){
		r=std::min(r,(*it));
	}
	return r;
}
double vectorsMaximum(std::vector<double> &d){
	double r(-1e10);
	for(std::vector<double>::iterator it(d.begin()); it!=d.end() ;++it){
		r=std::max(r,(*it));
	}
	return r;
}
double vectorsAverage(std::vector<double> &d){
	double r(0);
	for(std::vector<double>::iterator it(d.begin()); it!=d.end() ;++it){
		r+=(*it);
	}
	return (r/(double)d.size());
}
double vectorsMaxRelativeDeviationFromOne(std::vector<double> &d){
	double r(0),absr(-1);
	for(std::vector<double>::iterator it(d.begin()); it!=d.end() ;++it){
		double delta=(*it) - 1.0;
		if(absr<fabs(delta)){
			absr=fabs(delta);
			r=delta;
		}
	}
	return (r+1.0);
}


