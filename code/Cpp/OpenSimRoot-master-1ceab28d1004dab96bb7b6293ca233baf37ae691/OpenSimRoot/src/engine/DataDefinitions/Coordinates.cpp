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

#include "Coordinates.hpp"
#include "../../cli/Messages.hpp"
#include <vector>
#include <math.h>
#define square(X) pow(X,2)

//==================coordinates==================
//constructors
Coordinate::Coordinate():x(0),y(0),z(0),estimate(false){}
Coordinate::Coordinate(const Coordinate& p):
	x(p.x),y(p.y),z(p.z),estimate(p.estimate)
{}
Coordinate::Coordinate(const MovingCoordinate& p):
	x(p.x),y(p.y),z(p.z),estimate(p.estimate)
{}
Coordinate::Coordinate(const double &newX, const double &newY, const double &newZ, const bool &newEst):
		x(newX),y(newY),z(newZ),estimate(newEst)
{}

//overloading aritmatic operators
Coordinate & Coordinate::operator= (const Coordinate &p) {
		x=p.x;
		y=p.y;
		z=p.z;
		estimate=p.estimate;
		return *this;
}
void Coordinate::operator+= (const Coordinate &p)  {
		x += p.x;
		y += p.y;
		z += p.z;
		estimate+=p.estimate;
}
void Coordinate::operator-= (const Coordinate &p)  {
		x -= p.x;
		y -= p.y;
		z -= p.z;
		estimate+=p.estimate;
}
void Coordinate::operator*= (const State &d) {
		x*=d.value;
		y*=d.value;
		z*=d.value;
		estimate+=d.estimate;
}
void Coordinate::operator/= (const State &d) {
		x/=d.value;
		y/=d.value;
		z/=d.value;
		estimate+=d.estimate;
}
Coordinate Coordinate::operator+ (const Coordinate &p) const {
		Coordinate r(*this);
		r.x += p.x;
		r.y += p.y;
		r.z += p.z;
		r.estimate+=p.estimate;
		return r;
}
Coordinate Coordinate::operator- (const Coordinate &p) const {
		//name of right side is maintained
		Coordinate r(*this);
		r.x -= p.x;
		r.y -= p.y;
		r.z -= p.z;
		r.estimate+=p.estimate;
		return r;
	}
Coordinate Coordinate::operator+ (const State &p) const {
		Coordinate r(*this);
		r.x += p.value;
		r.y += p.value;
		r.z += p.value;
		r.estimate+=p.estimate;
		return r;
}
Coordinate Coordinate::operator- (const State &p) const {
		//name of right side is maintained
		Coordinate r(*this);
		r.x -= p.value;
		r.y -= p.value;
		r.z -= p.value;
		r.estimate+=p.estimate;
		return r;
	}
Coordinate Coordinate::operator* (const State &d)const  {
		Coordinate r(*this);
		r.x *= d.value  ;
		r.y *= d.value ;
		r.z *= d.value ;
		r.estimate+=d.estimate+estimate;
		return r;
	}
Coordinate Coordinate::operator* (const Coordinate &p)const{
	Coordinate r(*this);
	r.x*=p.x;
	r.y*=p.y;
	r.z*=p.z;
	r.estimate+=p.estimate;
	return r;
}
Coordinate Coordinate::operator/ (const Coordinate &p)const{
	Coordinate r(*this);
	r.x/=p.x;
	r.y/=p.y;
	r.z/=p.z;
	r.estimate+=p.estimate;
	return r;
}

State Coordinate::operator^ (const Coordinate &p)const{
	State r(square(x-p.x) + square(y-p.y) + square(z-p.z),estimate+p.estimate);
	return r;
}


Coordinate Coordinate::operator/ (const State &d) const {
		Coordinate r(*this);
		if(d==0) {
			msg::error("Coordinate::operator/ : can't divide through zero");
		}else{
			r.x /=d.value;
			r.y /=d.value;
			r.z /=d.value;
		};
		r.estimate+=d.estimate;
		return r;
	}

bool  Coordinate::operator== (const Coordinate &d)const {
	if(x!=d.x) return false;
	if(y!=d.y) return false;
	if(z!=d.z) return false;
	return true;
}
bool  Coordinate::operator!= (const Coordinate &d)const {
	if(x!=d.x) return true;
	if(y!=d.y) return true;
	if(z!=d.z) return true;
	return false;
}

bool  Coordinate::operator> (const Coordinate &d)const{
	if(y>d.y){
		return true;
	}else if(y<d.y){
		return false;
	}else if(x>d.x){
		return true;
	}else if(x<d.x){
		return false;
	}else if(z>d.z){
		return true;
	}else{
		return false;
	}
}
bool  Coordinate::operator< (const Coordinate &d)const{
	if(y<d.y){
		return true;
	}else if(y>d.y){
		return false;
	}else if(x<d.x){
		return true;
	}else if(x>d.x){
		return false;
	}else if(z<d.z){
		return true;
	}else{
		return false;
	}
}
bool  Coordinate::operator>= (const Coordinate &d)const{
	if(*this<d){
		return false;
	}else{
		return true;
	}
}
bool  Coordinate::operator<= (const Coordinate &d)const{
	if(*this>d){
		return false;
	}else{
		return true;
	}
}


void Coordinate::split(State & xs, State & ys, State & zs)const{
	xs.value=x;
	xs.estimate=estimate;
	ys.value=y;
	ys.estimate=estimate;
	zs.value=z;
	zs.estimate=estimate;
}
void Coordinate::join(const State & xs, const State & ys, const State & zs){
	x=xs.value;
	y=ys.value;
	z=zs.value;
	estimate=xs.estimate+ys.estimate+zs.estimate;
}
void Coordinate::normalize(){
	State l=this->length();
	if(l!=0) *this/=l;
}
double Coordinate::length()const{
	double r(sqrt( square(x) + square(y) + square(z)));
	return r;
}
void Coordinate::setLength(const State & l){
	this->normalize();
	*this*=l;
}
double Coordinate::sum()const{
	double r(x + y + z);
	return r;
}
double Coordinate::squareSum()const{
	double r(square(x) + square(y) + square(z));
	return r;
}


std::ostream &operator<<(std::ostream &os, const Coordinate &obj){
	os	<<obj.x<<'\t'<<obj.y<<'\t'<<obj.z;
	return os;
}



MovingCoordinate::MovingCoordinate():
	Coordinate(),state(*this),rate(0,0,0,true){}
MovingCoordinate::MovingCoordinate(const MovingCoordinate& p):
	Coordinate(p.state),state(*this),rate(p.rate){}
MovingCoordinate::MovingCoordinate(const Coordinate& p):
	Coordinate(p),state(*this),rate(0,0,0,true){}
MovingCoordinate::MovingCoordinate(const double &newX, const double &newY, const double &newZ):
	Coordinate(newX,newY,newZ),state(*this),rate(0,0,0,true){}


bool  MovingCoordinate::operator== (const MovingCoordinate &d)const{
	if(state!=d.state) return false;
	if(rate!=d.rate) return false;
	return true;
}
bool  MovingCoordinate::operator!= (const MovingCoordinate &d)const{
	if(state==d.state) return false;
	if(rate==d.rate) return false;
	return true;
}

	//overloading aritmatic operators
MovingCoordinate & MovingCoordinate::operator= (const MovingCoordinate &p) {
	state=p.state;
	rate=p.rate;
//	estimate=p.estimate;
//	rate.estimate=p.rate.estimate;
	return *this;
}
MovingCoordinate MovingCoordinate::operator+ (const MovingCoordinate &p) const {
	MovingCoordinate r(*this);
	r.state+=p.state;
	r.rate+=p.rate;
//	r.estimate+=p.estimate;
//	rate.estimate+=p.rate.estimate;
	return r;
}
MovingCoordinate MovingCoordinate::operator- (const MovingCoordinate &p) const {
	MovingCoordinate r(*this);
	r.state-=p.state;
	r.rate-=p.rate;
//	r.estimate+=p.estimate;
//	rate.estimate+=p.rate.estimate;
	return r;
}
void MovingCoordinate::operator+= (const MovingCoordinate &p){
	state+=p.state;
	rate+=p.rate;
//	estimate+=p.estimate;
//	rate.estimate+=p.rate.estimate;
}
void MovingCoordinate::operator-= (const MovingCoordinate &p){
	state-=p.state;
	rate-=p.rate;
//	estimate+=p.estimate;
//	rate.estimate+=p.rate.estimate;
}
MovingCoordinate MovingCoordinate::operator* (const State &d)const{
	MovingCoordinate r(*this);
	r.state*=d;
	r.rate*=d;
	return r;
}
Coordinate MovingCoordinate::operator/ (const State &d) const{
	MovingCoordinate r(*this);
	r.state/=d;
	r.rate/=d;
	return r;
}

void MovingCoordinate::split(StateRate & xs, StateRate & ys, StateRate & zs)const{
	xs.state=x;
	xs.estimate=estimate;
	xs.rate=rate.x;
	xs.rate.estimate=estimate;
	ys.state=y;
	ys.estimate=estimate;
	ys.rate=rate.y;
	ys.rate.estimate=estimate;
	zs.state=z;
	zs.estimate=estimate;
	zs.rate=rate.z;
	zs.rate.estimate=estimate;
}
void MovingCoordinate::join(const StateRate & xs, const StateRate & ys, const StateRate & zs){
	x=xs.value;
	y=ys.value;
	z=zs.value;
	rate.x=xs.rate.value;
	rate.y=ys.rate.value;
	rate.z=zs.rate.value;
	estimate=xs.estimate+ys.estimate+zs.estimate;
	rate.estimate=xs.rate.estimate+ys.rate.estimate+zs.rate.estimate;
}
State MovingCoordinate::length()const{
	State r(sqrt( square(x) + square(y) + square(z)));
	r.estimate=this->estimate;
	return r;
}


std::ostream &operator<<(std::ostream &os, const MovingCoordinate &obj){
	os	<<obj.state<<"\t"<<obj.rate<<"\t";
	return os;
}





