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

#ifndef COORDINATES_HPP_
#define COORDINATES_HPP_

#include <iostream>
#include "StateRate.hpp"

//Structure for storing coordinates
struct MovingCoordinate;

struct Coordinate{
	//data
	double x, y, z;
	bool estimate;

	//constructors
	Coordinate();
	Coordinate(const Coordinate& p);
	Coordinate(const MovingCoordinate& p);
	Coordinate(const double &newX, const double &newY, const double &newZ, const bool &newEst=false);

	//TODO overloading =

	//overloading aritmatic operators
	Coordinate & operator= (const Coordinate &p) ;
	void operator+= (const Coordinate &p);
	void operator-= (const Coordinate &p);
	void operator*= (const State &d);
	void operator/= (const State &d);
	Coordinate operator+ (const Coordinate &p) const ;
	Coordinate operator- (const Coordinate &p) const ;
	Coordinate operator+ (const State &p) const ;
	Coordinate operator- (const State &p) const ;
	Coordinate operator* (const State &d)const;
	Coordinate operator* (const Coordinate &p)const;
	State operator^ (const Coordinate &p)const;
	Coordinate operator/ (const State &d) const;
	Coordinate operator/ (const Coordinate &d) const;
	bool  operator== (const Coordinate &d)const;
	bool  operator!= (const Coordinate &d)const;
	bool  operator> (const Coordinate &d)const;
	bool  operator< (const Coordinate &d)const;
	bool  operator>= (const Coordinate &d)const;
	bool  operator<= (const Coordinate &d)const;

	//methods
	void split(State & x, State & y, State & z)const;
	void join(const State & x, const State & y, const State & z);
	void normalize();
	void setLength(const State & l);
	double length()const;
	double sum()const;
	double squareSum()const;
};
std::ostream &operator<<(std::ostream &os, const Coordinate &obj);

struct MovingCoordinate:public Coordinate{
	//data
	Coordinate &state;
	Coordinate rate;

	//constructors
	MovingCoordinate();
	MovingCoordinate(const Coordinate& p);
	MovingCoordinate(const MovingCoordinate& p);
	MovingCoordinate(const double &newX, const double &newY, const double &newZ);
	//MovingCoordinate(const StateRate &newX, const StateRate &newY, const StateRate &newZ);

	//overloading aritmatic operators
	MovingCoordinate & operator= (const MovingCoordinate &p) ;
	MovingCoordinate operator+ (const MovingCoordinate &p) const ;
	MovingCoordinate operator- (const MovingCoordinate &p) const ;
	void operator+= (const MovingCoordinate &p);
	void operator-= (const MovingCoordinate &p);
	MovingCoordinate operator* (const State &d)const;
	Coordinate operator/ (const State &d) const;
	bool  operator== (const MovingCoordinate &d)const;
	bool  operator!= (const MovingCoordinate &d)const;


	void split(StateRate & x, StateRate & y, StateRate & z)const;
	void join(const StateRate & x, const StateRate & y, const StateRate & z);
	State length()const;

};
std::ostream &operator<<(std::ostream &os, const MovingCoordinate &obj);


#endif /*COORDINATES_HPP_*/
