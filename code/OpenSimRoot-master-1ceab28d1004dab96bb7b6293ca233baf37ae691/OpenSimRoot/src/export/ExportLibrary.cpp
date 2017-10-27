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

#include "../export/ExportLibrary.hpp"

#include <fstream>
#include "../cli/Messages.hpp"
#include "../cli/Info.hpp"
#include "../engine/SimulaBase.hpp"
#include "../engine/Origin.hpp"
#include "../math/MathLibrary.hpp"
#include <vector>
#include <map>
#include <iomanip>
#include "../export/3dimage/rawimage.hpp"
#include "../export/ExportBaseClass.hpp"
#include "../export/General/GarbageCollection.hpp"
#include "../export/General/PrimeModel.hpp"
#include "../export/General/ProbeAllObjects.hpp"
#include "../export/Text/TabledOutput.hpp"
#include "../export/Text/WriteModel2File.hpp"
#include "../export/VTK/VTU.hpp"


#if _WIN32 || _WIN64
#include "windows.h"
#include "psapi.h"
//#include "Wincon.h"
#else
   static double maxMemUse(0);
#endif

ModuleList modules;///todo runExportMoudules needs to be converted to a proper class which exposes this.
static std::string inputfoldername;
std::string getCurrentInputFolderName(){//hack solang cli does not provide a standard way to retrieve this
	return inputfoldername;
}

void runExportModules(int &exitCode, const std::string &dir){
	try{
		inputfoldername=dir;
		while(!inputfoldername.empty()) {
			if(*(inputfoldername.rbegin())=='/') break;
			inputfoldername.erase(--inputfoldername.end());
		}

		std::cout<<ANSI_Black<<"Running modules: "<<std::flush;
		///@todo we want to loop through list of registered functions.
		//instantiation of export modules,
		///@todo muliple instantiations of rawimage under different names is now possible. Instantiations should really be based on a loop through objects in the simulacontrolparameter list, but we need to recognize the type
		//@todo when not enabled, should be destructed again?
		PrimeModel pm; ProbeAllObjects pao; VTU vtu("VTU"); VTU vtp("vtp"); VTU vtuwire("wireframe.vtu"); RawImage raw1("rasterImage"); Table tbl;GarbageCollection grb;ModelDump mdp;
		if(pm.enabled()) modules.push_back(&pm);
		if(pao.enabled()) modules.push_back(&pao);
		if(tbl.enabled()) modules.push_back(&tbl);
		if(vtu.enabled()) modules.push_back(&vtu);
		if(vtp.enabled()) modules.push_back(&vtp);
		if(vtuwire.enabled()) modules.push_back(&vtuwire);
		if(raw1.enabled()) modules.push_back(&raw1);
		/*if(raw2.enabled()) modules.push_back(&raw2);
		if(raw3.enabled()) modules.push_back(&raw3);
		if(raw4.enabled()) modules.push_back(&raw4);
		if(raw5.enabled()) modules.push_back(&raw5);*/
		//make sure this runs last?
		if(grb.enabled()) modules.push_back(&grb);
		if(mdp.enabled()) modules.push_back(&mdp);

		//initiation of export modules
		for(ModuleList::iterator it(modules.begin()); it!=modules.end(); ++it){
			(*it)->initialize();
		}

		//time loop for export modules
		std::cout<<std::fixed;

		//max endtime
		Time endTime(0);
		for(ModuleList::iterator it(modules.begin()); it!=modules.end(); ++it){
			endTime=maximum(endTime,(*it)->getEndTime());
		}

		//we will push the model in small blocks
		for (Time t(0); t-endTime<=0.05; t+=0.01){
			bool doneSomething(false);
			for(ModuleList::iterator it(modules.begin()); it!=modules.end(); ++it){
				//run modules
				if((*it)->enabled() && (*it)->getCurrentTime()-t<TIMEERROR) {
					(*it)->run((*it)->getCurrentTime());doneSomething=true;//run
					if((*it)->getCurrentTime()>=(*it)->getEndTime()) (*it)->enabled()=false;
					(*it)->getCurrentTime()+=(*it)->getIntervalTime();
					(*it)->getCurrentTime()=minimum((*it)->getCurrentTime(),(*it)->getEndTime());
				}
			}
			if(doneSomething){
//note that on mingw this does not link, but the "linux" code does.
#if _WIN32 || _WIN64
				
				MEMORYSTATUSEX memInfo;
				memInfo.dwLength = sizeof(MEMORYSTATUSEX);
				GlobalMemoryStatusEx(&memInfo);
				//DWORDLONG totalPhysMem = memInfo.ullTotalPhys;
				PROCESS_MEMORY_COUNTERS_EX pmc;
				pmc.cb = sizeof(pmc);
				GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, pmc.cb);
				size_t mem = pmc.WorkingSetSize;
				double m = mem/1000000.;
				std::cout << std::fixed << std::setprecision(1);
				std::cout << " " << t << "/" << endTime << " days. Mem " << (int)floor(m) << " mB. #obj.=" << SimulaBase::numberOfInstantiations
					<< " x64b/obj.=" << (int)(m / ((double)SimulaBase::numberOfInstantiations) * 131072) << "\r" << std::flush;
#else
				/*memory usage status*/
				/*TODO maybe call something smart like a garbage collector or so to reduce memory usage when high?*/
				/*TODO discussion on internet suggests statm is not correct - status is better*/
				std::ifstream statm("/proc/self/statm");
				/*TODO this assumes a page size of 4096, for the correct page size run "getconf PAGE_SIZE" */
				double m;
				statm>>m;statm>>m;m*=4.;m/=1024.;
				statm.close();
				if(m>maxMemUse)maxMemUse=m;
				/*progress indicator*/
				std::cout<<std::fixed<<std::setprecision(1);
				std::cout<<ANSI_save<<t<<"/"<<endTime<<" days. Mem "<<(int)floor(m)<<" mB. #obj.="<<SimulaBase::numberOfInstantiations
					<<" x64b/obj.="<<(int)(m/((double)SimulaBase::numberOfInstantiations)*131072)<<ANSI_restore<<std::flush;
#endif
			}
		}

#if _WIN32_WINNT_WIN7
		std::cout << std::endl;
#endif 
		ANSI_OKmessage
		std::cout<<std::endl<<"Finalizing output:";

		for(ModuleList::iterator it(modules.begin()); it!=modules.end(); ++it){
			(*it)->finalize();
		}

		ANSI_OKmessage



	}	catch (std::exception& error){
		ANSI_FAILmessage
		std::cout << ANSI_Red <<error.what()<< ANSI_Black;
		exitCode=1;
	}
}


//TODO move this to more sensible place
void terminate(){
	// Tell about possible warnings
	writeWarnings();

	//write run time
	writeExecutionTime();

	//reset terminal
	resetTerminal();
}


void writeExecutionTime(){
	std::size_t seconds(clock()/CLOCKS_PER_SEC);
	std::size_t minutes(seconds/60);
	std::size_t hours(minutes/60);
	seconds=seconds-60*minutes;
	minutes=minutes-60*hours;
	std::cout<<std::endl<<"Simulation took (hours:minutes:seconds): "<<hours<<":"<<minutes<<":"<<seconds;
}
void resetTerminal(){
	std::cout<<ANSI_ResetTerminal<<std::endl<<std::flush;
}
