/* * Copyright (c) 2016 Regents of the University of California. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* 3. The names of its contributors may not be used to endorse or promote
*    products derived from this software without specific prior written
*    permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* *********************************************************************************************** *
* CARLsim
* created by: (MDR) Micah Richert, (JN) Jayram M. Nageswaran
* maintained by:
* (MA) Mike Avery <averym@uci.edu>
* (MB) Michael Beyeler <mbeyeler@uci.edu>,
* (KDC) Kristofor Carlson <kdcarlso@uci.edu>
* (TSC) Ting-Shuo Chou <tingshuc@uci.edu>
* (HK) Hirak J Kashyap <kashyaph@uci.edu>
*
* CARLsim v1.0: JM, MDR
* CARLsim v2.0/v2.1/v2.2: JM, MDR, MA, MB, KDC
* CARLsim3: MB, KDC, TSC
* CARLsim4: TSC, HK
* CARLsim5: HK, JX, KC
* CARLsim6: LN, JX, KC, KW
*
* CARLsim available from http://socsci.uci.edu/~jkrichma/CARLsim/
* Ver 12/31/2016
*/

// include CARLsim user interface

#include <carlsim.h>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

void printToFile(vector<vector<int>> &spikes, vector<vector<float>> &voltages, vector<vector<float>> &currents) {
	int dt = 1;
	ofstream outfile;
	// spikes
	outfile.open("/home/nmsutton/Dropbox/CompNeuro/gmu/research/sim_project/code/moradi_code/moradi_stp_compare/stp_compare_spikes.csv", ios::out | ios::trunc);
	for (int i = 0; i < spikes[0].size(); i++) {
		if (spikes[0][i] != 0) {
			outfile << spikes[0][i] << ", ";
		}
	}
	outfile.close();
	// volts, current
	outfile.open("/home/nmsutton/Dropbox/CompNeuro/gmu/research/sim_project/code/moradi_code/moradi_stp_compare/stp_compare_results.csv", ios::out | ios::trunc);
	for (int i = 0; i < voltages[0].size(); i++) {
		for (int j = 0; j < dt; j++) {
			outfile << voltages[0][i] << endl << endl << endl;
			outfile << currents[0][i] << endl << endl << endl;
		}
	}
	outfile.close();
}

int main() {
	// ---------------- CONFIG STATE -------------------
	int numGPUs = 1;	
	int randSeed = 42;
	CARLsim sim("moradi_stp_compare", GPU_MODE, USER, numGPUs, randSeed);
	#include "../generate_config_state.cpp" 	// configure the network
	sim.setIntegrationMethod(RUNGE_KUTTA4, 40);
	NeuronMonitor* NrnMonAA = sim.setNeuronMonitor(Axo_Axonic,"DEFAULT");
	NeuronMonitor* NrnMonPyr = sim.setNeuronMonitor(Pyramidal,"DEFAULT");	
	// ---------------- SETUP STATE -------------------
	sim.setupNetwork();
	SpikeMonitor* SpkMonAA = sim.setSpikeMonitor(Axo_Axonic,"DEFAULT");
	SpikeMonitor* SpkMonPyr = sim.setSpikeMonitor(Pyramidal,"DEFAULT");
	// ---------------- RUN STATE -------------------
	sim.setExternalCurrent(Axo_Axonic, 150);
	SpkMonAA->startRecording(); SpkMonPyr->startRecording(); 
	NrnMonAA->startRecording(); NrnMonPyr->startRecording();
	// run
	for (int i=0; i<10; i++) {sim.runNetwork(0,100, true);}
	SpkMonAA->stopRecording(); SpkMonPyr->stopRecording(); 
	NrnMonAA->stopRecording(); NrnMonPyr->stopRecording();
	SpkMonAA->print(false); SpkMonPyr->print(false);
	vector<vector<int>> spikes = SpkMonAA->getSpikeVector2D();
	vector<vector<float>> voltages = NrnMonPyr->getVectorV();
	vector<vector<float>> currents = NrnMonPyr->getVectorI();
	printToFile(spikes, voltages, currents);
	return 0;
}
