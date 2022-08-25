#include "TMath.h"

#include <iostream>
#include <vector>

double bias_voltage = 1000; // volt
double gap_length = 0.01; // meter
double electron_drift_velocity = 3.0E4; // meter / second (source: https://www.star.bnl.gov/public/tpc/hard/tpcrings/p10_magboltz1b.html)
double electron_mobility = electron_drift_velocity * gap_length / bias_voltage;

void signal_single_pair(double x, double charge, double mobility_factor, std::vector<double>& output) {
  double tcn = x / (electron_mobility * bias_voltage / gap_length);
  double tcp = (gap_length - x) / (mobility_factor * electron_mobility * bias_voltage / gap_length);

  double icnp = charge * electron_mobility * (1 + mobility_factor) * bias_voltage / TMath::Power(gap_length, 2);
  double icp = charge * electron_mobility * mobility_factor * bias_voltage / TMath::Power(gap_length, 2);
  output.push_back(tcn / 1E-9);
  output.push_back(tcp / 1E-9);
  output.push_back(icnp);
  output.push_back(icp);
}

void signal() {
/* 
    double x = 0.1 * gap_length;
    double charge = 1.6E-19;
    double mobility_factor = 2. / 1000.;
    std::vector<double> output;
    signal_single_pair(x, charge, mobility_factor, output);

    for (auto& elem : output)
        std::cout << elem << std::endl;

    std::cout << output.at(2) * output.at(0) + output.at(3) * (output.at(1) - output.at(0))  << std::endl;
*/


    double total_Edep = 3.913E-3; // MeV
    double ar_ionE = 26.4; // eV
    int total_pair = (int) (total_Edep * 1E6 / ar_ionE);
    double average_pair_distance = gap_length / (double) total_pair;


    std::vector<std::vector<double>> single_pair_current;
    for (unsigned int i = 0; i < total_pair; i++) {
        double x = average_pair_distance * (0.5 + (double) i);
        double charge = 1.6E-19;
        // double mobility_factor = 2. / 1000.;
        double mobility_factor = 2. / 100.;
        std::vector<double> output;
        signal_single_pair(x, charge, mobility_factor, output);
        single_pair_current.push_back(output);
    }

    std::vector<std::pair<double, double>> total_current;
    for (unsigned int i = 0; i < 10001; i++) {
        double current = 0;
        for (auto& elem : single_pair_current) {
            if (i <= elem.at(0)) {
                current += elem.at(2);
            } else if (i > elem.at(0) && i <= elem.at(1)) {
                current += elem.at(3);
            }
        }
        if (current != 0) std::cout << Form("%05i : %15.4f", i, current) << std::endl;
        total_current.push_back(std::make_pair((double)i, current));
    }

    TGraph* gr = new TGraph();
    for (auto& elem : total_current) {
        gr->AddPoint(elem.first, elem.second);
    }
    gr->Draw("APL");
}
