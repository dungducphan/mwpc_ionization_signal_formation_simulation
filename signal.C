double bias_voltage = 1000; // volt
double gap_length = 0.01; // meter
double electron_drift_velocity = 3.0E4; // meter / second (source: https://www.star.bnl.gov/public/tpc/hard/tpcrings/p10_magboltz1b.html)
double electron_mobility = electron_drift_velocity * gap_length / bias_voltage;

void signal_single_pair(double x, double charge, double mobility_factor, std::vector<double>& output) {
  double tcn = x / (electron_mobility * bias_voltage / gap_length);
  double tcp = (gap_length - x) / (mobility_factor * electron_mobility * bias_voltage / gap_length);
  double icnp = charge * electron_mobility * (1 + mobility_factor) * bias_voltage / TMath::Power(gap_length, 2);
  double icp = charge * electron_mobility * mobility_factor * bias_voltage / TMath::Power(gap_length, 2);
}

void signal() {
  double x = 0.5 * gap_length;

}
