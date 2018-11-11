
#ifndef SPHERICALHARMONICS_H
#define SPHERICALHARMONICS_H

namespace SH
{
	// Sample the SH func at position (theta, phi) with band l and parameter m.
	double SampleSHFunc(const int l, const int m, const double theta, const double phi);
}

#endif
