
// RA, 2017-11-19

#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef vector<int> Vec;
typedef vector<Vec> VEC;

VEC read(istream& in) {
	VEC J2I;
	
	string line;
	while (std::getline(in, line)) {
		istringstream buffer(line);
		J2I.push_back(Vec(istream_iterator<int>(buffer), istream_iterator<int>()));
	}
	
	return J2I;
}

VEC transpose(const VEC& A) {
	VEC B;
	
	for (int j = 0; j != A.size(); ++j) {
		for (auto i : A[j]) {
			while (!(B.size() > i)) B.push_back(Vec());
			B[i].push_back(j);
		}
	}
	
	return B;
}

// Symmetric difference
void SD(const Vec& A, const Vec& B, Vec& C) {
	C.clear();
	set_symmetric_difference(A.begin(), A.end(), B.begin(), B.end(), back_inserter(C));
}


int main() {
	string input_file_J2I = "./D/0.txt";
	
	if (!cin.eof()) cin >> input_file_J2I;

	VEC J2I, I2J;
	{
		fstream is(input_file_J2I.c_str(), ios_base::in);

		J2I = read(is);
		I2J = transpose(J2I);
	}
	
	int rank = 0;

	for (int p = 0; p != I2J.size(); ++p) {
		Vec J(I2J[p]);
		
		if (J.empty()) continue;
		
		Vec I(J2I[J[0]]);
		
		for (auto i : I) SD(J, I2J[i], I2J[i]);
		for (auto j : J) SD(I, J2I[j], J2I[j]);
		
		rank++;
	}
	
	cout << rank << endl;
	
	return 0;
}
