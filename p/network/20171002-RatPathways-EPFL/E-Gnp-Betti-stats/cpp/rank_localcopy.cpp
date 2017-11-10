
// RA, 2017-11-09

#include <algorithm>
#include <iostream>
#include <iterator>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <set>

// https://stackoverflow.com/questions/10750057/how-to-print-out-the-contents-of-a-vector
template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
	out << '[';
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
	if (!v.empty()) out << "\b\b";
	out << "]";
	return out;
}


using namespace std;

typedef vector<int> Vec;
typedef vector<Vec> map;

map read(istream& in) {
	map J2I;
	
	string line;
	while (std::getline(in, line)) {
		istringstream buffer(line);
		
		J2I.push_back(
			Vec(
				istream_iterator<int>(buffer),
				istream_iterator<int>()
			)
		);
		
		auto& I = *J2I.rbegin();
		sort(I.begin(), I.end());
	}
	
	return J2I;
}

map transpose(const map& A) {
	map B;
	
	for (int j = 0; j != A.size(); ++j) {
		for (auto i : A[j]) {
			while (!(B.size() > i))
				B.push_back(Vec());
			B[i].push_back(j);
		}
	}
	
	for (auto J : B) sort(J.begin(), J.end());
	
	return B;
}

// Symmetric difference assuming sorted vectors
void SD2(const Vec& A, const Vec& B, Vec& C) {
	auto a = A.begin();
	auto b = B.begin();
	
	Vec tmp;
	
	while ((a != A.end()) && (b != B.end())) {
		if (*b < *a) {
			tmp.push_back(*b++);
			continue;
		}
		
		if (*a < *b) {
			tmp.push_back(*a++);
			continue;
		}
		
		if (*a == *b) {
			a++;
			b++;
			continue;
		}
	}
	
	while (a != A.end()) tmp.push_back(*a++);
	while (b != B.end()) tmp.push_back(*b++);
	
	tmp.swap(C);
}

// Symmetric difference
void SD(const Vec& A, const Vec& B, Vec& C) {
	assert(is_sorted(B.begin(), B.end()));
	
	Vec tmp;
	set_symmetric_difference(A.begin(), A.end(), B.begin(), B.end(), back_inserter(tmp));
	tmp.swap(C);
}

int main() {
	string input_file_J2I = "./J2I_tmp.txt";
	//cout << input_file_J2I << endl;
	
	if (!cin.eof()) cin >> input_file_J2I;

	auto is = fstream(input_file_J2I.c_str(), ios_base::in);

	map J2I = read(is);
	map I2J = transpose(J2I);
	//cout << I2J[12] << endl;
	
	//for (auto I : J2I) assert(is_sorted(I.begin(), I.end()));
	//for (auto J : I2J) assert(is_sorted(J.begin(), J.end()));
	
	int Q = J2I.size();
	int rank = 0;
	
	//*/
	for (int q = 0; q != Q; ++q) {
		Vec I(J2I[q]);
		
		if (I.empty()) continue;
		
		Vec J(I2J[I[0]]);
		
		for (auto i : I) SD2(J, I2J[i], I2J[i]);
		for (auto j : J) SD2(I, J2I[j], J2I[j]);
		
		rank++;
	}
	/*/
	for (int p = 0; p != P; ++p) {
		Vec J(I2J[p]);
		
		if (J.empty()) continue;
		
		int q = *J.begin();
		Vec I(J2I[q]);
		
		for (auto i : I) SD(J, I2J[i], I2J[i]);
		
		for (auto j : J) SD(I, J2I[j], J2I[j]);
		
		rank++;
	}
	//*/
	
	cout << rank << endl;
	
	return 0;
}

/*/
int main() {
	string input_file_J2I = "./J2I_tmp.txt";
	//cout << input_file_J2I << endl;
	
	if (!cin.eof()) cin >> input_file_J2I;

	auto is = fstream(input_file_J2I.c_str(), ios_base::in);

	map J2I = read(is);
	map I2J = transpose(J2I);
	//cout << I2J[12] << endl;
	
	set<int> nzj;
	for (int j = 0; j < J2I.size(); ++j)
		if (!J2I[j].empty())
			nzj.insert(j);
	
	set<int> nzi;
	for (int i = 0; i < I2J.size(); ++i)
		if (!I2J[i].empty())
			nzi.insert(i);
	
	int rank = 0;
	
	while (!nzj.empty() && !nzi.empty()) {
		int p = -1;
		int q = -1;
		
		{
			int qa = *nzj.begin();
			Vec& IA = J2I[qa];
			int pa = *IA.begin();
			Vec& JA = I2J[pa];
			
			int pb = *nzi.begin();
			Vec& JB = I2J[pb];
			int qb = *JB.begin();
			Vec& IB = J2I[qb];
			
			if ((IA.size() + JA.size()) < (IB.size() + JB.size())) {
				p = pa;
				q = qa;
			} else {
				p = pb;
				q = qb;
			}
		}
		
		Vec J(I2J[p]);
		Vec I(J2I[q]);
		
		for (auto i : I) {
			SD(J, I2J[i], I2J[i]);
			if (I2J[i].empty()) nzi.erase(i);
		}
		
		for (auto j : J) {
			SD(I, J2I[j], J2I[j]);
			if (J2I[j].empty()) nzj.erase(j);
		}
		
		rank++;
	}
	
	cout << rank << endl;
	
	return 0;
}
//*/