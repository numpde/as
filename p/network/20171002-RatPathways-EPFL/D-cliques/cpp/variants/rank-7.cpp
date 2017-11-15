
// RA, 2017-11-09

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
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
typedef vector<Vec> VEC;
typedef set < pair<int, int> > Len;

VEC read(istream& in) {
	VEC J2I;
	
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

VEC transpose(const VEC& A) {
	VEC B;
	
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

int shortest(const Vec& I, const VEC& I2J) {
	int k = I[0];
	int L = I2J[k].size();
	for (auto i : I) {
		if (I2J[i].size() >= L) continue;
		L = I2J[i].size();
		k = i;
	}
	return k;
}

int gauss_step(const Vec& I, const Vec& J, VEC& J2I, Len& lenj, bool update = true) {
	int dnnz = 0;

	Vec tmp;
	for (auto j : J) {
		Vec& result = (update ? J2I[j] : tmp);
		
		int sz0 = J2I[j].size();
		SD2(I, J2I[j], result);
		int sz1 = result.size();
		
		dnnz += (sz1 - sz0);

		if (!update) continue;
		
		if (sz0 == sz1) continue;
		lenj.erase(make_pair(sz0, j));
		if (!sz1) continue;
		lenj.insert(make_pair(sz1, j));
	}
	
	return dnnz;
}

void sparsify(Vec P, VEC& I2J, VEC& J2I) {
	if (P.size() <= 1) return;
	
// 	int nnz = 0;
// 	for (auto& J : I2J) nnz += J.size();
// 	cerr << P.size() << " / " << nnz << endl;
	
	int p = P.back();
	P.pop_back();
	
	Vec Q(I2J[p]);
	
	Vec I;
	for (auto i : P) {
		if (I2J[i].empty()) continue;
		
		Vec R;
		SD2(Q, I2J[i], R);
		if (R.size() < I2J[i].size()) I.push_back(i);
		
// 		cerr << "Q:      " << Q      << endl;
// 		cerr << "I2J[i]: " << I2J[i] << endl;
// 		cerr << "R:      " << R      << endl;
	}
	
	for (auto i : I) SD(Q, I2J[i], I2J[i]);
	for (auto j : Q) SD(I, J2I[j], J2I[j]);
	
	sparsify(P, I2J, J2I);
}

int main() {
	VEC J2I, I2J;
	{
		string input_file_J2I = "./J2I_tmp.txt";
		
		if (!cin.eof()) cin >> input_file_J2I;
		cerr << "Input file: " << input_file_J2I << endl;

		fstream is(input_file_J2I.c_str(), ios_base::in);

		cerr << "Reading file..." << endl;
		J2I = read(is);
// 		cerr << "Shuffling..." << endl;
// 		random_shuffle(J2I.begin(), J2I.end());
		cerr << "Transposing..." << endl;
		I2J = transpose(J2I);
	}
	
// 	{
// 		cerr << "Sparsifying..." << endl;
// 		
// // 		Vec J;
// // 		for (int j = 0; j != J2I.size(); ++j) J.push_back(j);
// // 		sparsify(J, J2I, I2J);
// 		
// 		Vec I;
// 		for (int i = 0; i != I2J.size(); ++i) I.push_back(i);
// 		sparsify(I, I2J, J2I);
// 		
// 	}
	
	//for (auto I : J2I) assert(is_sorted(I.begin(), I.end()));
	//for (auto J : I2J) assert(is_sorted(J.begin(), J.end()));
	
	int rank = 0;
	
	for (auto& I : J2I) {
		int nnz = 0;
		for (auto& J : I2J) nnz += J.size();
		
		if ((rank % 1000) == 0) {
			cerr << "size: " << I2J.size() << "x" << J2I.size() << " | ";
			cerr << "nnz: "  << nnz  << " | ";
			cerr << "rank: " << rank << endl;
		}
		
		if (I.empty()) continue;
		
		// Make a copy
		Vec J(I2J[*I.begin()]);
		
		for (auto i : I) SD(J, I2J[i], I2J[i]);
		for (auto j : J) SD(I, J2I[j], J2I[j]);
		
		
		int dnnz = 0;
		for (auto& J : I2J) dnnz -= J.size();
		sparsify(I, I2J, J2I);
		sparsify(J, J2I, I2J);
		for (auto& J : I2J) dnnz += J.size();
// 		if (dnnz != 0) cerr << "dnnz: " << dnnz << endl;
		
		rank++;
	}
	
	cout << rank << endl;
	
	return 0;
}
