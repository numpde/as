
// RA, 2017-11-16

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

#include <omp.h>

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

int main() {
	string input_file_J2I = "./J2I_tmp.txt";
	//cout << input_file_J2I << endl;
	
	if (!cin.eof()) cin >> input_file_J2I;

	VEC J2I, I2J;
	{
		fstream is(input_file_J2I.c_str(), ios_base::in);

		cerr << "Reading file..." << endl;
		J2I = read(is);
		cerr << "Transposing..." << endl;
		I2J = transpose(J2I);
	}
	
	// 
	omp_set_num_threads(2);
	
	cerr << "Constructing len-j / len-i sets..." << endl;
	//
	set < pair<int, int> > lenj, leni;
	//
	#pragma omp parallel
	{
		int thread = omp_get_thread_num();
		
		if (thread == 0) {
			for (auto I : J2I) assert(is_sorted(I.begin(), I.end()));
			
			for (int j = 0; j != J2I.size(); ++j) {
				if (!J2I[j].empty())
					lenj.insert( make_pair(J2I[j].size(), j) );
			}
		}

		if (thread == 1) {
			for (auto J : I2J) assert(is_sorted(J.begin(), J.end()));
			
			for (int i = 0; i != I2J.size(); ++i) {
				if (!I2J[i].empty())
					leni.insert( make_pair(I2J[i].size(), i) );
			}
		}
	}

	cerr << "Computing rank..." << endl;
	//
	int rank = 0;
	//
	while (!lenj.empty() && !leni.empty()) {
		if ((rank % 1000) == 0)
			cerr << lenj.size() << " / " << rank << endl;

		int p = -1;
		int q = -1;
		{
			int i, j;
			#pragma omp parallel
			{
				int thread = omp_get_thread_num();
				
				if (thread == 0) {
					// Index to the shortest J2I[q]
					q = lenj.begin()->second;
					i = shortest(J2I[q], I2J);
				}
				
				if (thread == 1) {
					// Index of the shortest I2J[p]
					p = leni.begin()->second;
					j = shortest(I2J[p], J2I);
				}
			}

			int leniq = I2J[i].size() + J2I[q].size();
			int lenpj = I2J[p].size() + J2I[j].size();

			if (lenpj < leniq) {
				p = i;
			} else {
				q = j;
			}
		}

		Vec I(J2I[q]);
		Vec J(I2J[p]);
		
		#pragma omp parallel
		{
			int thread = omp_get_thread_num();
			
			if (thread == 0) {
				for (int i : I) {
					int sz0 = I2J[i].size();
					SD2(J, I2J[i], I2J[i]);
					int sz1 = I2J[i].size();

					if (sz0 == sz1) continue;
					leni.erase(make_pair(sz0, i));
					if (!sz1) continue;
					leni.insert(make_pair(sz1, i));
				}
			}

			if (thread == 1) {
				for (int j : J) {
					int sz0 = J2I[j].size();
					SD2(I, J2I[j], J2I[j]);
					int sz1 = J2I[j].size();

					if (sz0 == sz1) continue;
					lenj.erase(make_pair(sz0, j));
					if (!sz1) continue;
					lenj.insert(make_pair(sz1, j));
				}
			}
		}
		
		rank++;
	}
	
	cout << rank << endl;
	
	return 0;
}

