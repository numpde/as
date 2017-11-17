
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
#include <chrono>
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

template <class T>
void clear_each(vector<T>& V) {
	for (auto& v : V) v.clear();
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
	
	// Use MAX_OMP_THREADS to determine the maximal number
	int max_omp_threads = omp_get_max_threads();
	int min_omp_threads = 2;
	// but use initially this number of threads
	int omp_threads = min_omp_threads;
	omp_set_num_threads(omp_threads);
	assert(omp_get_max_threads() == omp_threads);
	
	cerr << "Using maximally " << max_omp_threads << " threads." << endl;
	
	
	cerr << "Constructing len-j / len-i sets..." << endl;
	//
	set < pair<int, int> > lenj, leni;
	//
	#pragma omp parallel
	{
		int thread = omp_get_thread_num();
		
		if (thread == (0 % omp_threads)) {
			for (auto I : J2I) assert(is_sorted(I.begin(), I.end()));
			
			for (int j = 0; j != J2I.size(); ++j) {
				if (!J2I[j].empty())
					lenj.insert( make_pair(J2I[j].size(), j) );
			}
		}

		if (thread == (1 % omp_threads)) {
			for (auto J : I2J) assert(is_sorted(J.begin(), J.end()));
			
			for (int i = 0; i != I2J.size(); ++i) {
				if (!I2J[i].empty())
					leni.insert( make_pair(I2J[i].size(), i) );
			}
		}
	}
	
	#pragma omp barrier

	cerr << "Computing rank..." << endl;
	
	int rank = 0;
	
	typedef vector< pair<int, int> > VOP;
	vector< VOP > erasei, inserti, erasej, insertj;
	
	for (int t = 0; t < max_omp_threads; ++t) {
		erasei.push_back(VOP()); inserti.push_back(VOP());
		erasej.push_back(VOP()); insertj.push_back(VOP());
	}
	
	// The initial length should be at least one
	// Suggested length: omp_threads
	struct C { int p, q, L; };
	vector<C> pqlen(max_omp_threads);
	
	typedef std::chrono::high_resolution_clock Time;
	typedef std::chrono::duration<double> Duration;

	vector<double> timings(1 + max_omp_threads, 0.0);
	
	auto t0 = Time::now();

	while (!lenj.empty() && !leni.empty()) {
		if ((rank % 1000) == 0) {
			cerr << lenj.size() << " / " << rank << endl;
			
			// Adjust the number of openmp threads
			timings[omp_threads] = 0.5 * (timings[omp_threads] + Duration(Time::now() - t0).count());
			for (int n = min(max_omp_threads, omp_threads + 1); n >= max(min_omp_threads, omp_threads - 1); --n) {
				if (timings[n] < timings[omp_threads]) {
					omp_threads = n;
					omp_set_num_threads(omp_threads);
					cerr << "Switching to " << omp_threads << " thread(s)" << endl;
					break;
				}
			}
			
			t0 = Time::now();
		}
		
		#pragma omp parallel for
		for (int k = 0; k < pqlen.size(); ++k) {
			pqlen[k].L = J2I.size() + I2J.size();
			
			if ((k < lenj.size()) && ((k % 2) == 0)) {
				// Index to the k-th shortest J2I[q]
				auto j = lenj.begin();
				advance(j, k);
				C c;
				c.q = j->second;
				c.p = shortest(J2I[c.q], I2J);
				c.L = I2J[c.p].size() + J2I[c.q].size();
				//
				if (c.L < pqlen[k].L) pqlen[k] = c;
			}
			
			if ((k < leni.size()) && ((k % 2) == 1)) {
				// Index of the k-th shortest I2J[p]
				auto i = leni.begin();
				advance(i, k);
				C c;
				c.p = i->second;
				c.q = shortest(I2J[c.p], J2I);
				c.L = I2J[c.p].size() + J2I[c.q].size();
				//
				if (c.L < pqlen[k].L) pqlen[k] = c;
			}
		}
		
		C bestc = pqlen[0];
		for (auto& c : pqlen) if (c.L < bestc.L) bestc = c;

		Vec I(J2I[bestc.q]);
		Vec J(I2J[bestc.p]);
		
		#pragma omp parallel
		{
			int thread = omp_get_thread_num();

			erasei[thread].clear(); inserti[thread].clear();
			erasej[thread].clear(); insertj[thread].clear();
			
			#pragma omp for
			for (int ii = 0; ii < I.size(); ++ii) {
				int i = I[ii];
				
				int sz0 = I2J[i].size();
				SD2(J, I2J[i], I2J[i]);
				int sz1 = I2J[i].size();

				if (sz0 == sz1) continue;
				erasei[thread].push_back(make_pair(sz0, i));
				if (!sz1) continue;
				inserti[thread].push_back(make_pair(sz1, i));
			}
			
			//sort(inserti[thread].begin(), inserti[thread].end());

			#pragma omp for
			for (int jj = 0; jj < J.size(); ++jj) {
				int j = J[jj];
				
				int sz0 = J2I[j].size();
				SD2(I, J2I[j], J2I[j]);
				int sz1 = J2I[j].size();

				if (sz0 == sz1) continue;
				erasej[thread].push_back(make_pair(sz0, j));
				if (!sz1) continue;
				insertj[thread].push_back(make_pair(sz1, j));
			}
			
			//sort(insertj[thread].begin(), insertj[thread].end());
		}
		
		// omp barrier
		
		#pragma omp parallel
		{
			int thread = omp_get_thread_num();
			
			if (thread == (0 % omp_threads)) {
				for (int t = 0; t < omp_threads; ++t) for (auto& p : erasei[t]) leni.erase(p);
				for (int t = 0; t < omp_threads; ++t) leni.insert(inserti[t].begin(), inserti[t].end());
			}
			
			if (thread == (1 % omp_threads)) {
				for (int t = 0; t < omp_threads; ++t) for (auto& p : erasej[t]) lenj.erase(p);
				for (int t = 0; t < omp_threads; ++t) lenj.insert(insertj[t].begin(), insertj[t].end());
			}
		}
		
		rank++;
	}
	
	cout << rank << endl;
	
	return 0;
}

