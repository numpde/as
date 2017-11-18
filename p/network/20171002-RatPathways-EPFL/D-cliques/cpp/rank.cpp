
// RA, 2017-11-17

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
typedef vector<Vec> VEC; // Mapped matrix 

// Read a stream of the form
//   4 5 7
//   7 8
//   1 3 9
//   ...
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

// Transpose the mapped matrix
VEC transpose(const VEC& A) {
	VEC B;
	
	for (int j = 0; j != A.size(); ++j) {
		for (auto i : A[j]) {
			while (!(B.size() > i))
				B.push_back(Vec());
			B[i].push_back(j);
		}
	}
	
	//#pragma omp for
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

// Find k in I such that I2J[k] is shortest
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
	
	// Use the environment variable OMP_NUM_THREADS, cap to 12
	int max_omp_threads = min(omp_get_max_threads(), 12);
	// Use at least two threads, unless only 1 allowed
	int min_omp_threads = min(2, max_omp_threads);
	// ...but use initially:
	int omp_threads = min_omp_threads;
	omp_set_num_threads(omp_threads);
	assert(omp_get_max_threads() == omp_threads);
	
	// Look for the optimal number of threads?
	bool adaptive_threading = true;
	
	cerr << "OpenMP threads min/max/adaptive: " << min_omp_threads << "/" << max_omp_threads << "/" << adaptive_threading << " threads." << endl;
	
	
	// The sets lenj/leni point to the shortest column/row:
	//
	// The k-th element p = (len, j) in lenj is such that
	// len is the length of the column j
	// and it is the k-th shortest column
	//
	// Similarly for leni
	//
	typedef set < pair<int, int> > LEN;
	LEN lenj, leni;
	//
	#pragma omp parallel
	{
		#pragma omp single nowait
		cerr << "Constructing len-j / len-i sets..." << endl;
		
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
	
	// Essential barrier: wait for lenj/i sets
	#pragma omp barrier

	cerr << "Computing rank..." << endl;
	
	int rank = 0;
	
	//
	vector< vector< LEN::iterator > > erasei(omp_threads), erasej(omp_threads);
	//
	typedef vector< pair< LEN::iterator, pair<int, int> > > VOP;
	vector< VOP > inserti(omp_threads), insertj(omp_threads);
	
	//
	struct C { int p = -1, q = -1, L = -1; };
	vector<C> pqlen(2); // Do not change the initial size here
	
	//
	typedef std::chrono::high_resolution_clock Time;
	typedef std::chrono::duration<double> Duration;
	//
	vector<double> timings(1 + max_omp_threads, 0.0);
	//
	auto t0 = Time::now();
	
	//
	Vec I, J;

	// OUTER LOOP
	while (!lenj.empty() && !leni.empty()) {
		
		// SECTION 1: Thread management
		
		// Adjust the number of openmp threads
		if (adaptive_threading) {
			// timings[n] estimates the runtime of SECTION 2
			// when the number of active threads is n
			timings[omp_threads] = 0.5 * (timings[omp_threads] + Duration(Time::now() - t0).count());
			// Tend to increase the number of threads
			for (int n = min(max_omp_threads, omp_threads + 1); n >= max(min_omp_threads, omp_threads - 1); --n) {
				if (timings[n] < timings[omp_threads]) {
					omp_threads = n;
					omp_set_num_threads(omp_threads);
					cerr << "Switching to " << omp_threads << " thread(s)" << endl;
					break;
				}
			}
			
			// How many pivot candidates to compute?
			// At least 2 in order to check leni and lenj
			pqlen.resize(max(2, omp_threads/2));
			
			erasei.resize(omp_threads); inserti.resize(omp_threads);
			erasej.resize(omp_threads); insertj.resize(omp_threads);
			
			t0 = Time::now();
		}
		
		// SECTION 2: Parallel computation
		
		#pragma omp parallel 
		{
			int thread = omp_get_thread_num();
			
			// INNER LOOP
			while (!lenj.empty() && !leni.empty()) {
				
				#pragma omp parallel for
				for (int k = 0; k < pqlen.size(); ++k) {
					pqlen[k].L = J2I.size() + I2J.size();
					
					// At even k
					if (((k % 2) == 0) && (k < lenj.size())) {
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
					
					// At odd k
					if (((k % 2) == 1) && (k < leni.size())) {
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
				
				// Essential barrier: wait for pqlen
				#pragma omp barrier
			
				// Cannot set "nowait": synchronize I and J
				#pragma omp single 
				{
					// Find the best pivot
					C bestc = pqlen[0];
					for (auto& c : pqlen) if (c.L < bestc.L) bestc = c;

					// Nonzeros in the pivot column/row
					I = J2I[bestc.q];
					J = I2J[bestc.p];
				}

				// Manipulate the matrix
				// Compute erasei/inserti for leni update
				{
					erasei[thread].clear(); inserti[thread].clear();
					
					#pragma omp for nowait
					for (int ii = 0; ii < I.size(); ++ii) {
						int i = I[ii];
						
						int sz0 = I2J[i].size();
						SD2(J, I2J[i], I2J[i]);
						int sz1 = I2J[i].size();

						if (sz0 == sz1) continue;
						erasei[thread].push_back(leni.find(make_pair(sz0, i)));
						if (!sz1) continue;
						auto p = make_pair(sz1, i);
						inserti[thread].push_back(make_pair(leni.upper_bound(p), p));
					}
				}

				// Manipulate the matrix-transpose
				// Compute erasej/insertj for lenj update
				{
					erasej[thread].clear(); insertj[thread].clear();
					
					#pragma omp for nowait
					for (int jj = 0; jj < J.size(); ++jj) {
						int j = J[jj];
						
						int sz0 = J2I[j].size();
						SD2(I, J2I[j], J2I[j]);
						int sz1 = J2I[j].size();

						if (sz0 == sz1) continue;
						erasej[thread].push_back(lenj.find(make_pair(sz0, j)));
						if (!sz1) continue;
						auto p = make_pair(sz1, j);
						insertj[thread].push_back(make_pair(lenj.upper_bound(p), p));
					}
				}
			
				// Essential barrier: wait for inserti/j, erasei/j
				#pragma omp barrier
				
				#pragma omp single nowait
				{
					// Cannot reverse the order of insertion and erasure
					for (auto& v : inserti) for (auto& p : v) leni.insert(p.first, p.second);
					for (auto& v : erasei)  for (auto& i : v) leni.erase(i);
				}
					
				#pragma omp single nowait
				{
					// Cannot reverse the order of insertion and erasure
					for (auto& v : insertj) for (auto& p : v) lenj.insert(p.first, p.second);
					for (auto& v : erasej)  for (auto& j : v) lenj.erase(j);
				}
				
				// Cannot set "nowait": essential (implicit) barrier
				#pragma omp single
				rank++;
				
				if ((rank % 1000) == 0) {
					#pragma omp single
					{
						cerr << "#I2J x #J2I: " << I2J.size() << "x" << J2I.size() << " | ";
						cerr << "#lenj: " << lenj.size() << " | " << "rank: " << rank << endl;
					}
					
					// leave inner loop
					break; 
				}
			} // while
		} // omp parallel
	} // while
	
	// Communicate the result
	cout << rank << endl;
	
	return 0;
}

