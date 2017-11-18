
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

//
struct Pivot { 
	int p = -1;  // row number
	int q = -1;  // col number
	
	int nI = -1; // nnz in pivot row
	int nJ = -1; // nnz in pivot col
	
	Pivot t() const {
		Pivot Q;
		
		Q.p = q; Q.nI = nJ;
		Q.q = p; Q.nJ = nI;
		
		return Q;
	}
	
	// What is a better (smaller) pivot?
	bool operator<(const Pivot& Q) { 
		return (nI * nJ) < (Q.nI * Q.nJ); 
	}
	
};
//
std::ostream& operator<<(std::ostream& out, const Pivot& P) {
	out << "Pivot[(" << P.p << ", " << P.q << "), " << P.nJ << "x" << P.nI << "]";
	return out;
}
	
// Find i in J2I[j] such that I2J[i] is shortest
// Make a Pivot with (p, q) = (i, j)
Pivot shortest(int j, const VEC& J2I, const VEC& I2J) {
	Pivot P;
	VEC::value_type I(J2I[j]);

	P.q = j;
	P.nI = I.size();
	
	P.p = I[0];
	P.nJ = I2J[P.p].size();
	
	for (auto i : I) {
		if (I2J[i].size() < P.nJ) {
			P.p = i;
			P.nJ = I2J[i].size();
		}
	}
	
	return P;
}

//

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
	
	// Use the environment variable OMP_NUM_THREADS, cap to 24
	int max_omp_threads = min(omp_get_max_threads(), 24);
	// Use at least two threads, unless only 1 is allowed
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
	#pragma omp parallel sections
	{
		#pragma omp section
		cerr << "Constructing len-j / len-i sets..." << endl;
		
		#pragma omp section
		{
			for (auto I : J2I) assert(is_sorted(I.begin(), I.end()));
		}
		
		#pragma omp section
		{
			for (auto J : I2J) assert(is_sorted(J.begin(), J.end()));
		}
		
		#pragma omp section
		{
			for (int j = 0; j != J2I.size(); ++j) {
				if (!J2I[j].empty())
					lenj.insert( make_pair(J2I[j].size(), j) );
			}
		}

		#pragma omp section
		{
			for (int i = 0; i != I2J.size(); ++i) {
				if (!I2J[i].empty())
					leni.insert( make_pair(I2J[i].size(), i) );
			}
		}
	}

	cerr << "Computing rank..." << endl;
	
	int rank = 0;
	
	//
	vector< vector< LEN::iterator > > erasei(omp_threads), erasej(omp_threads);
	//
	typedef vector< pair< LEN::iterator, pair<int, int> > > VOP;
	vector< VOP > inserti(omp_threads), insertj(omp_threads);
	
	// Do not change the default initial size here:
	vector<Pivot> pivots(2);
	
	//
	typedef std::chrono::high_resolution_clock Time;
	typedef std::chrono::duration<double> Duration;
	//
	vector<double> timings(1 + max_omp_threads, 0.0);
	//
	// Used for selecting the number of threads:
	auto t0 = Time::now();
	
	double time_a_pivoting = 0;
	double time_b_flipping = 0;
	double time_c_ordering = 0;
	
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
			
			{
				cerr << "Timings. ";
				cerr << "A: " << time_a_pivoting << "s" << ", "; 
				cerr << "B: " << time_b_flipping << "s" << ", "; 
				cerr << "C: " << time_c_ordering << "s" << endl; 
			}
			
			// How many pivot candidates to compute?
			// At least 2 in order to check leni and lenj
			// (Comment out to use the default defined above)
			//pivots.resize(max(2, omp_threads/2));
			
			assert(pivots.size());
			
			erasei.resize(omp_threads); inserti.resize(omp_threads);
			erasej.resize(omp_threads); insertj.resize(omp_threads);
			
			cerr << endl;
			
			// Used for selecting the number of threads:
			t0 = Time::now();
		}
		
		// SECTION 2: Parallel computation
		
		#pragma omp parallel 
		{
			int thread_id = omp_get_thread_num();
			
			// INNER LOOP
			while (!lenj.empty() && !leni.empty()) {
				
				//
				// Step 1: Find a pivot
				//
				
				auto ta = Time::now();
	
				// Need to find at least one pivot candidate
				assert(pivots.size());
				
				// Essential implicit barrier: wait for 'pivots'
				#pragma omp for 
				for (int k = 0; k < pivots.size(); ++k) {
					// At even k
					if (((k % 2) == 0) && (k < lenj.size())) {
						// Find index to the k-th shortest J2I[q]
						auto j = lenj.begin();
						advance(j, k);
						// Best pivot in this column
						pivots[k] = shortest(j->second, J2I, I2J);
					}
					
					// At odd k
					if (((k % 2) == 1) && (k < leni.size())) {
						// Find index of the k-th shortest I2J[p]
						auto i = leni.begin();
						advance(i, k);
						// Best pivot in this row
						pivots[k] = shortest(i->second, I2J, J2I).t(/* (!) */);
					}
				}
				
				// Essential implicit barrier: wait for 'I' and 'J'
				#pragma omp single
				{
					// Find the best pivot
					Pivot best = pivots[0];
					for (auto& P : pivots) if (P < best) best = P;
					
					//cout << best << endl;

					// Nonzeros in the pivot ...
					I = J2I[best.q]; //  ... column
					J = I2J[best.p]; //  ... row
				}
				
				#pragma omp single
				time_a_pivoting += Duration(Time::now() - ta).count();
				
				//
				// Step 2: Flipping matrix
				//

				auto tb = Time::now();
				
				// Manipulate the matrix
				// Compute erasei/inserti for leni update
				{
					erasei[thread_id].clear(); inserti[thread_id].clear();
					
					#pragma omp for nowait
					for (int ii = 0; ii < I.size(); ++ii) {
						int i = I[ii];
						
						int sz0 = I2J[i].size();
						SD2(J, I2J[i], I2J[i]);
						int sz1 = I2J[i].size();

						if (sz0 == sz1) continue;
						erasei[thread_id].push_back(leni.find(make_pair(sz0, i)));
						if (!sz1) continue;
						auto p = make_pair(sz1, i);
						inserti[thread_id].push_back(make_pair(leni.upper_bound(p), p));
					}
				}

				// Manipulate the matrix-transpose
				// Compute erasej/insertj for lenj update
				{
					erasej[thread_id].clear(); insertj[thread_id].clear();
					
					#pragma omp for nowait
					for (int jj = 0; jj < J.size(); ++jj) {
						int j = J[jj];
						
						int sz0 = J2I[j].size();
						SD2(I, J2I[j], J2I[j]);
						int sz1 = J2I[j].size();

						if (sz0 == sz1) continue;
						erasej[thread_id].push_back(lenj.find(make_pair(sz0, j)));
						if (!sz1) continue;
						auto p = make_pair(sz1, j);
						insertj[thread_id].push_back(make_pair(lenj.upper_bound(p), p));
					}
				}
				
				// Essential barrier: wait for 'inserti/j' and 'erasei/j'
				#pragma omp barrier
				
				#pragma omp single
				time_b_flipping += Duration(Time::now() - tb).count();
				
				//
				// Step 3: Ordering lenghts
				//
				
				auto tc = Time::now();
				
				#pragma omp single nowait
				{
					// Note: Cannot reorder insertion and erasure!
					// (iterator invalidation possible after 'erase')
					for (auto& v : inserti) for (auto& p : v) leni.insert(p.first, p.second);
					for (auto& v : erasei)  for (auto& i : v) leni.erase(i);
				}
						
				#pragma omp single nowait
				{
					// Note: Cannot reorder insertion and erasure!
					// (iterator invalidation possible after 'erase')
					for (auto& v : insertj) for (auto& p : v) lenj.insert(p.first, p.second);
					for (auto& v : erasej)  for (auto& j : v) lenj.erase(j);
				}
				
				#pragma omp single nowait
				rank++;
				
				// Essential barrier: wait for 'leni/j' and 'rank'
				#pragma omp barrier
				
				#pragma omp single
				time_c_ordering += Duration(Time::now() - tc).count();
				
				
				// Note: no thread can reach here before 'rank'
				//       is updated by virtue of the omp barrier
				if ((rank % 1000) == 0) {
					#pragma omp single
					{
						cerr << "#I2J x #J2I: " << I2J.size() << "x" << J2I.size() << " | ";
						cerr << "#lenj: " << lenj.size() << " | " << "rank: " << rank << endl;
					}
					
					// leave inner loop and the parallel section
					break; 
				}
				
			} // while
		} // omp parallel
	} // while
	
	// Communicate the result
	cout << rank << endl;
	
	return 0;
}

