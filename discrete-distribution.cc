// C++ implementation of a fast algorithm for generating samples from a
// discrete distribution.
//
// David Pal, December 2015
//
// To compile the program run:
//
//   g++ -Wall -Wextra -Werror -std=c++11 discrete-distribution.cc

#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

class DiscreteDistribution {
	public:
		DiscreteDistribution(const vector<double>& weights) {
			PreProcess(weights);
		}


	private:
		struct Bucket {
			size_t color1 = 0;
			size_t color2 = 0;
			double threshold = 0.0;
		};

		void PreProcess(const vector<double>& weights) {
			const double sum = std::accumulate(weights.cbegin(), weights.cend(), 0.0);
			cout << "weights = " << weights.size() << endl;
			cout << "sum = " << sum << endl;

			// Reserve capacity to avoid reallocations.
			// This might take twice as much memory as is actually needed.
			const size_t N = weights.size();
			vector<double> small, big;
			small.reserve(N);
			big.reserve(N);

			// Normalize the probabilities.
			for (size_t i = 0; i < weights.size(); ++i) {
				const double probability = weights[i] / sum;
				if (probability < (1.0 / N))
					small.push_back(probability);
				else
					big.push_back(probability);
			}

			cout << "small.size() = " << small.size() << endl;
			cout << "big.size() = " << big.size() << endl;

			for (size_t i = 0; i < weights.size(); ++i) {
				//
			}
		}

		vector<Bucket> buckets;
};

int main() {
	DiscreteDistribution d({1, 2, 3});
}
