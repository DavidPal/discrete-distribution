// C++ implementation of a fast algorithm for generating samples from a
// discrete distribution.
//
// David Pal, December 2015
//
// To compile the program run:
//
//   g++ -Wall -Wextra -Werror -std=c++11 discrete-distribution.cc

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

using std::cout;
using std::endl;

// Stack that does not own the underlying storage.
template<typename T, typename BidirectionalIterator>
class stack_view {
  public:
    stack_view(const BidirectionalIterator base)
      : base_(base), top_(base) { };

    void push(const T& element) {
      *top_ = element;
      ++top_;
    }

    T pop() {
      --top_;
      return *top_;
    }

    bool empty() {
      return top_ == base_;
    }

  private:
    const BidirectionalIterator base_;
    BidirectionalIterator top_;
};

template<typename IntType = int>
class fast_discrete_distribution {
  public:
    typedef IntType result_type;

    fast_discrete_distribution(const std::vector<double>& weights)
      : uniform_distribution_(0.0, 1.0) {
      normalize_weights(weights);
      create_buckets();
    }

    result_type operator()(std::default_random_engine& generator) {
      const double number = uniform_distribution_(generator);
      size_t index = floor(buckets_.size() * number);

      // Fix index
      if (index >= buckets_.size()) index = buckets_.size() - 1;

      const Bucket& bucket = buckets_[index];
      if (number < std::get<2>(bucket))
        return std::get<0>(bucket);
      else
        return std::get<1>(bucket);
    }

    result_type min() const {
      return 0;
    }

    result_type max() const {
      return probabilities_.size();
    }

    std::vector<double> probabilities() const {
      return probabilities_;
    }

    void reset() {
      // Empty
    }

    void PrintBuckets() {
      cout << "buckets.size() = " << buckets_.size() << endl;
      for (auto bucket : buckets_) {
        cout << std::get<0>(bucket) << "  "
             << std::get<1>(bucket) << "  "
             << std::get<2>(bucket) << "  "
             << endl;
      }
    }

  private:
    typedef std::pair<double, size_t> Segment;
    typedef std::tuple<result_type, result_type, double> Bucket;

    void normalize_weights(const std::vector<double>& weights) {
      const double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
      probabilities_.reserve(weights.size());
      for (auto weight : weights) {
        probabilities_.push_back(weight / sum);
      }
    }

    void create_buckets() {
      // Two stacks in one vector.  First stack grows from the begining of the
      // vector. The second stack grow from the end of the vector.
      const size_t N = probabilities_.size();
      std::vector<Segment> segments(N);
      stack_view<Segment, std::vector<Segment>::iterator>
        small(segments.begin());
      stack_view<Segment, std::vector<Segment>::reverse_iterator>
        large(segments.rbegin());

      // Split probabilities into small and large.
      result_type i = 0;
      for (auto probability : probabilities_) {
        if (probability < (1.0 / N)) {
          small.push(Segment(probability, i));
        } else {
          large.push(Segment(probability, i));
        }
        ++i;
      }

      buckets_.reserve(N);

      i = 0;
      while (!small.empty()) {
        const Segment s = small.pop();
        const Segment l = large.pop();

        // Create a mixed bucket
        buckets_.emplace_back(s.second,
                             l.second,
                             s.first + static_cast<double>(i) / N);

        // Calculate the length of the left-over segment
        const double left_over = s.first + l.first - static_cast<double>(i) / N;

        // Re-insert the left-over segment
        if (left_over < (1.0 / N))
          small.push(Segment(left_over, l.second));
        else
          large.push(Segment(left_over, l.second));

        ++i;
      }

      // Create pure buckets
      while (!large.empty()) {
        const Segment l = large.pop();
        buckets_.emplace_back(l.second, l.second, 0.0);
      }
    }

    // Uniform distribution over interval [0,1].
    std::uniform_real_distribution<double> uniform_distribution_;

    // List of probabilities
    std::vector<double> probabilities_;
    std::vector<Bucket> buckets_;
};

void Test(const std::vector<double>& weights, const size_t nrolls) {
  std::default_random_engine generator;
  fast_discrete_distribution<int> distribution(weights);
  distribution.PrintBuckets();

  std::vector<size_t> counts(weights.size(), 0);
  for (size_t i = 0; i < nrolls; ++i) {
    const int number = distribution(generator);
    ++counts[number];
  }

  std::cout << "counts:" << std::endl;
  for (size_t i = 0; i < weights.size(); ++i)
    cout << i << " (" << weights[i] << ") : "
         << std::string(counts[i], '*') << endl;

  cout << endl;
}

int main() {
  Test({0}, 100);
  Test({0, 1e-20, 0}, 100);
  Test({1, 1, 1}, 300);
  Test({1, 1}, 200);
  Test({1}, 100);
  Test({1, 1, 2}, 300);
  Test({1, 0, 2}, 300);

  std::discrete_distribution<int> distribution({10.0, 20.0, 30.0});
  cout << distribution << endl;
  return 0;
}
