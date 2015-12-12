// C++ implementation of a fast algorithm for generating samples from a
// discrete distribution.
//
// David Pal, December 2015
//
// To compile the program run:
//
//   g++ -Wall -Wextra -Werror -std=c++11 discrete-distribution.cc

#include <cmath>
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

using std::cout;
using std::endl;

namespace {
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
}

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

      // Fix index.  TODO: This probably not necessary?
      if (index >= buckets_.size()) index = buckets_.size() - 1;

      const Bucket& bucket = buckets_[index];
      if (number < std::get<2>(bucket))
        return std::get<0>(bucket);
      else
        return std::get<1>(bucket);
    }

    result_type min() const {
      return static_cast<result_type>(0);
    }

    result_type max() const {
      return probabilities_.empty()
             ? static_cast<result_type>(0)
             : static_cast<result_type>(probabilities_.size() - 1);
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
    // TODO: Figure out how to replace size_t in Segment with result_type.
    // GCC 4.8.4 refuses to compile it.
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
      const size_t N = probabilities_.size();
      if (N <= 0) {
        buckets_.emplace_back(0, 0, 0.0);
        return;
      }

      // Two stacks in one vector.  First stack grows from the begining of the
      // vector. The second stack grows from the end of the vector.
      std::vector<Segment> segments(N);
      stack_view<Segment, std::vector<Segment>::iterator>
        small(segments.begin());
      stack_view<Segment, std::vector<Segment>::reverse_iterator>
        large(segments.rbegin());

      // Split probabilities into small and large
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
      while (!small.empty() && !large.empty()) {
        const Segment s = small.pop();
        const Segment l = large.pop();

        // Create a mixed bucket
        buckets_.emplace_back(s.second, l.second,
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
        // The last argument is irrelevant as long it's not a NaN.
        buckets_.emplace_back(l.second, l.second, 0.0);
      }

      // This loop can be executed only due to numerical inaccuracies.
      // TODO: Find an example when it actually happens.
      while (!small.empty()) {
        const Segment s = small.pop();
        cout << "Here" << endl;
        // The last argument is irrelevant as long it's not a NaN.
        buckets_.emplace_back(s.second, s.second, 0.0);
      }
    }

    // Uniform distribution over interval [0,1].
    std::uniform_real_distribution<double> uniform_distribution_;

    // List of probabilities
    std::vector<double> probabilities_;
    std::vector<Bucket> buckets_;
};

void Test(const std::vector<double>& weights, const size_t num_samples) {
  std::default_random_engine generator;
  fast_discrete_distribution<int> distribution(weights);
  distribution.PrintBuckets();

  std::vector<size_t> counts(weights.size(), 0);
  for (size_t i = 0; i < num_samples; ++i) {
    const int number = distribution(generator);
    assert(number >= 0);
    assert(number < static_cast<int>(weights.size()));
    ++counts[number];
  }

  std::cout << "counts:" << std::endl;
  for (size_t i = 0; i < weights.size(); ++i)
    cout << i << " (" << weights[i] << ") : "
         << std::string(counts[i], '*') << endl;

  cout << endl;
}

void TestEmpty(const size_t num_samples) {
  std::default_random_engine generator;
  fast_discrete_distribution<int> distribution({});
  distribution.PrintBuckets();

  for (size_t i = 0; i < num_samples; ++i) {
    const int number = distribution(generator);
    assert(number == 0);
  }
}

int main() {
  TestEmpty(100);
  Test({0}, 100);
  Test({1}, 100);
  Test({1, 1}, 200);
  Test({1, 1, 1}, 300);
  Test({1, 1, 2}, 300);
  Test({1, 0, 2}, 300);
  Test({20, 10, 30}, 300);
  Test({0, 1e-20, 0}, 100);
  Test({1 - 1e-10, 1 - 1e-10, 1 - 1e-10}, 100);

  std::discrete_distribution<int> distribution({10.0, 20.0, 30.0});
  cout << distribution << endl;
  return 0;
}
