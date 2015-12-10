// C++ implementation of a fast algorithm for generating samples from a
// discrete distribution.
//
// David Pal, December 2015
//
// To compile the program run:
//
//   g++ -Wall -Wextra -Werror -std=c++11 discrete-distribution.cc

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

using std::cout;
using std::endl;

// Simple stack with fixed capacity.
template<typename T>
class FixedSizeStack {
  public:
    FixedSizeStack(const size_t capacity) : top(0) {
      data.reserve(capacity);
    }

    void push(const T& element) {
      if (data.size() > top) {
        data[top] = element;
      } else {
        data.push_back(element);
      }
      ++top;
    }

    T& pop() {
      return data[--top];
    }

    size_t size() {
      return top;
    }

  private:
    size_t top;
    std::vector<T> data;
};


class DiscreteDistribution {
  public:
    DiscreteDistribution(const std::vector<double>& weights)
      : generator(),
        distribution(0.0, 1.0) {
      Preprocess(weights);
    }

    size_t GetSample() {
      const double number = distribution(generator);
      size_t index = floor(buckets.size() * number);

      // Fix index
      if (index >= buckets.size()) index = buckets.size() - 1;

      const Bucket& bucket = buckets[index];
      if (number < std::get<2>(bucket))
        return std::get<0>(bucket);
      else
        return std::get<1>(bucket);
    }

    void PrintBuckets() {
      cout << "buckets.size() = " << buckets.size() << endl;
      for (auto bucket : buckets) {
        cout << std::get<0>(bucket) << "  "
             << std::get<1>(bucket) << "  "
             << std::get<2>(bucket) << "  "
             << endl;
      }
    }

  private:
    typedef std::pair<double, size_t> Segment;
    typedef std::tuple<size_t, size_t, double> Bucket;

    void Preprocess(const std::vector<double>& weights) {
      const double sum = std::accumulate(weights.cbegin(), weights.cend(), 0.0);

      // Reserve capacity to avoid reallocations.
      // This might take twice as much memory as is actually needed.
      const size_t N = weights.size();
      FixedSizeStack<Segment> small(N);
      FixedSizeStack<Segment> big(N);

      // Normalize the probabilities.
      for (size_t i = 0; i < weights.size(); ++i) {
        const double probability = weights[i] / sum;
        if (probability < (1.0 / N))
          small.push(Segment(probability, i));
        else
          big.push(Segment(probability, i));
      }

      buckets.reserve(N);

      size_t i = 0;
      while (small.size() > 0) {
        const Segment s = small.pop();
        const Segment b = big.pop();

        // Create a mixed bucket
        buckets.emplace_back(s.second,
                             b.second,
                             s.first + static_cast<double>(i) / N);

        // Calculate the length of the left-over segment
        const double left_over = s.first + b.first - static_cast<double>(i) / N;

        // Re-insert the left-over segment
        if (left_over < (1.0 / N))
          small.push(Segment(left_over, b.second));
        else
          big.push(Segment(left_over, b.second));

        ++i;
      }

      // Create pure buckets
      while (big.size() > 0) {
        const Segment b = big.pop();
        buckets.emplace_back(b.second, b.second, 0.0);
      }
    }

    std::vector<Bucket> buckets;

    // Random number generator and uniform distribution.
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
};

void Test(const std::vector<double>& weights, const size_t nrolls) {
  DiscreteDistribution d(weights);
  d.PrintBuckets();

  std::vector<size_t> counts(weights.size(), 0);
  for (size_t i = 0; i < nrolls; ++i) {
    const int number = d.GetSample();
    ++counts[number];
  }

  std::cout << "counts:" << std::endl;
  for (size_t i = 0; i < weights.size(); ++i)
    cout << i << " (" << weights[i] << ") : "
         << std::string(counts[i], '*') << endl;

  cout << endl;
}

int main() {
  Test({1, 1, 1}, 300);
  Test({1, 1}, 200);
  Test({1}, 100);
  Test({1, 1, 2}, 300);
  Test({1, 0, 2}, 300);
  return 0;
}
