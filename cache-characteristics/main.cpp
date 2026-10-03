#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <random>
#include <ratio>
#include <sstream>
#include <utility>
#include <vector>

using std::size_t;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::chrono::steady_clock;
using params_t = std::pair<size_t, size_t>;

static const size_t TRIALS = 1000000;
static const size_t WARMUP = 1000000;
static const size_t MEASURES = 5;
static const size_t MAX_ASSOC = 20;
static const size_t MAX_STRIDE = 256 * 1024;
static const size_t MIN_HIGHER_STRIDE = 16;
static const size_t MIN_LOWER_STRIDE = 2;
static const size_t PAGE_SIZE = 1024 * 4;
static const double JUMP_BOUND = 1.25;
static const double CONF_THRESHOLD = 0.8;
static const size_t SPOTS_CL_FIND = 4 * 1024;

static const std::string_view BOLD = "\x1b[1m";
static const std::string_view RESET = "\x1b[0m";
static const std::string_view RED = "\x1b[31m";
static const std::string_view BLUE = "\x1b[34m";
static const size_t COL0_WIDTH = 14;
static const size_t NUM_WIDTH = 10;
static const size_t RESULT_TABLE_COL_WIDHT = 16;

struct Measurement {
  bool is_jump;
  double time;
  size_t stride;
  size_t spots;
};

class Table {
  std::vector<std::vector<Measurement>> measurements;

public:
  Table() : measurements(MAX_ASSOC) {}

  void push(size_t stride, size_t spots, double time) {
    measurements[spots - 1].push_back(Measurement(false, time, stride, spots));
  }

  std::optional<double> get_measurement(size_t spots, size_t stride) const {
    if (spots > measurements.size())
      return {};

    for (size_t i = 0; i < measurements[spots - 1].size(); i++) {
      const auto &measurement = measurements[spots - 1][i];
      if (measurement.stride == stride)
        return measurement.time;
    }

    return {};
  }

  std::pair<std::vector<params_t>, double> detect_jumps() {
    std::vector<double> sample;
    for (size_t s = 1; s < measurements.size() / 3; s++)
      for (size_t m = 0; m < measurements[s].size() / 3; m++)
        sample.push_back(measurements[s][m].time);
    std::sort(sample.begin(), sample.end());

    double l1_median = 0;
    if (!sample.empty())
      l1_median = sample[sample.size() / 2];

    std::vector<params_t> jumps;
    std::vector<size_t> was_jump_on_stride(measurements[0].size());

    for (size_t s = 2; s < measurements.size(); s++) {
      for (size_t m = 0; m < measurements[s].size(); m++) {
        if (was_jump_on_stride[m])
          continue;

        if (measurements[s][m].time > l1_median * JUMP_BOUND) {
          jumps.push_back(
              {measurements[s - 1][m].stride, measurements[s - 1][m].spots});

          measurements[s - 1][m].is_jump = true;
          was_jump_on_stride[m] = true;
        }
      }
    }

    std::sort(jumps.begin(), jumps.end(),
              [](auto &j1, auto &j2) { return j1.first < j2.first; });

    return {jumps, l1_median};
  }

  std::string compile() const {
    std::stringstream buffer;

    buffer << BOLD << std::format("{:<{}}", "spots/stride", COL0_WIDTH);
    for (size_t stride = MIN_HIGHER_STRIDE; stride <= MAX_STRIDE; stride *= 2)
      buffer << std::format("{:>{}}", stride, NUM_WIDTH);
    buffer << RESET << "\n";

    for (size_t spots = 0; spots < measurements.size(); spots++) {
      buffer << BOLD << std::format("{:<{}}", spots + 1, COL0_WIDTH) << RESET;
      for (const auto &measurement : measurements[spots]) {
        auto formatted_num = std::format("{:.3g}", measurement.time);
        if (measurement.is_jump)
          buffer << BOLD << RED;
        buffer << std::format("{:>{}}", formatted_num, NUM_WIDTH) << RESET;
      }
      buffer << "\n";
    }

    return buffer.str();
  }
};

size_t DUMMY;

double measure_once(size_t stride, size_t spots) {
  /* stride = 32, spots = 4;
   * sizeof(size_t) = 8 => stride = 4 elements
   * [*, _, _, _, _, *, _, _, _, _, *, _, _, _, _, *]: *size_t
   */
  const size_t elem_size = sizeof(params_t);
  const size_t elem_stride = stride / elem_size;
  params_t *array = new params_t[elem_stride * spots + PAGE_SIZE / elem_size];

  size_t start_el = 0;
  while ((uintptr_t)(array + start_el) % PAGE_SIZE != 0)
    start_el++;

  std::vector<size_t> points;
  for (size_t i = 0; i < spots; i++)
    points.push_back(start_el + i * elem_stride);

  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(points.begin(), points.end(), g);

  for (size_t i = 0; i < points.size() - 1; i++) {
    array[points[i]] = {points[i + 1], 0};
  }
  array[points[points.size() - 1]] = {points[0], 0};

  size_t current = points[0];

  for (size_t trial = 0; trial < WARMUP; trial++) {
    array[current].second++;
    current = array[current].first;
  }

  auto start = steady_clock::now();
  for (size_t trial = 0; trial < TRIALS; trial++) {
    array[current].second++;
    current = array[current].first;
  }
  auto end = steady_clock::now();

  delete[] array;
  DUMMY = current;
  return duration<double, std::nano>(end - start).count() / TRIALS;
}

double measure(size_t stride, size_t spots) {
  std::vector<double> times;
  for (size_t i = 0; i < MEASURES; i++)
    times.push_back(measure_once(stride, spots));

  return times[times.size() / 2];
}

std::optional<params_t> detect_L1(const std::vector<params_t> &jumps) {
  for (size_t i = 0; i < jumps.size() - 1; i++) {
    const auto &[stride, spots] = jumps[i];
    const auto &[next_stride, next_spots] = jumps[i + 1];

    bool valid = false;

    size_t spots_norm = spots;
    size_t next_spots_norm = next_spots;

    if (stride * 2 == next_stride && spots_norm == next_spots_norm * 2) {
      valid = true;

      // first jump is skipped (its a noise)
      if (i != 0) {
        for (size_t j = i - 1; j > 0; j--)
          if (jumps[j].second < spots)
            valid = false;
      }
    }
    if (valid)
      return {{next_stride, next_spots_norm}};
  }
  return {};
}

size_t binary_search(size_t l_spots, size_t h_spots, size_t stride,
                     double median) {
  while (l_spots < h_spots) {
    size_t spots = (l_spots + h_spots) / 2;
    double time = measure(stride, spots);

    if (time - median > 0.1 * median)
      h_spots = spots;
    else
      l_spots = spots;
  }
  return l_spots;
}

int main() {
  size_t stride = MIN_HIGHER_STRIDE;
  size_t spots = 1;

  Table table;

  for (size_t i = 0; stride <= MAX_STRIDE; stride *= 2, spots = 1, i++) {
    for (size_t j = 0; spots <= MAX_ASSOC; spots++, j++) {
      double time = measure(stride, spots);
      table.push(stride, spots, time);
    }
  }

  const auto [jumps, median] = table.detect_jumps();
  std::cout << "\n" << table.compile() << std::endl;

  const auto &[l1_stride, l1_assoc] =
      detect_L1(jumps).value_or(std::pair(0, 0));

  bool last_was_dec = false;
  bool waiting_inc = false;
  size_t line_size = 0;
  size_t waiting_stride = 0;

  for (size_t stride = 16; stride <= MAX_STRIDE; stride *= 2) {
    size_t spots1 = 0;
    size_t spots2 = 0;

    for (size_t spots = 2; spots <= SPOTS_CL_FIND; spots *= 2) {
      double time = measure(stride, spots);
      double min_time = table.get_measurement(l1_assoc, stride).value_or(0.0);

      if (time - min_time > 0.1 * min_time) {
        spots1 = binary_search(spots / 2, spots, stride, min_time);
        break;
      }
    }

    if (spots1 > 0)
      for (size_t spots = 2;
           spots <= std::max<size_t>(SPOTS_CL_FIND, spots1 * 2); spots *= 2) {
        double time = measure(stride + stride / 2, spots);
        double min_time = table.get_measurement(l1_assoc, stride).value_or(0.0);

        if (time - min_time > 0.1 * min_time) {
          spots2 =
              binary_search(spots / 2, spots, stride + stride / 2, min_time);
          break;
        }
      }

    if (spots1 > 0 && spots2 > 0) {
      if (waiting_inc && spots2 > spots1) {
        line_size = waiting_stride;
        break;
      }

      waiting_inc = false;

      if (last_was_dec && spots2 > spots1) {
        waiting_inc = true;
        waiting_stride = stride;
      }

      if (spots1 > spots2)
        last_was_dec = true;
    }
  }

  std::stringstream buffer;
  buffer << BOLD << std::format("{:<{}}", "Capacity", RESULT_TABLE_COL_WIDHT)
         << std::format("{:<{}}", "Associativity", RESULT_TABLE_COL_WIDHT)
         << std::format("{:<{}}", "Line size", RESULT_TABLE_COL_WIDHT) << RESET
         << "\n";

  buffer << std::format("{:<{}}", l1_stride * l1_assoc, RESULT_TABLE_COL_WIDHT)
         << std::format("{:<{}}", l1_assoc, RESULT_TABLE_COL_WIDHT)
         << std::format("{:<{}}", line_size, RESULT_TABLE_COL_WIDHT);

  std::cout << buffer.str() << "\n" << std::endl;
}
