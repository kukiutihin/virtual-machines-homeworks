#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <ratio>
#include <set>
#include <utility>
#include <vector>

using std::size_t;
using std::chrono::duration;
using std::chrono::high_resolution_clock;
using std::chrono::steady_clock;
using params_t = std::pair<size_t, size_t>;
using jumps_t = std::map<size_t, std::set<size_t>>;

static const size_t TRIALS = 1000000;
static const size_t WARMUP = 1000000;
static const size_t MEASURES = 5;
static const size_t PAGE_SIZE = 1024 * 4;

static const size_t MAX_ASSOC = 25;
static const size_t MAX_STRIDE = 256 * 1024;
static const double JUMP_BOUND = 1.25;

static const size_t MIN_HIGHER_STRIDE = 16;
static const size_t MIN_LOWER_STRIDE = 2;
static const double CONF_THRESHOLD = 0.8;
static const size_t SPOTS_CL_FIND = 4 * 1024;

static const std::string_view BOLD = "\x1b[1m";
static const std::string_view RESET = "\x1b[0m";
static const std::string_view RED = "\x1b[31m";
static const std::string_view BLUE = "\x1b[34m";

static const size_t COL0_WIDTH = 14;
static const size_t NUM_WIDTH = 10;
static const size_t RESULT_TABLE_COL_WIDHT = 16;

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

struct Measurement {
  bool is_jump;
  double time;
  size_t stride;
  size_t spots;
};

using Measurements = std::map<size_t, std::vector<Measurement>>;

std::optional<double> get_measurement(size_t stride, size_t spots,
                                      const Measurements measurements) {
  if (!measurements.contains(stride))
    return {};

  if (measurements.at(stride).size() < spots)
    return {};

  return measurements.at(stride)[spots - 1].time;
}

jumps_t detect_jumps(const Measurements &measurements) {
  jumps_t jumps;

  size_t prev_spots = 0;
  for (const auto &[stride, ms] : measurements) {
    double prev_time = 0;
    for (const auto &m : ms) {
      if (prev_spots > 0 && prev_time > 0 && m.time > prev_time * JUMP_BOUND)
        jumps[stride].insert(prev_spots);
      prev_time = m.time;
      prev_spots = m.spots;
    }
  }

  return jumps;
}

std::string compile_table(const Measurements &measurements,
                          const jumps_t &jumps) {
  std::stringstream buffer;

  buffer << BOLD << std::format("{:<{}}", "spots/stride", COL0_WIDTH);
  for (const auto &[stride, _] : measurements)
    buffer << std::format("{:>{}}", stride, NUM_WIDTH);
  buffer << RESET << "\n";

  size_t max_spots = 0;
  for (const auto &[_, ms] : measurements)
    for (const auto &m : ms)
      max_spots = std::max(max_spots, m.spots);

  std::vector<std::map<size_t, double>> spots_to_strides(max_spots + 1);

  for (const auto &[stride, ms] : measurements)
    for (const auto &m : ms)
      spots_to_strides[m.spots - 1][stride] = m.time;

  for (size_t spots = 0; spots < max_spots; spots++) {
    buffer << BOLD << std::format("{:<{}}", spots + 1, COL0_WIDTH) << RESET;
    for (const auto &[stride, time] : spots_to_strides[spots]) {
      auto formatted_num = std::format("{:.3g}", time);
      if (jumps.contains(stride)) {
        if (jumps.at(stride).contains(spots + 1))
          buffer << BOLD << RED;
      }
      buffer << std::format("{:>{}}", formatted_num, NUM_WIDTH) << RESET;
    }
    buffer << "\n";
  }

  return buffer.str();
}

std::optional<params_t> detect_L1(const jumps_t &jumps) {
  std::vector<params_t> entities;

  for (const auto &[stride, spots_vec] : jumps)
    for (const auto spots : spots_vec)
      if (spots % 2 == 0 && (spots / 2 % 2) == 0)
        if (jumps.contains(stride * 2))
          if (jumps.at(stride * 2).contains(spots / 2))
            entities.push_back({stride * 2, spots / 2});

  std::sort(entities.begin(), entities.end(),
            [](const auto &ent1, const auto &ent2) {
              return ent1.first * ent1.second < ent2.first * ent2.second;
            });

  return entities.size() > 0 ? std::optional(entities[0]) : std::nullopt;
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

  Measurements measurements;

  for (size_t i = 0; stride <= MAX_STRIDE; stride *= 2, spots = 1, i++) {
    for (size_t j = 0; spots <= MAX_ASSOC; spots++, j++) {
      double time = measure(stride, spots);
      measurements[stride].push_back(Measurement(false, time, stride, spots));
    }
  }

  jumps_t jumps = detect_jumps(measurements);
  std::cout << "\n" << compile_table(measurements, jumps) << std::endl;

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
      double min_time =
          get_measurement(stride, l1_assoc, measurements).value_or(0.0);

      if (time - min_time > 0.1 * min_time) {
        spots1 = binary_search(spots / 2, spots, stride, min_time);
        break;
      }
    }

    if (spots1 > 0)
      for (size_t spots = 2;
           spots <= std::max<size_t>(SPOTS_CL_FIND, spots1 * 2); spots *= 2) {
        double time = measure(stride + stride / 2, spots);
        double min_time =
            get_measurement(stride, l1_assoc, measurements).value_or(0.0);

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
