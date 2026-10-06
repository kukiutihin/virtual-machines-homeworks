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
static const double JUMP_BOUND = 1.1;

static const size_t MIN_HIGHER_STRIDE = 16;
static const size_t MIN_LOWER_STRIDE = 8;
static const size_t MAX_STRIDE_CL = 128;
static const size_t MIN_SPOTS_CL = 8;
static const size_t MAX_SPOTS_CL = 4 * 1024;
static const size_t MAX_LINE_SIZE = 256;

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
  const size_t elem_size = sizeof(size_t);
  const size_t elem_stride = stride / elem_size;
  size_t *array = new size_t[elem_stride * spots + PAGE_SIZE / elem_size];

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
    array[points[i]] = points[i + 1];
  }
  array[points[points.size() - 1]] = points[0];

  size_t current = points[0];

  for (size_t trial = 0; trial < WARMUP; trial++) {
    current = array[current];
  }

  auto start = steady_clock::now();
  for (size_t trial = 0; trial < TRIALS; trial++) {
    current = array[current];
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

  std::sort(times.begin(), times.end());
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
                                      const Measurements &measurements) {
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
                     double target) {
  while (h_spots - l_spots > 1) {
    size_t spots = (l_spots + h_spots) / 2;
    double time = measure(stride, spots);

    if (time > target * JUMP_BOUND)
      h_spots = spots;
    else
      l_spots = spots;
  }
  return l_spots;
}

std::pair<size_t, std::string>
detect_line_size(size_t cache_capacity, const Measurements &measurements,
                 size_t assoc) {
  std::stringstream buffer;

  size_t line_size = 0;
  double best = 100;

  for (size_t h_stride = MIN_HIGHER_STRIDE; h_stride <= MAX_STRIDE_CL;
       h_stride *= 2) {
    size_t a = cache_capacity / h_stride;
    double target = get_measurement(h_stride, assoc, measurements).value_or(0);

    buffer << BOLD
           << std::format("H_STRIDE = {:<6} Expected spots (A) = {:<6} "
                          "Target time = {:.3g} ns",
                          h_stride, a, target)
           << RESET << "\n";
    buffer << std::string(62, '-') << "\n";
    buffer << std::format("  {:<10} {:<14} {:<12} {:<12} {:<8}\n", "L_stride",
                          "Total_stride", "Spots", "Ratio", "Trend");
    buffer << std::string(62, '-') << "\n";

    std::vector<size_t> results;
    for (size_t l_stride = MIN_LOWER_STRIDE; l_stride < h_stride;
         l_stride *= 2) {
      size_t total_stride = h_stride + l_stride;
      size_t s =
          binary_search(MIN_SPOTS_CL, MAX_SPOTS_CL, total_stride, target);
      results.push_back(s);

      double ratio = static_cast<double>(s) / static_cast<double>(a);
      std::string_view trend_str = (s < a) ? "DEC" : "INC";
      std::string_view color = (s < a) ? BLUE : RED;

      buffer << std::format("  {:<10} {:<14} {:<12} {:<12.4f} {}{}{}\n",
                            l_stride, total_stride, s, ratio, color, trend_str,
                            RESET);
    }
    buffer << "\n";

    std::vector<double> rations;
    for (const auto &x : results)
      rations.push_back((double)x / (double)a);

    std::sort(rations.begin(), rations.end());
    double median = rations[results.size() / 2];

    if (std::abs(1 - median) < best) {
      line_size = h_stride;
      best = std::abs(1 - median);
    }
  }

  return {line_size, buffer.str()};
}

int main() {
  size_t stride = MIN_HIGHER_STRIDE;
  size_t spots = 1;

  for (size_t i = 0; stride <= MAX_STRIDE; stride *= 2, spots = 1, i++) {
    for (size_t j = 0; spots <= MAX_ASSOC; spots++, j++) {
      double time = measure(stride, spots);
      DUMMY = time;
    }
  }

  Measurements measurements;
  stride = MIN_HIGHER_STRIDE;
  spots = 1;

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

  const auto &[line_size, log] =
      detect_line_size(l1_assoc * l1_stride, measurements, l1_assoc);

  std::cout << "\n" << log << "\n" << std::endl;

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
