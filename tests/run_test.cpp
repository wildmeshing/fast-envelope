#include <fastenvelope/FastEnvelope.h>
#include <fastenvelope/Types.hpp>

#include "utils/csv_reader.h"

#include <catch2/catch_all.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace fastEnvelope;

void read_off(
    const std::string& filename,
    std::vector<Vector3>& vertices,
    std::vector<Vector3i>& faces)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Failed to open OFF file: " + filename);
    }

    std::string magic;
    size_t num_vertices = 0;
    size_t num_faces = 0;
    size_t num_edges = 0;
    input >> magic >> num_vertices >> num_faces >> num_edges;
    if (!input || magic != "OFF") {
        throw std::runtime_error("Invalid OFF header: " + filename);
    }

    vertices.resize(num_vertices);
    for (auto& vertex : vertices) {
        input >> vertex[0] >> vertex[1] >> vertex[2];
    }

    faces.resize(num_faces);
    for (auto& face : faces) {
        int face_size = 0;
        input >> face_size >> face[0] >> face[1] >> face[2];
        if (face_size != 3) {
            throw std::runtime_error("Only triangular OFF faces are supported");
        }
    }

    if (!input) {
        throw std::runtime_error("Invalid OFF body: " + filename);
    }
}

std::vector<bool> read_reference_results(const std::string& filename)
{
    std::ifstream input(filename);
    std::string header;
    input >> header;
    if (!input || header != "results") {
        throw std::runtime_error("Invalid reference result file: " + filename);
    }

    std::vector<bool> results;
    bool value = false;
    while (input >> value) {
        results.push_back(value);
    }
    return results;
}

std::vector<bool> run_envelope_queries(
    const std::string& query_filename,
    const std::string& model_filename,
    Scalar relative_epsilon)
{
    std::vector<int> input_labels;
    const auto triangles = read_CSV_triangle(query_filename, input_labels);

    std::vector<Vector3> env_vertices;
    std::vector<Vector3i> env_faces;
    read_off(model_filename, env_vertices, env_faces);

    Vector3 min, max;
    algorithms::get_bb_corners(env_vertices, min, max);
    const Scalar epsilon = (max - min).norm() * relative_epsilon;

    const FastEnvelope fast_envelope(env_vertices, env_faces, epsilon);

    std::vector<bool> results(triangles.size());
    for (size_t i = 0; i < triangles.size(); ++i) {
        results[i] = fast_envelope.is_outside(triangles[i]);
    }
    return results;
}

TEST_CASE("matches main reference output", "[integration][regression]")
{
    const std::string data_path = ENVELOPE_TEST_DATA_DIR;
    const auto actual = run_envelope_queries(
        data_path + "63465.stl_envelope_log.csv",
        data_path + "63465.off",
        1e-3);
    const auto expected = read_reference_results(data_path + "63465_result.csv");

    REQUIRE(actual.size() == 100000);
    REQUIRE(actual.size() == expected.size());

    size_t mismatch_count = 0;
    size_t first_mismatch = actual.size();
    for (size_t i = 0; i < actual.size(); ++i) {
        if (actual[i] != expected[i]) {
            if (first_mismatch == actual.size()) first_mismatch = i;
            ++mismatch_count;
        }
    }
    CAPTURE(first_mismatch);
    REQUIRE(mismatch_count == 0);
    REQUIRE(std::count(actual.begin(), actual.end(), false) == 97688);
}
