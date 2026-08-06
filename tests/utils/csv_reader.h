#pragma once

#include <fastenvelope/Types.hpp>

#include <array>
#include <string>
#include <vector>

namespace fastEnvelope {

std::vector<std::array<Vector3, 3>> read_CSV_triangle(
    const std::string& input_filename,
    std::vector<int>& in_envelope);

} // namespace fastEnvelope
