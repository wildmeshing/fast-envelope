#include "csv_reader.h"

#include <array>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fastEnvelope {
std::vector<std::array<Vector3, 3>> read_CSV_triangle(
    const std::string& input_filename,
    std::vector<int>& in_envelope)
{
    std::vector<std::array<Vector3, 3>> triangle;

    std::ifstream input(input_filename);
    if (!input) {
        throw std::runtime_error("Failed to open CSV file: " + input_filename);
    }

    std::size_t line_number = 0;
    std::string line;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line.front() == '#') continue;

        std::array<double, 10> record {};
        std::istringstream row(line);
        std::string field;
        std::size_t column = 0;
        while (std::getline(row, field, ',')) {
            if (column == record.size()) {
                throw std::runtime_error(
                    "Too many fields in " + input_filename + " at line "
                    + std::to_string(line_number));
            }

            std::size_t parsed_characters = 0;
            try {
                record[column] = std::stod(field, &parsed_characters);
            } catch (const std::exception&) {
                throw std::runtime_error(
                    "Invalid number in " + input_filename + " at line "
                    + std::to_string(line_number));
            }
            if (field.find_first_not_of(" \t\r\n", parsed_characters) != std::string::npos) {
                throw std::runtime_error(
                    "Invalid number in " + input_filename + " at line "
                    + std::to_string(line_number));
            }
            ++column;
        }

        if (column != record.size()) {
            throw std::runtime_error(
                "Expected 10 fields in " + input_filename + " at line "
                + std::to_string(line_number));
        }

        triangle.push_back(
            {{Vector3(record[0], record[1], record[2]),
              Vector3(record[3], record[4], record[5]),
              Vector3(record[6], record[7], record[8])}});
        in_envelope.push_back(static_cast<int>(record[9]));
    }

    if (!input.eof()) {
        throw std::runtime_error("Failed while reading CSV file: " + input_filename);
    }

    return triangle;
}

} // namespace fastEnvelope
