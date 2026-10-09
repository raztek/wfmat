// Result and Error types [API-03].
#pragma once

#include <cstddef>
#include <optional>
#include <string>

#include <tl/expected.hpp>

namespace wfmat {

// [API-03] unsupported: valid input that needs a later milestone (arcs before M3, clusters before M5).
enum class ErrorCode { invalid_input, numerical_failure, invariant_violation, unsupported };

struct Error {
    ErrorCode code = ErrorCode::invalid_input;
    std::string requirement;               // the spec requirement that failed, e.g. "IN-04"
    std::string message;
    std::optional<std::size_t> loop;       // 0 is the outer loop
    std::optional<std::size_t> segment;    // input segment index (segment i runs from vertex i to i + 1)
    std::string dump_path;                 // diagnostic dump, when one was written
};

template <class T>
using Result = tl::expected<T, Error>;

inline tl::unexpected<Error> make_error(ErrorCode code, std::string requirement, std::string message,
                                        std::optional<std::size_t> loop = std::nullopt,
                                        std::optional<std::size_t> segment = std::nullopt)
{
    return tl::unexpected<Error>(Error{code, std::move(requirement), std::move(message), loop, segment, {}});
}

// One-line description: "IN-04: loop 0, segment 3: ...".
std::string to_string(const Error& e);

} // namespace wfmat
