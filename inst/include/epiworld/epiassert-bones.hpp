#ifndef EPIWORLD_EPIASSERT_HPP
#define EPIWORLD_EPIASSERT_HPP

/**
 * @defgroup epiassert Argument Assertion Utilities
 * @brief Runtime validation of function arguments.
 *
 * The `EpiAssert` class provides static template methods for validating
 * function arguments. Each method works with both scalar and iterable
 * (container) types via `if constexpr` dispatch. On failure the methods
 * throw `std::range_error` or `std::invalid_argument` with a descriptive
 * message that includes the variable name, the offending value, and an
 * optional caller location string.
 *
 * @par Example
 * @code
 * EpiAssert::check_bounds(pop_size, 0, 10000, "pop_size", "Model::init");
 * EpiAssert::check_probability(rate, "rate", "ModelSIR::ModelSIR");
 * EpiAssert::check_sum(probs, 1.0, "probabilities", "create_init_function");
 * @endcode
 *
 * @note Requires C++17 or later (`if constexpr`, `std::void_t`).
 * @{
 */

namespace epiassert_detail {

    /**
     * @brief Type trait that detects iterable types (types with `begin()`
     *        and `end()`).
     *
     * `std::string` is explicitly excluded so that it is treated as a
     * scalar value rather than a character container.
     */
    template<typename T, typename = void>
    struct is_iterable : std::false_type {};

    template<typename T>
    struct is_iterable<T, std::void_t<
        decltype(std::begin(std::declval<const T &>())),
        decltype(std::end(std::declval<const T &>()))
    >> : std::true_type {};

    template<>
    struct is_iterable<std::string> : std::false_type {};

    /**
     * @brief Convert a numeric value to its string representation.
     *
     * Uses `std::ostringstream` so that floating-point values are printed
     * with enough precision to be useful in error messages.
     */
    template<typename T>
    inline std::string to_str(const T & v)
    {
        std::ostringstream oss;
        oss << v;
        return oss.str();
    }

} // namespace epiassert_detail

/**
 * @brief Argument-assertion utilities for epiworld.
 *
 * All methods are `static` and accept an optional @p varname (for the
 * variable being checked) and @p caller (for the function or method that
 * called the check). Both are used only in the error message.
 */
class EpiAssert {
private:
    /** @brief Format the caller location for error messages. */
    static std::string fmt_location(const std::string & caller)
    {
        if (caller.empty())
            return "";
        return " (in '" + caller + "')";
    }

public:
    // -----------------------------------------------------------------
    //  check_bounds – value(s) in [lower, upper]
    // -----------------------------------------------------------------

    /**
     * @brief Assert that @p value is in `[lower, upper]`.
     *
     * For container types every element is checked individually.
     *
     * @throws std::range_error if any value violates the bounds.
     */
    template<typename T, typename BoundT>
    static void check_bounds(
        const T      value,
        const BoundT lower,
        const BoundT upper,
        const std::string varname = "value",
        const std::string caller  = ""
    )
    {
        if (lower > upper)
        {
            throw std::invalid_argument(
                "check_bounds: 'lower' (" +
                epiassert_detail::to_str(lower) +
                ") must be <= 'upper' (" +
                epiassert_detail::to_str(upper) + ")" +
                fmt_location(caller)
            );
        }

        if constexpr (epiassert_detail::is_iterable<T>::value)
        {
            size_t idx = 0;
            for (const auto & v : value)
            {
                if ((v < lower) || (v > upper))
                {
                    throw std::range_error(
                        "'" + varname + "[" +
                        std::to_string(idx) + "]' must be in [" +
                        epiassert_detail::to_str(lower) + ", " +
                        epiassert_detail::to_str(upper) + "], but got " +
                        epiassert_detail::to_str(v) +
                        fmt_location(caller)
                    );
                }
                ++idx;
            }
        }
        else
        {
            if ((value < lower) || (value > upper))
            {
                throw std::range_error(
                    "'" + varname + "' must be in [" +
                    epiassert_detail::to_str(lower) + ", " +
                    epiassert_detail::to_str(upper) + "], but got " +
                    epiassert_detail::to_str(value) +
                    fmt_location(caller)
                );
            }
        }
    }

    // -----------------------------------------------------------------
    //  check_non_negative – value(s) >= 0
    // -----------------------------------------------------------------

    /**
     * @brief Assert that @p value is non-negative.
     *
     * For container types every element is checked individually.
     *
     * @throws std::range_error if any value is negative.
     */
    template<typename T>
    static void check_non_negative(
        const T value,
        const std::string varname = "value",
        const std::string caller  = ""
    )
    {
        if constexpr (epiassert_detail::is_iterable<T>::value)
        {
            size_t idx = 0;
            for (const auto & v : value)
            {
                if (v < 0)
                {
                    throw std::range_error(
                        "'" + varname + "[" +
                        std::to_string(idx) + "]' must be non-negative"
                        ", but got " + epiassert_detail::to_str(v) +
                        fmt_location(caller)
                    );
                }
                ++idx;
            }
        }
        else
        {
            if (value < 0)
            {
                throw std::range_error(
                    "'" + varname + "' must be non-negative"
                    ", but got " + epiassert_detail::to_str(value) +
                    fmt_location(caller)
                );
            }
        }
    }

    // -----------------------------------------------------------------
    //  check_probability – value(s) in [0, 1]
    // -----------------------------------------------------------------

    /**
     * @brief Assert that @p value is a valid probability in `[0, 1]`.
     *
     * For container types every element is checked individually.
     *
     * @throws std::range_error if any value is outside `[0, 1]`.
     */
    template<typename T>
    static void check_probability(
        const T value,
        const std::string varname = "value",
        const std::string caller  = ""
    )
    {
        if constexpr (epiassert_detail::is_iterable<T>::value)
        {
            size_t idx = 0;
            for (const auto & v : value)
            {
                if (v < 0.0 || v > 1.0)
                {
                    throw std::range_error(
                        "'" + varname + "[" +
                        std::to_string(idx) +
                        "]' must be a probability in [0, 1]"
                        ", but got " + epiassert_detail::to_str(v) +
                        fmt_location(caller)
                    );
                }
                ++idx;
            }
        }
        else
        {
            if (value < 0.0 || value > 1.0)
            {
                throw std::range_error(
                    "'" + varname +
                    "' must be a probability in [0, 1]"
                    ", but got " + epiassert_detail::to_str(value) +
                    fmt_location(caller)
                );
            }
        }
    }

    // -----------------------------------------------------------------
    //  check_sum – container elements sum to target ± tolerance
    // -----------------------------------------------------------------

    /**
     * @brief Assert that the elements of @p values sum to @p target
     *        within ± @p tolerance.
     *
     * @throws std::invalid_argument if the sum is outside
     *         `[target − tolerance, target + tolerance]`.
     */
    template<typename T>
    static void check_sum(
        const T      values,
        double         target,
        const std::string varname   = "values",
        const std::string caller    = "",
        double         tolerance = 1e-8
    )
    {
        static_assert(
            epiassert_detail::is_iterable<T>::value,
            "check_sum requires an iterable type."
        );

        if (tolerance < 0.0)
        {
            throw std::invalid_argument(
                "check_sum: 'tolerance' must be non-negative, but got " +
                epiassert_detail::to_str(tolerance) +
                fmt_location(caller)
            );
        }

        double s = 0.0;
        for (const auto & v : values)
            s += static_cast<double>(v);

        if (std::abs(s - target) > tolerance)
        {
            throw std::invalid_argument(
                "'" + varname + "' elements must sum to " +
                epiassert_detail::to_str(target) +
                " (tolerance " + epiassert_detail::to_str(tolerance) +
                "), but got " + epiassert_detail::to_str(s) +
                fmt_location(caller)
            );
        }
    }

    // -----------------------------------------------------------------
    //  check_size – container has expected number of elements
    // -----------------------------------------------------------------

    /**
     * @brief Assert that @p values has exactly @p expected elements.
     *
     * @throws std::invalid_argument if the size differs.
     */
    template<typename T>
    static void check_size(
        const T values,
        size_t    expected,
        const std::string varname = "values",
        const std::string caller  = ""
    )
    {
        static_assert(
            epiassert_detail::is_iterable<T>::value,
            "check_size requires an iterable type."
        );

        const auto actual_size = static_cast<size_t>(
            std::distance(std::begin(values), std::end(values))
        );

        if (actual_size != expected)
        {
            throw std::invalid_argument(
                "'" + varname + "' must have " +
                std::to_string(expected) + " elements, but got " +
                std::to_string(actual_size) +
                fmt_location(caller)
            );
        }
    }

    // -----------------------------------------------------------------
    //  check – custom predicate validation
    // -----------------------------------------------------------------

    /**
     * @brief Assert that @p pred returns `true` for @p value.
     *
     * The predicate receives the value by `const` reference and should
     * return `true` when the value is valid.
     *
     * @throws std::invalid_argument with the caller-provided @p message
     *         when the predicate returns `false`.
     */
    template<typename T, typename Predicate>
    static void check(
        const T     value,
        Predicate     pred,
        const std::string message,
        const std::string caller = ""
    )
    {
        if (!pred(value))
        {
            throw std::invalid_argument(
                message + fmt_location(caller)
            );
        }
    }
};

/** @} */  // end of epiassert group

#endif
