#!/usr/bin/env python3

import sys


################################################################################
# Helpers
################################################################################


def repeat(value, size):
    """Build a C++ SIMD expression containing size copies of value.

    Value represents one Complex32 [real, imag] pair, so the result contains
    2 * size float lanes.

    Args:
        value: C++ expression for one Complex32 value.
        size: Number of copies to place in the SIMD expression.
    """
    if size == 1:
        return value

    half = size // 2
    repeated = repeat(value, half)
    return f"SIMD<float,{2 * size}>({repeated}, {repeated})"


def add_offset(first, offset):
    """Build the C++ index expression first + offset.

    Dynamic expressions are parenthesized, while a zero offset is omitted.

    Args:
        first: Base C++ index expression.
        offset: Constant offset added to the base expression.
    """
    if offset == 0:
        return str(first)
    return f"({first})+{offset}"


def load_complex32_b_gather(row, first, size):
    """Gather size Complex32 entries into one C++ SIMD expression.

    This is used when adjacent logical columns of B are not contiguous.

    Args:
        row: C++ expression selecting the B row.
        first: C++ expression selecting the first B column.
        size: Number of Complex32 entries to gather.
    """
    if size == 1:
        return f"SIMD<float,2>((float*)(void*)b.Addr({row},{first}))"

    half = size // 2
    return (
        f"SIMD<float,{2 * size}>({load_complex32_b_gather(row, first, half)}, "
        f"{load_complex32_b_gather(row, add_offset(first, half), half)})"
    )


def load_complex32_b(row, first, size, ordering_b):
    """Build a SIMD load expression for one group of B entries.

    Row-major B uses one contiguous vector load; column-major B is gathered.

    Args:
        row: C++ expression selecting the B row.
        first: C++ expression selecting the first B column.
        size: Number of Complex32 entries to load.
        ordering_b: C++ matrix-ordering name for B.
    """
    if ordering_b == "RowMajor":
        return (
            f"SIMD<float,{2 * size}>"
            f"((float*)(void*)b.Addr({row},{first}))"
        )
    return load_complex32_b_gather(row, first, size)


def load_float_b_gather(row, first, size):
    """Gather size float entries from logical columns of B.

    Adjacent columns are strided for column-major B, so the expression is
    assembled recursively from scalar loads.

    Args:
        row: C++ expression selecting the B row.
        first: C++ expression selecting the first B column.
        size: Number of float entries to gather.
    """
    if size == 1:
        return f"SIMD<float,1>(b({row},{first}))"

    half = size // 2
    return (
        f"SIMD<float,{size}>({load_float_b_gather(row, first, half)}, "
        f"{load_float_b_gather(row, add_offset(first, half), half)})"
    )


def load_float_b(row, first, size, ordering_b):
    """Build the B-vector load used by a float outer-product kernel.

    Row-major B is loaded contiguously; column-major B uses the same logical
    SIMD lanes, assembled with scalar gathers.

    Args:
        row: C++ expression selecting the B row.
        first: C++ expression selecting the first B column.
        size: Number of float entries to load.
        ordering_b: C++ matrix-ordering name for B.
    """
    if ordering_b == "RowMajor":
        return f"SIMD<float,{size}>((const float*)b.Addr({row},{first}))"
    return load_float_b_gather(row, first, size)


def bool_name(value):
    """Convert a Python boolean to its lowercase C++ literal.

    The returned string is inserted directly into generated template arguments.

    Args:
        value: Boolean value to format.
    """
    return "true" if value else "false"


################################################################################
# Float kernels
################################################################################


def generate_float_kernel(
    out, simd_size, height, width, add, pos, ordering_a, ordering_b
):
    """Emit one height-by-width float GEMM microkernel specialization.

    The kernel accumulates across K with SIMD lanes following output columns.
    A and B ordering are fixed in the specialization; ADD/POS controls C.

    Args:
        out: Destination text stream.
        simd_size: Native SIMD width measured in doubles.
        height: Number of output rows computed by the kernel.
        width: Number of output columns computed by the kernel.
        add: Whether the result is accumulated into C.
        pos: Whether the product is added with a positive sign.
        ordering_a: C++ matrix-ordering name for A.
        ordering_b: C++ matrix-ordering name for B.
    """
    group_size = min(width, 2 * simd_size)
    groups = width // group_size
    add_name = bool_name(add)
    pos_name = bool_name(pos)

    out.write(
        f"""template <>
INLINE void MatKernelFloat<{height},{width},{add_name},{pos_name},{ordering_a},{ordering_b}>
(size_t nk,
 BareSliceMatrix<float,{ordering_a}> a,
 BareSliceMatrix<float,{ordering_b}> b,
 BareSliceMatrix<float,RowMajor> c)
{{
"""
    )

    for row in range(height):
        for group in range(groups):
            out.write(f"    SIMD<float,{group_size}> sum{row}_{group}(0.0f);\n")

    out.write("    for (size_t k = 0; k < nk; k++)\n      {\n")

    for row in range(height):
        out.write(f"        float a{row} = a({row},k);\n")

    for group in range(groups):
        out.write(
            f"        auto b{group} = "
            f"{load_float_b('k', group * group_size, group_size, ordering_b)};\n"
        )

    for group in range(groups):
        for row in range(height):
            out.write(
                f"        sum{row}_{group} = "
                f"FMA(a{row}, b{group}, sum{row}_{group});\n"
            )

    out.write("      }\n")

    for row in range(height):
        for group in range(groups):
            first = group * group_size
            target = f"c.Addr({row},{first})"
            value = f"sum{row}_{group}" if pos else f"-sum{row}_{group}"
            if add:
                operation = "+" if pos else "-"
                value = (
                    f"SIMD<float,{group_size}>({target}) {operation} "
                    f"sum{row}_{group}"
                )
            out.write(f"    ({value}).Store({target});\n")

    out.write("}\n\n")


def generate_float_inner_product_kernel(
    out, simd_size, height, width, add, pos, ordering_a, ordering_b
):
    """Emit a float inner-product tile for a K-contiguous A/B layout.

    Both operands are loaded contiguously along K. Each SIMD accumulator holds
    a partial dot product for one output entry and is reduced before storing C.

    Args:
        out: Destination text stream.
        simd_size: Native SIMD width measured in doubles.
        height: Number of output rows computed by the kernel.
        width: Number of output columns computed by the kernel.
        add: Whether the result is accumulated into C.
        pos: Whether the product is added with a positive sign.
        ordering_a: C++ matrix-ordering name for A.
        ordering_b: C++ matrix-ordering name for B.
    """
    vector_size = 2 * simd_size
    add_name = bool_name(add)
    pos_name = bool_name(pos)

    out.write(
        f"""template <>
INLINE void MatKernelFloat<{height},{width},{add_name},{pos_name},{ordering_a},{ordering_b}>
(size_t nk,
 BareSliceMatrix<float,{ordering_a}> a,
 BareSliceMatrix<float,{ordering_b}> b,
 BareSliceMatrix<float,RowMajor> c)
{{
"""
    )

    for row in range(height):
        for col in range(width):
            out.write(
                f"    SIMD<float,{vector_size}> sum{row}_{col}(0.0f);\n"
            )

    out.write("    size_t k = 0;\n")
    out.write(
        f"    for ( ; k+{vector_size} <= nk; k += {vector_size})\n"
        "      {\n"
    )

    for row in range(height):
        out.write(
            f"        SIMD<float,{vector_size}> a{row}(a.Addr({row},k));\n"
        )
    for col in range(width):
        out.write(
            f"        SIMD<float,{vector_size}> b{col}(b.Addr(k,{col}));\n"
        )

    for row in range(height):
        for col in range(width):
            out.write(
                f"        sum{row}_{col} = "
                f"FMA(a{row}, b{col}, sum{row}_{col});\n"
            )

    out.write("      }\n\n")

    for row in range(height):
        for col in range(width):
            out.write(
                f"    float value{row}_{col} = HSum(sum{row}_{col});\n"
            )

    out.write("\n    for ( ; k < nk; k++)\n      {\n")
    for row in range(height):
        out.write(f"        float a{row} = a({row},k);\n")
    for col in range(width):
        out.write(f"        float b{col} = b(k,{col});\n")
    for row in range(height):
        for col in range(width):
            out.write(f"        value{row}_{col} += a{row}*b{col};\n")
    out.write("      }\n")

    for row in range(height):
        for col in range(width):
            out.write(
                f"    Func<{add_name},{pos_name}>"
                f"(c({row},{col}), value{row}_{col});\n"
            )

    out.write("}\n\n")


def generate_float_small_k_tile(out, simd_size, width, height):
    """Emit one width-outer float tile loop for a fixed small K.

    B values are loaded once per column panel and reused across several rows.
    The first product initializes each accumulator directly; later products
    use FMA. The generated K operations disappear when AW is smaller.

    Args:
        out: Destination text stream.
        simd_size: Native SIMD width measured in doubles.
        width: Number of output columns in this generated tile loop.
        height: Number of output rows in each full tile.
    """
    group_size = min(width, 2 * simd_size)
    groups = width // group_size

    out.write(f"    for ( ; j+{width} <= bw; j += {width})\n      {{\n")

    for k in range(8):
        for group in range(groups):
            first = add_offset("j", group * group_size)
            direct = (
                f"SIMD<float,{group_size}>"
                f"((const float*)b.Addr({k},{first}))"
            )
            gathered = load_float_b_gather(str(k), first, group_size)
            out.write(f"        SIMD<float,{group_size}> b{k}_{group};\n")
            out.write(f"        if constexpr (AW > {k})\n          {{\n")
            out.write("            if constexpr (OB == RowMajor)\n")
            out.write(f"              b{k}_{group} = {direct};\n")
            out.write("            else\n")
            out.write(f"              b{k}_{group} = {gathered};\n")
            out.write("          }\n")

    out.write("\n        size_t i = 0;\n")
    out.write(f"        for ( ; i+{height} <= ah; i += {height})\n          {{\n")

    for row in range(height):
        for group in range(groups):
            out.write(
                f"            SIMD<float,{group_size}> "
                f"sum{row}_{group};\n"
            )

    out.write("            if constexpr (AW == 0)\n              {\n")
    for row in range(height):
        for group in range(groups):
            out.write(f"                sum{row}_{group} = 0.0f;\n")
    out.write("              }\n")

    out.write("            if constexpr (AW > 0)\n              {\n")
    for row in range(height):
        out.write(f"                float a{row} = a(i+{row},0);\n")
    for group in range(groups):
        for row in range(height):
            out.write(
                f"                sum{row}_{group} = "
                f"a{row} * b0_{group};\n"
            )
    out.write("              }\n")

    for k in range(1, 8):
        out.write(f"            if constexpr (AW > {k})\n              {{\n")
        for row in range(height):
            out.write(f"                float a{row} = a(i+{row},{k});\n")
        for group in range(groups):
            for row in range(height):
                out.write(
                    f"                sum{row}_{group} = "
                    f"FMA(a{row}, b{k}_{group}, sum{row}_{group});\n"
                )
        out.write("              }\n")

    for row in range(height):
        for group in range(groups):
            first = group * group_size
            target = f"c.Addr(i+{row},j+{first})"
            name = f"sum{row}_{group}"
            out.write("            if constexpr (ADD)\n              {\n")
            out.write("                if constexpr (POS)\n")
            out.write(
                f"                  {name} = SIMD<float,{group_size}>"
                f"({target}) + {name};\n"
            )
            out.write("                else\n")
            out.write(
                f"                  {name} = SIMD<float,{group_size}>"
                f"({target}) - {name};\n"
            )
            out.write("              }\n")
            out.write("            else if constexpr (!POS)\n")
            out.write(f"              {name} = -{name};\n")
            out.write(f"            {name}.Store({target});\n")

    out.write("          }\n")
    out.write("\n        for ( ; i < ah; i++)\n          {\n")

    for group in range(groups):
        out.write(f"            SIMD<float,{group_size}> sum{group};\n")

    out.write("            if constexpr (AW == 0)\n              {\n")
    for group in range(groups):
        out.write(f"                sum{group} = 0.0f;\n")
    out.write("              }\n")

    out.write("            if constexpr (AW > 0)\n              {\n")
    out.write("                float av = a(i,0);\n")
    for group in range(groups):
        out.write(f"                sum{group} = av * b0_{group};\n")
    out.write("              }\n")

    for k in range(1, 8):
        out.write(f"            if constexpr (AW > {k})\n              {{\n")
        out.write(f"                float av = a(i,{k});\n")
        for group in range(groups):
            out.write(
                f"                sum{group} = "
                f"FMA(av, b{k}_{group}, sum{group});\n"
            )
        out.write("              }\n")

    for group in range(groups):
        first = group * group_size
        target = f"c.Addr(i,j+{first})"
        name = f"sum{group}"
        out.write("            if constexpr (ADD)\n              {\n")
        out.write("                if constexpr (POS)\n")
        out.write(
            f"                  {name} = SIMD<float,{group_size}>"
            f"({target}) + {name};\n"
        )
        out.write("                else\n")
        out.write(
            f"                  {name} = SIMD<float,{group_size}>"
            f"({target}) - {name};\n"
        )
        out.write("              }\n")
        out.write("            else if constexpr (!POS)\n")
        out.write(f"              {name} = -{name};\n")
        out.write(f"            {name}.Store({target});\n")

    out.write("          }\n")
    out.write("      }\n")


def generate_float_small_k(out, simd_size):
    """Emit the complete fixed-K float multiplication kernel.

    Wider tiles are used when K is at most two and the architecture provides
    32 SIMD registers. Smaller tile loops cover every remaining output column.

    Args:
        out: Destination text stream.
        simd_size: Native SIMD width measured in doubles.
    """
    out.write(
        """template <size_t AW, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
inline void MatKernelFloatSmallK
(size_t ah, size_t bw,
 BareSliceMatrix<float,OA> a,
 BareSliceMatrix<float,OB> b,
 BareSliceMatrix<float,RowMajor> c)
{
  static_assert(AW <= 8);
  size_t j = 0;

#if defined(__AVX512F__) || defined(__arm64__) || defined(__aarch64__)
  if constexpr (AW <= 2)
    {
"""
    )
    generate_float_small_k_tile(out, simd_size, 8 * simd_size, 4)
    out.write("    }\n")

    width = 4 * simd_size
    while width >= 1:
        generate_float_small_k_tile(out, simd_size, width, 4)
        width //= 2

    out.write("#else\n")

    width = 2 * simd_size
    while width >= 1:
        generate_float_small_k_tile(out, simd_size, width, 2)
        width //= 2

    out.write("#endif\n}\n\n")


################################################################################
# Complex32 kernels
################################################################################


def generate_complex32_kernel(
    out, simd_size, height, width, add, pos, ordering_a, ordering_b
):
    """Emit one height-by-width Complex32 GEMM microkernel specialization.

    The A/B orderings and ADD/POS operation are fixed template arguments, while
    complex products are accumulated with FMAComplex before storing C.

    Args:
        out: Destination text stream.
        simd_size: Native SIMD width measured in doubles.
        height: Number of output rows computed by the kernel.
        width: Number of output columns computed by the kernel.
        add: Whether the result is accumulated into C.
        pos: Whether the product is added with a positive sign.
        ordering_a: C++ matrix-ordering name for A.
        ordering_b: C++ matrix-ordering name for B.
    """
    group_size = min(width, simd_size)
    groups = width // group_size
    lanes = 2 * group_size
    add_name = bool_name(add)
    pos_name = bool_name(pos)

    out.write(
        f"""template <>
INLINE void MatKernelComplex32<{height},{width},{add_name},{pos_name},{ordering_a},{ordering_b}>
(size_t nk,
 BareSliceMatrix<Complex32,{ordering_a}> a,
 BareSliceMatrix<Complex32,{ordering_b}> b,
 BareSliceMatrix<Complex32,RowMajor> c)
{{
"""
    )

    for row in range(height):
        for group in range(groups):
            out.write(f"    SIMD<float,{lanes}> sum{row}_{group}(0.0f);\n")

    out.write("    for (size_t k = 0; k < nk; k++)\n      {\n")

    for row in range(height):
        out.write(
            f"        SIMD<float,2> a{row}_((float*)(void*)a.Addr({row},k));\n"
        )
        out.write(f"        auto a{row} = {repeat(f'a{row}_', group_size)};\n")

    for group in range(groups):
        out.write(
            f"        auto b{group} = "
            f"{load_complex32_b('k', group * group_size, group_size, ordering_b)};\n"
        )

    for group in range(groups):
        for row in range(height):
            out.write(f"        FMAComplex(a{row}, b{group}, sum{row}_{group});\n")

    out.write("      }\n")

    for row in range(height):
        for group in range(groups):
            for lane in range(group_size):
                col = group * group_size + lane
                out.write(
                    f"    Func<{add_name},{pos_name}>(c({row},{col}), "
                    f"Complex32(sum{row}_{group}[{2 * lane}], "
                    f"sum{row}_{group}[{2 * lane + 1}]));\n"
                )

    out.write("}\n\n")


def generate_complex32_small_k_tile(out, simd_size, width, height):
    """Emit one width-outer tile loop for a fixed small inner dimension.

    Each B panel is loaded once, reused across height-row tiles, and also used
    for any remaining individual rows.

    Args:
        out: Destination text stream.
        simd_size: Native SIMD width measured in doubles.
        width: Number of output columns in this generated tile loop.
        height: Number of output rows in each full tile.
    """
    group_size = min(width, simd_size)
    groups = width // group_size
    lanes = 2 * group_size

    out.write(f"    for ( ; j+{width} <= bw; j += {width})\n      {{\n")

    for k in range(8):
        for group in range(groups):
            first = add_offset("j", group * group_size)
            direct = (
                f"SIMD<float,{lanes}>"
                f"((float*)(void*)b.Addr({k},{first}))"
            )
            gathered = load_complex32_b_gather(str(k), first, group_size)
            out.write(f"        SIMD<float,{lanes}> b{k}_{group};\n")
            out.write(f"        if constexpr (AW > {k})\n          {{\n")
            out.write("            if constexpr (OB == RowMajor)\n")
            out.write(f"              b{k}_{group} = {direct};\n")
            out.write("            else\n")
            out.write(f"              b{k}_{group} = {gathered};\n")
            out.write("          }\n")

    out.write("\n        size_t i = 0;\n")
    out.write(f"        for ( ; i+{height} <= ah; i += {height})\n          {{\n")

    for row in range(height):
        for group in range(groups):
            out.write(
                f"            SIMD<float,{lanes}> sum{row}_{group}(0.0f);\n"
            )

    for k in range(8):
        out.write(f"            if constexpr (AW > {k})\n              {{\n")
        for row in range(height):
            out.write(
                f"                SIMD<float,2> a{row}_"
                f"((float*)(void*)a.Addr(i+{row},{k}));\n"
            )
            out.write(
                f"                auto a{row} = "
                f"{repeat(f'a{row}_', group_size)};\n"
            )
        for group in range(groups):
            for row in range(height):
                out.write(
                    f"                FMAComplex(a{row}, b{k}_{group}, "
                    f"sum{row}_{group});\n"
                )
        out.write("              }\n")

    for row in range(height):
        for group in range(groups):
            for lane in range(group_size):
                col = group * group_size + lane
                out.write(
                    f"            Func<ADD,POS>(c(i+{row},j+{col}), "
                    f"Complex32(sum{row}_{group}[{2 * lane}], "
                    f"sum{row}_{group}[{2 * lane + 1}]));\n"
                )

    out.write("          }\n")
    out.write("\n        for ( ; i < ah; i++)\n          {\n")

    for group in range(groups):
        out.write(f"            SIMD<float,{lanes}> sum{group}(0.0f);\n")

    for k in range(8):
        out.write(f"            if constexpr (AW > {k})\n              {{\n")
        out.write(
            f"                SIMD<float,2> a_"
            f"((float*)(void*)a.Addr(i,{k}));\n"
        )
        out.write(f"                auto av = {repeat('a_', group_size)};\n")
        for group in range(groups):
            out.write(
                f"                FMAComplex(av, b{k}_{group}, sum{group});\n"
            )
        out.write("              }\n")

    for group in range(groups):
        for lane in range(group_size):
            col = group * group_size + lane
            out.write(
                f"            Func<ADD,POS>(c(i,j+{col}), "
                f"Complex32(sum{group}[{2 * lane}], "
                f"sum{group}[{2 * lane + 1}]));\n"
            )

    out.write("          }\n")
    out.write("      }\n")


def generate_complex32_small_k(out, simd_size):
    """Emit the complete fixed-K Complex32 multiplication kernel.

    It selects wider tiles on 32-register targets and emits progressively
    smaller widths to cover the remaining columns.

    Args:
        out: Destination text stream.
        simd_size: Native SIMD width measured in doubles.
    """
    out.write(
        """template <size_t AW, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
inline void MatKernelComplex32SmallK
(size_t ah, size_t bw,
 BareSliceMatrix<Complex32,OA> a,
 BareSliceMatrix<Complex32,OB> b,
 BareSliceMatrix<Complex32,RowMajor> c)
{
  static_assert(AW <= 8);
  size_t j = 0;

#if defined(__AVX512F__) || defined(__arm64__) || defined(__aarch64__)
  if constexpr (AW <= 2)
    {
"""
    )
    generate_complex32_small_k_tile(out, simd_size, 4 * simd_size, 4)
    out.write("    }\n")

    width = 2 * simd_size
    while width >= 1:
        generate_complex32_small_k_tile(out, simd_size, width, 4)
        width //= 2

    out.write("#else\n")

    width = simd_size
    while width >= 1:
        generate_complex32_small_k_tile(out, simd_size, width, 2)
        width //= 2

    out.write("#endif\n}\n\n")


################################################################################
# Header generation
################################################################################


def generate(output_name, simd_size):
    """Write the complete generated kernel header for one SIMD size.

    The header contains small-K code, all float and Complex32 tile
    specializations, and the architecture-dependent maximum tile widths.

    Args:
        output_name: Path of the generated C++ header.
        simd_size: Native SIMD width measured in doubles.
    """
    with open(output_name, "w", encoding="utf-8", newline="\n") as out:
        out.write(
            f"""// generated by generate_mat_kernels32.py
static_assert(SIMD<double>::Size() == {simd_size}, "inconsistent compile flags for generate_mat_kernels32.py and matkernel32.hpp");

template <size_t H, size_t W, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
INLINE void MatKernelFloat
(size_t nk,
 BareSliceMatrix<float,OA> a,
 BareSliceMatrix<float,OB> b,
 BareSliceMatrix<float,RowMajor> c);

template <size_t H, size_t W, bool ADD, bool POS,
          ORDERING OA, ORDERING OB>
INLINE void MatKernelComplex32
(size_t nk,
 BareSliceMatrix<Complex32,OA> a,
 BareSliceMatrix<Complex32,OB> b,
 BareSliceMatrix<Complex32,RowMajor> c);

"""
        )

        ### Float fixed small-K kernels
        generate_float_small_k(out, simd_size)

        ### Complex32 fixed small-K kernels
        generate_complex32_small_k(out, simd_size)

        ### Float outer-product kernels for all layouts except RowMajor x ColMajor
        for height in range(1, 5):
            width = 1
            while width <= 8 * simd_size:
                for add in (False, True):
                    for pos in (False, True):
                        for ordering_a, ordering_b in (
                            ("ColMajor", "RowMajor"),
                            ("RowMajor", "RowMajor"),
                            ("ColMajor", "ColMajor"),
                        ):
                            generate_float_kernel(
                                out,
                                simd_size,
                                height,
                                width,
                                add,
                                pos,
                                ordering_a,
                                ordering_b,
                            )
                width *= 2

        ### Float inner-product kernels for RowMajor A x ColMajor B
        for height in range(1, 5):
            width = 1
            while width <= 4:
                for add in (False, True):
                    for pos in (False, True):
                        for ordering_a, ordering_b in (
                            ("RowMajor", "ColMajor"),
                        ):
                            generate_float_inner_product_kernel(
                                out,
                                simd_size,
                                height,
                                width,
                                add,
                                pos,
                                ordering_a,
                                ordering_b,
                            )
                width *= 2

        ### General Complex32 kernels for all A and B layouts
        for height in range(1, 5):
            width = 1
            while width <= 4 * simd_size:
                for add in (False, True):
                    for pos in (False, True):
                        for ordering_a in ("ColMajor", "RowMajor"):
                            for ordering_b in ("ColMajor", "RowMajor"):
                                generate_complex32_kernel(
                                    out,
                                    simd_size,
                                    height,
                                    width,
                                    add,
                                    pos,
                                    ordering_a,
                                    ordering_b,
                                )
                width *= 2

        out.write(
            f"""#if defined(__AVX512F__) || defined(__arm64__) || defined(__aarch64__)
constexpr size_t NGBLAS32_FLOAT_MAX_NR = {8 * simd_size};
constexpr size_t NGBLAS32_FLOAT_INNER_MAX_MR = 4;
constexpr size_t NGBLAS32_FLOAT_INNER_MAX_NR = 4;
constexpr size_t NGBLAS32_COMPLEX_MAX_NR = {4 * simd_size};
#else
constexpr size_t NGBLAS32_FLOAT_MAX_NR = {4 * simd_size};
constexpr size_t NGBLAS32_FLOAT_INNER_MAX_MR = 2;
constexpr size_t NGBLAS32_FLOAT_INNER_MAX_NR = 2;
constexpr size_t NGBLAS32_COMPLEX_MAX_NR = {2 * simd_size};
#endif
"""
        )


################################################################################
# Command-line interface
################################################################################


def main():
    """Run the command-line interface for the kernel generator.

    The command expects an output header path and a supported SIMD size.
    """
    print(sys.argv)
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} output.hpp simd-size", file=sys.stderr)
        return 1

    try:
        simd_size = int(sys.argv[2])
    except ValueError:
        print(f"invalid SIMD size {sys.argv[2]}", file=sys.stderr)
        return 1

    if simd_size not in (1, 2, 4, 8):
        print(f"unsupported SIMD size {simd_size}", file=sys.stderr)
        return 1

    try:
        generate(sys.argv[1], simd_size)
    except OSError:
        print(f"cannot open {sys.argv[1]}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
