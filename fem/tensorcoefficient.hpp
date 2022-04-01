#ifndef FILE_TENSORCOEFFICIENT_HPP
#define FILE_TENSORCOEFFICIENT_HPP

/*********************************************************************/
/* File:   tensorcoefficient.hpp                                     */
/* Author: Matthias Rambausek                                        */
/* Date:   01. Apr. 2022                                             */
/*********************************************************************/

#include "coefficient.hpp"


namespace ngfem {

    namespace tensor_internal {

        struct Index {
            char symbol;
            size_t pos;
            size_t dim;
        };

        bool operator==(const Index &a, const Index &b) {
            return a.symbol == b.symbol && a.pos == b.pos && a.dim == b.dim;
        }

        bool operator!=(const Index &a, const Index &b) { return !(a == b); }

        class MultiIndex {
            Array <Index> indices{};
            Array <size_t> strides{};
            size_t total_dim = 1;

        public:
            size_t Append(Index idx) {
                for (int i = strides.Size() - 1; i >= 0; --i)
                    strides[i] *= idx.dim;
                strides.Append(1);
                total_dim *= idx.dim;
                indices.Append(idx);
                return total_dim;
            }

            const Index &operator[](size_t i) const { return indices[i]; }

            void RelabelIndex(size_t i, char symbol) { indices[i].symbol = symbol; }

            void ReplaceIndex(size_t i, Index idx) { indices[i] = idx; }

            void RepositionIndex(size_t i, size_t pos) { indices[i].pos = pos; }

            void MoveIndex(size_t i, int diff) {
                if (indices[i].pos + diff < 0)
                    throw NG_EXCEPTION("index position must not be smaller than 0");
                indices[i].pos += diff;
            }

            const Array <size_t> &Strides() const { return strides; }

            size_t Size() const { return indices.Size(); };

            size_t TotalDim() const { return total_dim; }
        };

        bool operator==(const MultiIndex &a, const MultiIndex &b) {
            if (a.Size() != b.Size())
                return false;

            for (size_t i: Range(a.Size()))
                if (a[i] != b[i])
                    return false;
            return true;
        }

        bool operator!=(const MultiIndex &a, const MultiIndex &b) { return !(a == b); }

        auto begin(const MultiIndex &mi) { return &(mi[0]); }

        auto end(const MultiIndex &mi) { return &(mi[0]) + mi.Size(); }

        size_t join(const FlatArray <size_t> full_indices, const MultiIndex &multi_idx);

        Array <size_t> split(size_t full_index, const MultiIndex &multi_idx);

        Array <size_t> split2(size_t full_index, const MultiIndex &multi_idx);

        Array <MultiIndex> compute_index_sets(const string &index_signature,
                                              const Array <shared_ptr<CoefficientFunction>> &cfs);

        template<class ForwardIt>
        bool is_even_iota_permutation(ForwardIt first, ForwardIt last, size_t base = 0) {
            auto base_it = find(first, last, base);
            if (base_it == last)
                return false;
            base = base_it - first;
            for (; base_it != last; ++base_it)
                if ((base_it - first) - base != *base_it)
                    return false;
            auto offset = (last - first) - base;
            for (base_it = first; base_it != first + base; ++base_it)
                if ((base_it - first) + offset != *base_it)
                    return false;
            return true;
        }

        template<class ForwardIt>
        bool is_odd_iota_permutation(ForwardIt first, ForwardIt last, size_t base = 0) {
            auto base_it = find(first, last, base);
            if (base_it == last)
                return false;
            base = base_it - first;
            for (base_it = first; base_it != first + base + 1; ++base_it)
                if (base + (first - base_it) != *base_it)
                    return false;
            for (; base_it != last; ++base_it)
                if ((last - base_it) + base != *base_it)
                    return false;
            return true;
        }

        string new_index_symbols(string existing, size_t nnew = 1);

        string tensor_product_signature(shared_ptr <CoefficientFunction> cf);

        bool is_tensor_product(shared_ptr <CoefficientFunction> cf);

        Vector<bool> nonzero_pattern(shared_ptr <CoefficientFunction> cf);
    }

} // namespace ngfem


#endif //FILE_TENSORCOEFFICIENT_HPP
