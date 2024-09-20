#ifndef FILE_TENSORCOEFFICIENT_HPP
#define FILE_TENSORCOEFFICIENT_HPP

/*********************************************************************/
/* File:   tensorcoefficient.hpp                                     */
/* Author: Matthias Rambausek                                        */
/* Date:   01. Apr. 2022                                             */
/*********************************************************************/

// #include "fem.hpp"

#include "coefficient.hpp"
#include "symbolicintegrator.hpp"

namespace ngfem {

    namespace tensor_internal {
        using namespace ngcore;
        using namespace ngstd;
        using namespace ngbla;

        class NGS_DLL_HEADER OutOfIndices : public Exception
        { public: using Exception::Exception; };
        
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

            void DoArchive(Archive &ar) {
                size_t isize;
                string tmp_symbol{};
                if(ar.Output()) {
                    isize = indices.Size();
                    ar & isize;
                    for (auto idx : indices) {
                      tmp_symbol = string{idx.symbol};
                      ar & tmp_symbol & idx.pos & idx.dim;
                    }
                }
                else {
                    ar & isize;
                    indices.SetSize(isize);
                    size_t tmp_pos;
                    size_t tmp_idx;
                    for (auto i : Range(isize)) {
                      ar & tmp_symbol & tmp_pos & tmp_idx;
                      indices[i] = Index{tmp_symbol.at(0), tmp_pos, tmp_idx};
                    }
                }
                ar & strides & total_dim;
            }
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

        Array <MultiIndex> compute_multi_indices(const string &index_signature,
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

        optional<string> substitute_id_index(
            string signature, pair<char, char> id_indices, size_t id_pos,
            const FlatArray<bool> marked_for_removal, bool recurse = false);

        vector<string> split_signature(string input);
        
        string form_index_signature(const vector<string>& parts);

        string expand_ellipses(const string& signature, const Array<shared_ptr<CoefficientFunction>>& cfs);

        Vector<bool> nonzero_pattern(shared_ptr <CoefficientFunction> cf);

        pair<string, Array<shared_ptr<CoefficientFunction>>>
        flatten_einsum(const string& signature,
                       const Array<shared_ptr<CoefficientFunction>>& cfs,
                       const map<string, bool> &options);
        
        shared_ptr<CoefficientFunction>
        optimize_legacy(const string& signature,
                      const Array<shared_ptr<CoefficientFunction>>& cfs,
                      const map<string, bool> &options);

        shared_ptr<CoefficientFunction>
        optimize_transpose(const string& signature,
                           const Array<shared_ptr<CoefficientFunction>>& cfs,
                           const map<string, bool> &aoptions);

        shared_ptr<CoefficientFunction>
        optimize_path(const string &signature,
                      const Array<shared_ptr<CoefficientFunction>>& input_cfs,
                      const map<string, bool> &aoptions);
        
        pair<string, Array<shared_ptr<CoefficientFunction>>>
        expand_higher_order_identities(string signature,
                       const Array<shared_ptr<CoefficientFunction>>& cfs,
                       [[maybe_unused]] const map<string, bool> &options);
        
        pair<string, Array<shared_ptr<CoefficientFunction>>>
        optimize_identities(string, const Array<shared_ptr<CoefficientFunction>>& cfs,
                            const map<string, bool> &options);


        // class LeviCivitaCoefficientFunction
        //     : public T_CoefficientFunction<LeviCivitaCoefficientFunction> {
        //   using BASE = T_CoefficientFunction<LeviCivitaCoefficientFunction>;

        //   int dim = 0;
        //   MultiIndex mi{};

        // public:
        //   LeviCivitaCoefficientFunction() = default;

        //   LeviCivitaCoefficientFunction(int adim);

        //   virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override {
        //     func(*this);
        //   }

        //   virtual string GetDescription() const override { return string("Levi-Civita Symbol"); }

        //   virtual void DoArchive(Archive &ar) override;

        //   virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
        //     GenerateCode(code, inputs, index, false);
        //   }

        //   virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index,
        //                             bool skip_zeroes = true) const;

        //   virtual void NonZeroPattern(const class ProxyUserData &ud,
        //                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override;

        //   virtual void NonZeroPattern(const class ProxyUserData &ud,
        //                               FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
        //                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override;

        //   using BASE::Evaluate;
        //   using typename BASE::T_DJC;

        //   virtual double Evaluate(const BaseMappedIntegrationPoint &ip) const override;

        //   template<typename MIR, typename T, ORDERING ORD>
        //   void T_Evaluate(const MIR &ir, BareSliceMatrix<T, ORD> values) const
        //   {
        //     auto val = T(0.0);
        //     values.AddSize(Dimension(), ir.Size()) = val;
        //     auto ir_size = ir.Size();
        //     for (size_t I: Range(Dimension())) {
        //       const auto I_array = split(I, mi);
        //       cout << "check I: " << I_array << endl;
        //       if (is_even_iota_permutation(I_array.begin(), I_array.end()))
        //         val = 1.0;
        //       else if (is_odd_iota_permutation(I_array.begin(), I_array.end()))
        //         val = -1.0;
        //       else
        //         continue;
        //       for (auto q: Range(ir_size))
        //         values(I, q) = val;
        //     }
        //   }

        //   template<typename MIR, typename T, ORDERING ORD>
        //   void T_Evaluate(const MIR &ir, FlatArray<BareSliceMatrix<T, ORD>> input,
        //                   BareSliceMatrix<T, ORD> values) const
        //   {
        //     T_Evaluate(ir, values);
        //   }


        //   shared_ptr<CoefficientFunction> Diff(const CoefficientFunction *var,
        //                                        shared_ptr<CoefficientFunction> dir) const override;

        //   shared_ptr<CoefficientFunction> DiffJacobi(const CoefficientFunction *var, T_DJC & cache) const override;
        // };

        class EinsumCoefficientFunction
            : public T_CoefficientFunction<EinsumCoefficientFunction>
        {
          using BASE = T_CoefficientFunction<EinsumCoefficientFunction>;
          using typename BASE::T_DJC;

          static constexpr bool sparse_evaluation_default = true;

          bool is_zero{false};
          size_t max_mem{0};
          map<string, bool> options{};

          Array<Vector<bool>> nz_inputs{};
          Vector<bool> nz_result{};
          Vector<bool> nz_all{};
          Matrix<int> index_maps{};
          Matrix<int> sparse_index_maps{};

          shared_ptr<CoefficientFunction> node{};

          string index_signature{};
          Array<shared_ptr<CoefficientFunction>> cfs{};

          string original_index_signature{};
          Array<shared_ptr<CoefficientFunction>> original_inputs{};

          string expanded_index_signature{};
          Array<shared_ptr<CoefficientFunction>> expanded_inputs{};

        public:
          EinsumCoefficientFunction() = default;

          EinsumCoefficientFunction(const string &aindex_signature,
                                    const Array<shared_ptr<CoefficientFunction>> &acfs,
                                    const map<string, bool> &aoptions);

        private:
          Matrix<int> build_index_maps(const Array<MultiIndex>& index_sets, const optional<Vector<bool>>& nz_pattern);

        public:

          virtual void DoArchive(Archive &ar) override;

          virtual shared_ptr<EinsumCoefficientFunction> Optimize(const map<string, bool> &aoptions) const;

          const string &IndexSignature() const { return index_signature; }

          const string &ExpandedIndexSignature() const {return expanded_index_signature;}

          const string &OriginalIndexSignature() const {return original_index_signature;}

          virtual Array<shared_ptr<CoefficientFunction>>
          InputCoefficientFunctions() const override
          {
            if (node)
              return node->InputCoefficientFunctions();
            return Array<shared_ptr<CoefficientFunction>>{cfs};
          }

          Array<shared_ptr<CoefficientFunction>>
          ExpandedInputCoefficientFunctions() const
          {
            return Array<shared_ptr<CoefficientFunction>>{expanded_inputs};
          }

          Array<shared_ptr<CoefficientFunction>>
          OriginalInputCoefficientFunctions() const
          {
            return Array<shared_ptr<CoefficientFunction>>{original_inputs};
          }

          shared_ptr<CoefficientFunction> OptimizedNode() const
          {
            return node;
          }

          virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override
          {
            const auto incfs = InputCoefficientFunctions();
            for_each(incfs.begin(), incfs.end(), [&](const auto &cf) { cf->TraverseTree(func); });
            func(*this);
          }

          virtual string GetDescription() const override;

          virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
            GenerateCode(code, inputs, index, false);
          }

          virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index,
                                    bool skip_zeroes = true) const;

          virtual void NonZeroPattern(const class ProxyUserData &ud,
                                      FlatVector<AutoDiffDiff<1,NonZero>> values) const override;

          virtual void NonZeroPattern(const class ProxyUserData &ud,
                                      FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                      FlatVector<AutoDiffDiff<1,NonZero>> values) const override;

          using BASE::Evaluate;

          virtual double Evaluate(const BaseMappedIntegrationPoint &ip) const override {
            if (Dimension() == 1)
              return BASE::Evaluate(ip);
            throw Exception("TensorProductCF scalar evaluate called for non-scalar result");
          }

          template<typename MIR, typename T, ORDERING ORD>
          void T_Evaluate(const MIR &ir, BareSliceMatrix<T, ORD> values) const
          {
            if (node) {
              node->Evaluate(ir, values);
              return;
            }

            ArrayMem<T, 1000> mem(max_mem * ir.Size());
            T *mem_pos = mem.Data();
            Array<FlatMatrix<T, ORD>> tmp_arrays(cfs.Size());
            for (size_t i: Range(cfs)) {
              tmp_arrays[i].AssignMemory(cfs[i]->Dimension(), ir.Size(), mem_pos);
              mem_pos += tmp_arrays[i].Height() * tmp_arrays[i].Width();
              cfs[i]->Evaluate(ir, tmp_arrays[i]);
            }

            values.AddSize(Dimension(), ir.Size()) = T(0.0);
            const auto cres = cfs.Size();

            const auto& I_maps = sparse_index_maps.Height() > 0 ? sparse_index_maps : index_maps;
            for (size_t I: Range(I_maps.Height())) {
              const auto& I_map = I_maps.Row(I);
              for (int q: Range(ir)) {
                T tmp(1.0);
                for (size_t i: Range(tmp_arrays))
                  tmp *= tmp_arrays[i](I_map(i), q);
                values(I_map(cres), q) += tmp;
              }
            }
          }

          template<typename MIR, typename T, ORDERING ORD>
          void T_Evaluate(const MIR &ir, FlatArray<BareSliceMatrix<T, ORD>> input,
                          BareSliceMatrix<T, ORD> values) const
          {
            if (node) {
              node->Evaluate(ir, input, values);
              return;
            }

            values.AddSize(Dimension(), ir.Size()) = T(0.0);
            const auto cres = cfs.Size();

            const auto& I_maps = sparse_index_maps.Height() > 0 ? sparse_index_maps : index_maps;
            for (size_t I: Range(I_maps.Height()))
            {
              const auto& I_map = I_maps.Row(I);
              for (int q: Range(ir)) {
                T tmp(1.0);
                for (size_t i: Range(input))
                  tmp *= input[i](I_map(i), q);
                values(I_map(cres), q) += tmp;
              }
            }
          }

          virtual shared_ptr<CoefficientFunction> Diff(
              const CoefficientFunction *var,
              shared_ptr<CoefficientFunction> dir) const override;

          virtual shared_ptr<CoefficientFunction> DiffJacobi(
              const CoefficientFunction *var, T_DJC & cache) const override;

          virtual bool IsZeroCF() const override { return is_zero; }
        };
    }

} // namespace ngfem


#endif //FILE_TENSORCOEFFICIENT_HPP
