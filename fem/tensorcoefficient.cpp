/*********************************************************************/
/* File:   tensorcoefficient.cpp                                     */
/* Author: Matthias Rambausek                                        */
/* Date:   01. Apr. 2022                                             */
/*********************************************************************/

#include <cmath>
#include <fem.hpp>
#include <ngstd.hpp>
#include <python_ngstd.hpp>

#include "tensorcoefficient.hpp"


namespace ngfem {

    namespace tensor_internal {

        bool get_option(const map<string, bool> &options, string key, bool or_val) {
            auto found = options.find(key);
            if (found != options.end())
                return found->second;
            else
                return or_val;
        };

        Vector<bool> nonzero_pattern(CoefficientFunction *cf) {
            Vector<AutoDiffDiff<1, bool>> nzvec_ad(cf->Dimension());
            Vector<bool> nzvec(cf->Dimension());
            nzvec = false;

            DummyFE<ET_TRIG> dummyfe;
            ProxyUserData ud;
            ud.fel = &dummyfe;

            Array<ProxyFunction *> trial_proxies, test_proxies;

            cf->TraverseTree([&](CoefficientFunction &nodecf) {
                auto proxy = dynamic_cast<ProxyFunction *>(&nodecf);
                if (proxy) {
                    if (proxy->IsTestFunction()) {
                        if (!test_proxies.Contains(proxy)) {
                            test_proxies.Append(proxy);
                        }
                    } else if (!trial_proxies.Contains(proxy)) {
                        trial_proxies.Append(proxy);
                    }
                }
            });

            cf->NonZeroPattern(ud, nzvec_ad);
            for (size_t i: Range(nzvec_ad))
                nzvec[i] = nzvec[i] || nzvec_ad[i].Value();

            for (int l1: trial_proxies.Range())
                for (int l2: Range(0, trial_proxies[l1]->Dimension())) {
                    ud.trialfunction = trial_proxies[l1];
                    ud.trial_comp = l2;
                    cf->NonZeroPattern(ud, nzvec_ad);
                    for (size_t i: Range(nzvec_ad))
                        nzvec[i] = nzvec[i] || nzvec_ad[i].Value();

                    for (int k1: test_proxies.Range())
                        for (int k2: Range(0, test_proxies[k1]->Dimension())) {
                            ud.testfunction = test_proxies[k1];
                            ud.test_comp = k2;
                            cf->NonZeroPattern(ud, nzvec_ad);
                            for (size_t i: Range(nzvec_ad))
                                nzvec[i] = nzvec[i] || nzvec_ad[i].Value();
                        }
                }
            return nzvec;
        }

        Vector<bool> nonzero_pattern(shared_ptr<CoefficientFunction> cf) {
            return nonzero_pattern(cf.get());
        }

        string new_index_symbols(string existing, size_t nnew) {
            char new_char = 'A';
            auto target_size = existing.size() + nnew;

            while (existing.size() < target_size) {
                if (new_char > 'Z')
                    break;
                if (find(begin(existing), end(existing), new_char) != end(existing))
                    ++new_char;
                else
                    existing += new_char;
            }

            new_char = 'a';
            while (existing.size() < target_size) {
                if (new_char > 'z')
                    throw NG_EXCEPTION("did not find any unused symbol in A-Z and a-z");
                if (find(begin(existing), end(existing), new_char) != end(existing))
                    ++new_char;
                else
                    existing += new_char;
            }

            return string{existing.end() - nnew, existing.end()};
        }

        string form_index_signature(FlatArray<MultiIndex> index_sets) {
            if (index_sets.Size() < 2)
                throw NG_EXCEPTION("there must at least be two index sets -- one input and one output");

            size_t n_inputs = index_sets.Size() - 1;

            stringstream str;
            for (size_t j: Range(index_sets[0].Size()))
                str << index_sets[0][j].symbol;

            for (size_t i: Range(size_t{1}, n_inputs)) {
                str << ',';
                for (size_t j: Range(index_sets[i].Size()))
                    str << index_sets[i][j].symbol;
            }

            str << "->";
            const auto result = index_sets[n_inputs];
            for (size_t j: Range(result.Size()))
                str << result[j].symbol;

            return str.str();
        }

        size_t join(const FlatArray<size_t> full_indices, const MultiIndex &multi_idx) {
            const auto &strides = multi_idx.Strides();
            size_t res = 0;
            for (size_t i: Range(multi_idx.Size()))
                res += full_indices[multi_idx[i].pos] * strides[i];
            return res;
        }

        Array<size_t> split(size_t full_index, const MultiIndex &multi_idx) {
            const auto &strides = multi_idx.Strides();
            Array<size_t> indices(multi_idx.Size());
            for (int i: Range(multi_idx.Size())) {
                indices[i] = full_index / strides[i];
                full_index -= indices[i] * strides[i];
            }
            return move(indices);
        }

        Array<size_t> split2(size_t full_index, const MultiIndex &multi_idx) {
            Array<size_t> indices(multi_idx.Size());
            for (int i = multi_idx.Size() - 1; i >= 0; i--) {
                indices[i] = full_index % multi_idx[i].dim;
                full_index /= multi_idx[i].dim;
            }
            return move(indices);
        }

        // index format: 'ij,jk->ik' (result after '->' is mandatory!); return: input data, result data,
        // full set
        Array<MultiIndex> compute_multi_indices(const string &index_signature,
                                                const Array<shared_ptr<CoefficientFunction>> &cfs) {
            Array<MultiIndex> idx_sets(cfs.Size() + 2);
            MultiIndex *current_idx_set = idx_sets.begin();
            MultiIndex &full_idx_set = idx_sets[cfs.Size() + 1];

            Array<char> chars;
            Array<int> char_counts;
            Array<char> free_symbols;

            int subset_cnt = 0;
            int subset_idx = 0;
            bool result_section = false;

            for_each(begin(index_signature), end(index_signature), [&](const char idx) {
                if (idx == '>' || idx == '\0')
                    return;

                if (idx == '-')
                    result_section = true;

                if (idx == ',' || idx == '-') {
                    ++current_idx_set;
                    ++subset_cnt;
                    subset_idx = 0;
                    return;
                }

                if ((idx < 'a' || idx > 'z') && (idx < 'A' || idx > 'Z'))
                    throw NG_EXCEPTION("invalid index character detected -- use only a-z, A-Z");

                if (!chars.Contains(idx)) {
                    if (result_section)
                        throw NG_EXCEPTION("detected 'result indices' not present in the inputs");

                    size_t subset_dim_i = 1;
                    if (cfs[subset_cnt]->Dimensions().Size() <= 1)
                        if (subset_idx > 1)
                            throw NG_EXCEPTION("tensor order mismatch");
                        else
                            subset_dim_i = cfs[subset_cnt]->Dimension();
                    else if (subset_idx + 1 > cfs[subset_cnt]->Dimensions().Size())
                        throw NG_EXCEPTION("tensor order mismatch");
                    else
                        subset_dim_i = cfs[subset_cnt]->Dimensions()[subset_idx];

                    full_idx_set.Append(Index{idx, chars.Size(), subset_dim_i});
                    current_idx_set->Append(full_idx_set[chars.Size()]);
                    chars.Append(idx);
                    char_counts.Append(1);
                } else if (!result_section && cfs[subset_cnt]->Dimensions()[subset_idx] !=
                                              full_idx_set[chars.Pos(idx)].dim) {
                    throw NG_EXCEPTION("dimensions of contracted indices do not match");
                } else {
                    current_idx_set->Append(full_idx_set[chars.Pos(idx)]);
                    if (!result_section)
                        char_counts[chars.Pos(idx)] += 1;
                }

                ++subset_idx;
            });

            return move(idx_sets);
        }

        Array<int> index_dimensions(const MultiIndex &mi) {
            Array<int> res_dims(mi.Size());
            for (size_t i: Range(res_dims.Size()))
                res_dims[i] = mi[i].dim;
            return res_dims;
        }

        bool substitute_id_index(FlatArray<MultiIndex> index_sets, Index id_idx_0, Index id_idx_1,
                                 size_t id_pos, const FlatArray<bool> marked_for_removal,
                                 bool recurse = false) {
            bool subs{false};
            MultiIndex &result_idx = index_sets[index_sets.Size() - 2];
            if (find_if(begin(result_idx), end(result_idx), [&](const auto &item) {
                return item.symbol == id_idx_0.symbol;
            }) != end(result_idx)) {
                if (recurse)
                    return substitute_id_index(index_sets, id_idx_1, id_idx_0, id_pos,
                                               marked_for_removal);
                else
                    return false;
            }

            for (size_t i2: Range(index_sets.Size() - 2)) {
                // do not touch current id tensor pos
                if (i2 == id_pos || marked_for_removal[i2])
                    continue;

                auto &idx2 = index_sets[i2];
                for (size_t k: Range(idx2.Size())) {
                    // substitute based on symbol
                    if (idx2[k].symbol == id_idx_0.symbol) {
                        idx2.ReplaceIndex(k, id_idx_1);
                        subs = true;
                    }
                }
            }

            if (subs) {
                // remove from full index set
                auto full_idx = index_sets[index_sets.Size() - 1];
                MultiIndex new_full_idx;
                for (size_t k: Range(new_full_idx.Size()))
                    if (full_idx[k].symbol != id_idx_0.symbol)
                        new_full_idx.Append(full_idx[k]);
                index_sets[index_sets.Size() - 1] = new_full_idx;

                // update index positions
                for (auto &idx2: index_sets)
                    for (size_t k: Range(idx2.Size()))
                        if (idx2[k].pos > id_idx_0.pos)
                            idx2.RepositionIndex(k, idx2[k].pos - 1);
            } else if (recurse)
                return substitute_id_index(index_sets, id_idx_1, id_idx_0, id_pos, marked_for_removal);

            return subs;
        }
    }

    using namespace tensor_internal;


    class LeviCivitaCoefficientFunction
            : public T_CoefficientFunction<LeviCivitaCoefficientFunction> {
        using BASE = T_CoefficientFunction<LeviCivitaCoefficientFunction>;

        int dim = 0;
        MultiIndex mi{};

    public:
        LeviCivitaCoefficientFunction() = default;

        LeviCivitaCoefficientFunction(int adim) : BASE(1), dim{adim} {

            Array<int> dims(dim);
            char symbol = 'a';
            for (size_t i: Range(dim)) {
                dims[i] = dim;
                mi.Append(Index{symbol++, i, static_cast<size_t>(dim)});
            }

            SetDimensions(dims.Part(0));
        }

        virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override {
            func(*this);
        }

        virtual string GetDescription() const override { return string("Levi-Civita Symbol"); }

        void DoArchive(Archive &ar) override { BASE::DoArchive(ar); }

        virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
            GenerateCode(code, inputs, index, false);
        }

        virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index,
                                  bool skip_zeroes = true) const {
            LocalHeap lh(100000, "Levi-Cevita code gen");

            auto dims = Dimensions();

            auto nzvec = nonzero_pattern(
                    const_cast<LeviCivitaCoefficientFunction *>(this)->shared_from_this());

            for (size_t I: Range(Dimension())) {

                if (skip_zeroes && !nzvec[I])
                    continue;

                const auto I_array = split(I, mi);

                if (is_even_iota_permutation(I_array.begin(), I_array.end()))
                    code.body += Var(index, I, dims).Assign(Var(1.0));
                else if (is_odd_iota_permutation(I_array.begin(), I_array.end()))
                    code.body += Var(index, I, dims).Assign(Var(-1.0));
                else
                    code.body += Var(index, I, dims).Assign(Var(0.0));
            }
        }

        virtual void NonZeroPattern(const class ProxyUserData &ud,
                                    FlatVector<AutoDiffDiff<1, bool>> values) const override {
            for (size_t I: Range(Dimension())) {
                const auto I_array = split(I, mi);

                if (is_even_iota_permutation(I_array.begin(), I_array.end()))
                    values(I) = true;
                else if (is_odd_iota_permutation(I_array.begin(), I_array.end()))
                    values(I) = true;
                else
                    values(I) = false;
            }
        }

        virtual void NonZeroPattern(const class ProxyUserData &ud,
                                    FlatArray<FlatVector<AutoDiffDiff<1, bool>>> input,
                                    FlatVector<AutoDiffDiff<1, bool>> values) const override {
            NonZeroPattern(ud, values);
        }

        using BASE::Evaluate;

        virtual double Evaluate(const BaseMappedIntegrationPoint &ip) const override {
            if (Dimension() == 1)
                return BASE::Evaluate(ip);
            throw Exception("LeviCivitaCF scalar evaluate called for non-scalar result");
        }

        template<typename MIR, typename T, ORDERING ORD>
        void T_Evaluate(const MIR &ir, BareSliceMatrix<T, ORD> values) const {
            auto val = T(0.0);
            values.AddSize(Dimension(), ir.Size()) = val;
            auto ir_size = ir.Size();
            for (size_t I: Range(Dimension())) {
                const auto I_array = split(I, mi);
                if (is_even_iota_permutation(I_array.begin(), I_array.end()))
                    val = 1.0;
                else if (is_odd_iota_permutation(I_array.begin(), I_array.end()))
                    val = -1.0;
                else
                    continue;
                for (auto q: Range(ir_size))
                    values(I, q) = val;
            }
        }

        template<typename MIR, typename T, ORDERING ORD>
        void T_Evaluate(const MIR &ir, FlatArray<BareSliceMatrix<T, ORD>> input,
                        BareSliceMatrix<T, ORD> values) const {
            T_Evaluate(ir, values);
        }

        shared_ptr<CoefficientFunction> Diff(const CoefficientFunction *var,
                                             shared_ptr<CoefficientFunction> dir) const override {
            // NOTE (MR): How to handle symmetry here?
            if (this == var)
                return dir;

            return ZeroCF(Array<int>{Dimensions()});
        }

        shared_ptr<CoefficientFunction> DiffJacobi(const CoefficientFunction *var) const override {
            // NOTE (MR): How to handle symmetry here?
            if (this == var)
                return IdentityCF(Array<int>{Dimensions()});

            Array<int> dims{Dimensions()};
            dims.Append(Dimensions());
            return ZeroCF(dims);
        }
    };

    shared_ptr<CoefficientFunction> LeviCivitaCF(int dimension) {
        return make_shared<LeviCivitaCoefficientFunction>(dimension);
    }

    pair<Array<MultiIndex>, Array<shared_ptr<CoefficientFunction>>>
    optimize_special_inputs(Array<MultiIndex> index_sets, Array<shared_ptr<CoefficientFunction>> cfs,
                            const map<string, bool> &options);

    shared_ptr<CoefficientFunction>
    optimize_blas(FlatArray<MultiIndex> index_sets, FlatArray<shared_ptr<CoefficientFunction>> cfs,
                  const map<string, bool> &options) {
        for (size_t i: Range(cfs))
            if (cfs[i]->IsZeroCF())
                return ZeroCF(index_dimensions(index_sets[cfs.Size()]));

        const auto identity_descr = IdentityCF(1)->GetDescription();
        cout << "TP: trying to detect some BLAS operations" << endl;

        // Trace
        if (cfs.Size() == 1 && index_sets[0].Size() == 2) {
            if (index_sets[0][0].symbol == index_sets[0][1].symbol && index_sets[1].Size() == 0)
                return TraceCF(cfs[0]);
        }

        // Transpose
        if (cfs.Size() == 1 && index_sets[0].Size() == 2 && index_sets[1].Size() == 2) {
            if (index_sets[0][0].symbol != index_sets[0][1].symbol &&
                index_sets[1][0].symbol == index_sets[0][1].symbol)
                return TransposeCF(cfs[0]);
        }

        // Mat * Vec
        if (cfs.Size() == 2 && index_sets[0].Size() == 2 && index_sets[1].Size() == 1 &&
            index_sets[2].Size() == 1 && index_sets[0][0].symbol != index_sets[0][1].symbol &&
            index_sets[1][0].symbol == index_sets[0][1].symbol) {
            // NOTE: no other way to detect identity!
            if (cfs[0]->GetDescription() == identity_descr) {
                cout << "EinsumCF: detected I * vec" << endl;
                return cfs[1];
            }
            cout << "EinsumCF: detected Mat * Vec" << endl;
            return cfs[0] * cfs[1];
        }

        // Mat * Mat
        if (cfs.Size() == 2 && index_sets[0].Size() == 2 && index_sets[1].Size() == 2 &&
            index_sets[2].Size() == 2 && index_sets[0][0].symbol != index_sets[0][1].symbol &&
            index_sets[1][0].symbol != index_sets[1][1].symbol) {

            bool ij_jk = index_sets[1][0].symbol == index_sets[0][1].symbol &&
                         index_sets[2][0].symbol == index_sets[0][0].symbol;

            bool jk_ij = index_sets[0][0].symbol == index_sets[1][1].symbol &&
                         index_sets[2][0].symbol == index_sets[1][0].symbol;

            if (ij_jk || jk_ij) {
                int i0 = ij_jk ? 0 : 1;
                int i1 = ij_jk ? 1 : 0;

                if (cfs[i0]->GetDescription() == identity_descr) {
                    cout << "EinsumCF: detected I * Mat" << endl;
                    return cfs[i1];
                }

                if (cfs[i1]->GetDescription() == identity_descr) {
                    cout << "EinsumCF: detected Mat * I" << endl;
                    return cfs[i0];
                }

                cout << "EinsumCF: detected Mat * Mat" << endl;
                return cfs[i0] * cfs[i1];
            }
        }

        // Mat.trans * Mat
        if (cfs.Size() == 2 && index_sets[0].Size() == 2 && index_sets[1].Size() == 2 &&
            index_sets[2].Size() == 2 && index_sets[0][0].symbol != index_sets[0][1].symbol &&
            index_sets[1][0].symbol != index_sets[1][1].symbol &&
            index_sets[1][0].symbol == index_sets[0][0].symbol &&
            index_sets[2][0].symbol == index_sets[0][1].symbol) {
            if (cfs[0]->GetDescription() == identity_descr) {
                cout << "EinsumCF: detected I * Mat" << endl;
                return cfs[1];
            }

            if (cfs[1]->GetDescription() == identity_descr) {
                cout << "EinsumCF: detected Mat.trans * I" << endl;
                return TransposeCF(cfs[0]);
            }

            cout << "EinsumCF: detected Mat.trans * Mat" << endl;
            return TransposeCF(cfs[0]) * cfs[1];
        }

        // Mat * Mat.trans
        if (cfs.Size() == 2 && index_sets[0].Size() == 2 && index_sets[1].Size() == 2 &&
            index_sets[2].Size() == 2 && index_sets[0][0].symbol != index_sets[0][1].symbol &&
            index_sets[1][0].symbol != index_sets[1][1].symbol &&
            index_sets[1][1].symbol == index_sets[0][1].symbol &&
            index_sets[2][0].symbol == index_sets[0][0].symbol) {
            if (cfs[0]->GetDescription() == identity_descr) {
                cout << "EinsumCF: detected I * Mat.trans" << endl;
                return TransposeCF(cfs[1]);
            }

            if (cfs[1]->GetDescription() == identity_descr) {
                cout << "EinsumCF: detected Mat * I" << endl;
                return cfs[0];
            }

            cout << "EinsumCF: detected Mat * Mat.trans" << endl;
            return cfs[0] * TransposeCF(cfs[1]);
        }
        return {};
    }

    shared_ptr<CoefficientFunction>
    optimize_path(const string &signature, FlatArray<shared_ptr<CoefficientFunction>> input_cfs,
                  const map<string, bool> &aoptions) {

        map<string, bool> options = aoptions;
        options["optimize_path"] = false;
        options["expand_einsum"] = false;

        namespace py = pybind11;

        auto np = py::module::import("numpy");

        using namespace pybind11::literals;
        py::object einsum_path = np.attr("einsum_path");

        py::list inputs{};
        for (auto icf: input_cfs) {
            py::array::ShapeContainer shape{};
            for (int dim: icf->Dimensions())
                shape->push_back(dim);
            inputs.append(py::array_t<double>(shape));
        }
        auto res = einsum_path(signature, *inputs, "einsum_call"_a = true);
        auto res_tuple = py::extract<py::tuple>(res)();
        auto path = py::extract<py::list>(res_tuple[1])();
        Array<shared_ptr<CoefficientFunction>> tp_inputs{input_cfs};
        for (size_t j: Range(path.size())) {
            auto op = py::extract<py::tuple>(path[j])();
            Array<shared_ptr<CoefficientFunction>> new_inputs;
            Array<shared_ptr<CoefficientFunction>> old_inputs;
            auto in_indices = py::extract<py::tuple>(op[0])();
            Array<bool> drop_from_old(tp_inputs.Size());
            drop_from_old = false;
            for (const auto &item: in_indices) {
                auto in_idx = py::extract<size_t>(item)();
                new_inputs.Append(tp_inputs[in_idx]);
                drop_from_old[in_idx] = true;
            }
            for (size_t i: Range(tp_inputs.Size()))
                if (!drop_from_old[i])
                    old_inputs.Append(tp_inputs[i]);

            string new_signature = py::extract<py::str>(op[2])();
            tp_inputs = old_inputs;
            tp_inputs.Append(EinsumCF(new_signature, move(new_inputs), options));
        }
        return tp_inputs.Last();
    }

    class EinsumCoefficientFunction
            : public T_CoefficientFunction<EinsumCoefficientFunction> {
        using BASE = T_CoefficientFunction<EinsumCoefficientFunction>;

        Array<MultiIndex> index_sets{};
        Array<shared_ptr<CoefficientFunction>> cfs{};
        shared_ptr<CoefficientFunction> node;
        string index_signature{};
        size_t max_mem{0};
        map<string, bool> options;
        Array<Vector<bool>> nz_inputs{};
        Vector<bool> nz_result{};
        Vector<bool> nz_all{};
        Array<Array<int>> input_indices{};
        Array<int> result_indices{};

        const MultiIndex &result_idx() const { return index_sets[cfs.Size()]; }

        const MultiIndex &full_idx() const { return index_sets[cfs.Size() + 1]; }

    public:
        EinsumCoefficientFunction() = default;

        EinsumCoefficientFunction(const string &aindex_signature,
                                  const Array<shared_ptr<CoefficientFunction>> &acfs,
                                  const map<string, bool> &aoptions)
                : BASE(1,
                       std::find_if(acfs.begin(), acfs.end(),
                                    [](const auto &item) { return item->IsComplex(); }) != acfs.end()),
                  index_sets{compute_multi_indices(index_signature, cfs)}, cfs{acfs},
                  index_signature{aindex_signature}, max_mem{0}, node{}, options{aoptions} {

            if (get_option(options, "use_blas_ops", false))
                node = optimize_blas(index_sets, cfs, options);

            if (!node) {
                auto optimized = optimize_special_inputs(move(index_sets), move(cfs), options);
                index_sets = optimized.first;
                cfs = optimized.second;

                if (cfs.Size() == 1 && index_sets[0] == index_sets[1])
                    node = cfs[0];
                else if (get_option(options, "optimize_path", true))
                    node = optimize_path(index_signature, cfs, options);
            }

            if (result_idx().Size() > 0)
                SetDimensions(index_dimensions(result_idx()).Part(0));

            for (size_t i: Range(cfs))
                max_mem += index_sets[i].TotalDim();

            if (!node && get_option(options, "precompute_index_operations", true)) {
                nz_inputs.SetSize(cfs.Size());
                for (size_t i: Range(cfs)) {
                    nz_inputs[i] = nonzero_pattern(cfs[i]);
                    input_indices.SetSize(full_idx().TotalDim());
                }
                nz_result = nonzero_pattern(this);
                nz_all.SetSize(full_idx().TotalDim());
                result_indices.SetSize(full_idx().TotalDim());

                const auto &RI = result_idx();
                for (size_t I: Range(full_idx().TotalDim())) {
                    const auto I_array = split(I, full_idx());
                    if (!nz_result(join(I_array, RI))) {
                        nz_all(I) = false;
                        continue;
                    }
                    result_indices[I] = join(I_array, RI);
                    for (size_t i: Range(cfs)) {
                        if (!nz_inputs[i](join(I_array, index_sets[i]))) {
                            nz_all(I) = false;
                            continue;
                        }
                        input_indices[i][I] = join(I_array, index_sets[i]);
                    }
                }
            }
        }

        const string &IndexSignature() const { return index_signature; }

        virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override {
            for_each(cfs.begin(), cfs.end(), [&](const auto &cf) { cf->TraverseTree(func); });
            func(*this);
        }

        virtual string GetDescription() const override {
            string descr = "EinsumCF ";
            descr += index_signature;
            if (node)
                descr += " [ optimized node (top): " + node->GetDescription() + " ]";
            return descr;
        }

        void DoArchive(Archive &ar) override {
            BASE::DoArchive(ar);
            for_each(cfs.begin(), cfs.end(), [&](auto cf) { ar.Shallow(cf); });
        }

        virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
            GenerateCode(code, inputs, index, false);
        }

        virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index,
                                  bool skip_zeroes = true) const {
            

            LocalHeap lh(100000, "TP code gen");
            const auto &RI = result_idx();

            Array<FlatArray<int>> input_dims(cfs.Size());
            Array<Vector<bool>> nz_vecs_mem;
            FlatArray<Vector<bool>> nz_vecs;
            if (nz_inputs.Size())
                nz_vecs.Assign(nz_inputs);
            else
                nz_vecs_mem.SetSize(cfs.Size());

            for (size_t i: Range(cfs)) {
                if (nz_vecs_mem.Size())
                    nz_vecs[i] = nonzero_pattern(cfs[i]);
                input_dims[i].Assign(index_sets[i].Size(), lh);
                for (size_t j: Range(input_dims[i]))
                    input_dims[i][j] = index_sets[i][j].dim;
            }
            if (nz_vecs_mem.Size())
                nz_vecs.Assign(nz_vecs);

            Vector<bool> result_nz_mem{};
            FlatVector<bool> result_nz_vec;

            if (nz_result.Size())
                result_nz_vec.AssignMemory(nz_result.Size(), nz_result.Data());
            else {
                auto this_cv = const_cast<EinsumCoefficientFunction *>(this);
                result_nz_mem = nonzero_pattern(this_cv->shared_from_this());
                result_nz_vec.AssignMemory(result_nz_mem.Size(), result_nz_mem.Data());
            }

            FlatArray<int> result_dims(RI.Size(), lh);
            for (size_t j: Range(RI.Size()))
                result_dims[j] = RI[j].dim;

            FlatArray<bool> declared(RI.TotalDim(), lh);
            declared = false;

            for (size_t I: Range(full_idx().TotalDim())) {
                const auto I_array = split(I, full_idx());
                const auto res_idx = join(I_array, RI);

                if (skip_zeroes && !result_nz_vec[res_idx])
                    continue;

                CodeExpr s;
                if (!result_nz_vec[res_idx] && !declared[res_idx])
                    s = Var(0.0);
                else
                    for (size_t i: Range(inputs)) {
                        const auto I_local = join(I_array, index_sets[i]);
                        if (!nz_vecs[i][I_local]) {
                            s.code.clear();
                            break;
                        }
                        if (s.code.empty())
                            s = Var(inputs[i], I_local, input_dims[i]);
                        else
                            s *= Var(inputs[i], I_local, input_dims[i]);
                    }

                if (s.code.empty())
                    continue;

                if (declared[res_idx])
                    code.body += Var(index, res_idx, result_dims)
                            .Assign(Var(index, res_idx, result_dims) + s, false);
                else {
                    code.body += Var(index, res_idx, result_dims).Assign(s);
                    declared[res_idx] = true;
                }
            }
        }

        virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override {
            return Array<shared_ptr<CoefficientFunction>>{cfs};
        }

        virtual void NonZeroPattern(const class ProxyUserData &ud,
                                    FlatVector<AutoDiffDiff<1, bool>> values) const override {

            Array<Vector<AutoDiffDiff<1, bool>>> vecs(cfs.Size());
            for (int i: Range(cfs)) {
                vecs[i].SetSize(cfs[i]->Dimension());
                cfs[i]->NonZeroPattern(ud, vecs[i]);
            }

            values = false;

            const auto &RI = result_idx();
            for (size_t I: Range(full_idx().TotalDim())) {
                if (!nz_all(I))
                    continue;
                const auto I_array = split(I, full_idx());
                AutoDiffDiff<1, bool> tmp(true);
                for (size_t i: Range(vecs))
                    tmp *= vecs[i](join(I_array, index_sets[i]));
                values(join(I_array, RI)) += tmp;
            }
        }

        virtual void NonZeroPattern(const class ProxyUserData &ud,
                                    FlatArray<FlatVector<AutoDiffDiff<1, bool>>> input,
                                    FlatVector<AutoDiffDiff<1, bool>> values) const override {
            values = false;

            const auto &RI = result_idx();
            for (size_t I: Range(full_idx().TotalDim())) {
                if (!nz_all(I))
                    continue;
                const auto I_array = split(I, full_idx());
                AutoDiffDiff<1, bool> tmp(true);
                for (size_t i: Range(input))
                    tmp *= input[i](join(I_array, index_sets[i]));
                values(join(I_array, RI)) += tmp;
            }
        }

        using BASE::Evaluate;

        virtual double Evaluate(const BaseMappedIntegrationPoint &ip) const override {
            if (Dimension() == 1)
                return BASE::Evaluate(ip);
            throw Exception("TensorProductCF scalar evaluate called for non-scalar result");
        }

        template<typename MIR, typename T, ORDERING ORD>
        void T_Evaluate(const MIR &ir, BareSliceMatrix<T, ORD> values) const {
            if (node) {
                node->Evaluate(ir, values);
                return;
            }

            ArrayMem<T, 1000> mem(max_mem);
            T *mem_pos = mem.Data();
            Array<FlatMatrix<T, ORD>> tmp_arrays(cfs.Size());
            for (size_t i: Range(cfs)) {
                tmp_arrays[i].AssignMemory(cfs[i]->Dimension(), ir.Size(), mem_pos);
                mem_pos += tmp_arrays[i].Height() * tmp_arrays[i].Width();
                cfs[i]->Evaluate(ir, tmp_arrays[i]);
            }

            values.AddSize(Dimension(), ir.Size()) = T(0.0);
            if (result_indices.Size()) {
                for (size_t I: Range(full_idx().TotalDim())) {
                    for (int q: Range(ir)) {
                        T tmp(1.0);
                        for (size_t i: Range(tmp_arrays))
                            tmp *= tmp_arrays[i](input_indices[i][I], q);
                        values(result_indices[I], q) += tmp;
                    }
                }
                return;
            }
            const auto &RI = result_idx();
            for (size_t I: Range(full_idx().TotalDim())) {
                const auto I_array = split(I, full_idx());
                for (int q: Range(ir)) {
                    T tmp(1.0);
                    for (size_t i: Range(tmp_arrays))
                        tmp *= tmp_arrays[i](join(I_array, index_sets[i]), q);
                    values(join(I_array, RI), q) += tmp;
                }
            }
        }

        template<typename MIR, typename T, ORDERING ORD>
        void T_Evaluate(const MIR &ir, FlatArray<BareSliceMatrix<T, ORD>> input,
                        BareSliceMatrix<T, ORD> values) const {
            if (node) {
                node->Evaluate(ir, input, values);
                return;
            }

            values.AddSize(Dimension(), ir.Size()) = T(0.0);
            if (result_indices.Size()) {
                for (size_t I: Range(full_idx().TotalDim())) {
                    for (int q: Range(ir)) {
                        T tmp(1.0);
                        for (size_t i: Range(input))
                            tmp *= input[i](input_indices[i][I], q);
                        values(result_indices[I], q) += tmp;
                    }
                }
                return;
            }

            const auto &RI = result_idx();
            for (size_t I: Range(full_idx().TotalDim())) {
                const auto I_array = split(I, full_idx());
                for (int q: Range(ir)) {
                    T tmp(1.0);
                    for (size_t i: Range(input))
                        tmp *= input[i](join(I_array, index_sets[i]), q);
                    values(join(I_array, RI), q) += tmp;
                }
            }
        }

        shared_ptr<CoefficientFunction> Diff(const CoefficientFunction *var,
                                             shared_ptr<CoefficientFunction> dir) const override {
            if (this == var)
                return dir;

            auto dres = ZeroCF(Array<int>{Dimensions()});
            for (size_t i: Range(cfs.Size())) {
                auto new_inputs = InputCoefficientFunctions();
                new_inputs[i] = cfs[i]->Diff(var, dir);
                if (new_inputs[i]->IsZeroCF())
                    continue;

                dres = dres + EinsumCF(index_signature, new_inputs, options);
            }
            return dres;
        }

        shared_ptr<CoefficientFunction> DiffJacobi(const CoefficientFunction *var) const override {
            if (this == var)
                return IdentityCF(Array<int>{Dimensions()});

            string new_symbols{new_index_symbols(index_signature, var->Dimensions().Size())};
            Array<Index> new_indices(var->Dimensions().Size());

            size_t pos = full_idx().Size();
            for (size_t i: Range(new_indices.Size()))
                new_indices[i] =
                        Index{new_symbols[i], pos, static_cast<size_t>(var->Dimensions()[i])};

            // new result index data
            MultiIndex new_result{index_sets[cfs.Size()]};
            for (size_t j: Range(new_indices.Size()))
                new_result.Append(new_indices[j]);

            Array<int> new_dims(new_result.Size());
            for (auto i: Range(new_dims.Size()))
                new_dims[i] = new_result[i].dim;

            auto dres = ZeroCF(new_dims);
            for (size_t i: Range(cfs.Size())) {
                auto new_inputs = InputCoefficientFunctions();
                new_inputs[i] = cfs[i]->DiffJacobi(var);
                if (new_inputs[i]->IsZeroCF())
                    continue;

                Array<MultiIndex> new_index_sets{index_sets};
                for (size_t j: Range(new_indices.Size()))
                    new_index_sets[i].Append(new_indices[j]);
                new_index_sets[cfs.Size()] = new_result;
                string new_index_signature{
                        form_index_signature(new_index_sets.Range(0, new_index_sets.Size() - 1))};

                dres = dres + EinsumCF(new_index_signature, new_inputs, options);
            }
            return dres;
        }
    };

    pair<Array<MultiIndex>, Array<shared_ptr<CoefficientFunction>>>
    optimize_special_inputs(Array<MultiIndex> index_sets, Array<shared_ptr<CoefficientFunction>> cfs,
                            const map<string, bool> &options) {

        const auto identity_descr = IdentityCF(1)->GetDescription();

        Array<shared_ptr<CoefficientFunction>> new_cfs;
        new_cfs.SetAllocSize(cfs.Size());
        Array<MultiIndex> new_index_sets;
        new_index_sets.SetAllocSize(cfs.Size() + 2);


        if (get_option(options, "expand_einsum", true)) {
            cout << "TP: expand einsum CFs (no recursion)" << endl;
            string signature = form_index_signature(index_sets);

            // split signature into sub strings
            vector<string> original_input_strings;
            original_input_strings.push_back("");
            string sub{};
            string result_sub{};
            int mode = 0;
            for (auto item: signature) {
                if (item == '\0')
                    break;
                if (item == '-')
                    continue;
                if (item == '>') {
                    mode++;
                    continue;
                }
                if (mode == 0) {
                    if (item == ',') {
                        original_input_strings.push_back(sub);
                        sub = "";
                    } else
                        sub += item;
                } else
                    result_sub += item;
            }

            string new_signature{};
            char new_symbol = 'A';

            for (size_t i: Range(cfs)) {
                if (!new_signature.empty())
                    new_signature += ",";

                if (auto cfi = dynamic_pointer_cast<EinsumCoefficientFunction>(cfs[i]); cfi) {
                    // expand the einsum object

                    // append inputs
                    auto input_cfs = cfi->InputCoefficientFunctions();
                    for (auto j: Range(input_cfs)) {
                        new_cfs.Append(input_cfs[i]);
                    }

                    // gather data
                    auto input_signature = cfi->IndexSignature();
                    auto input_multi_indices = compute_multi_indices(input_signature, input_cfs);
                    auto input_result_multi_index = input_multi_indices[input_cfs.Size()];
                    auto input_all_multi_index = input_multi_indices[input_cfs.Size() + 1];

                    set<char> input_all_symbols{};
                    for (auto j: Range(input_all_multi_index))
                        input_all_symbols.insert(input_all_multi_index[j].symbol);

                    set<char> input_result_symbols{};
                    for (auto j: Range(input_result_multi_index))
                        input_result_symbols.insert(input_result_multi_index[j].symbol);

                    /// create a mapping for substitution of already used dummy indices
                    map<char, char> index_map;

                    for (auto item: input_all_symbols) {
                        // already found another symbol
                        if (index_map.count(item) > 0)
                            continue;

                        // no need to find another symbol
                        if (input_result_symbols.count(item) > 0 ||
                            (signature.find(item) == std::string::npos &&
                             new_signature.find(new_symbol) == std::string::npos))
                            continue;

                        // as long as symbol candidate is already in use...
                        while (signature.find(new_symbol) != std::string::npos ||
                               new_signature.find(new_symbol) != std::string::npos) {
                            if (new_symbol == 'Z')
                                new_symbol = 'a';
                            else if (new_symbol == 'z')
                                throw NG_EXCEPTION("out of index symbols");
                            else
                                new_symbol++;
                        }
                        index_map[item] = new_symbol;
                    }

                    // add new signature part
                    for (auto item: input_signature) {
                        if (item == '-')
                            break;
                        if (index_map.count(item) > 0)
                            new_signature.push_back(index_map[item]);
                        else
                            new_signature.push_back(item);
                    }
                } else {
                    new_cfs.Append(cfs[i]);
                    new_signature += original_input_strings[i];
                }
            }
            // put everyting together
            new_signature += string("->") + result_sub;
            index_sets = compute_multi_indices(new_signature, new_cfs);
            cfs = new_cfs;
        }

        if (get_option(options, "expand_higher_order_identity", true) ||
            get_option(options, "eliminate_identity_tensors", false)) {
            cout << "TP: expand higher order identity tensors" << endl;

            for (size_t i: Range(cfs))
                if (cfs[i]->GetDescription() == identity_descr) {
                    auto dims = cfs[i]->Dimensions();
                    for (size_t j: Range(dims.Size() / 2)) {
                        MultiIndex id_mi;
                        id_mi.Append(index_sets[i][j]);
                        id_mi.Append(index_sets[i][j + dims.Size() / 2]);
                        new_cfs.Append(IdentityCF(Array<int>{dims.Range(j, j + 1)}));
                        new_index_sets.Append(id_mi);
                    }
                } else {
                    new_cfs.Append(cfs[i]);
                    new_index_sets.Append(index_sets[i]);
                }
            new_index_sets.Append(index_sets.Part(cfs.Size()));
            index_sets = new_index_sets;
            cfs = new_cfs;
        }

        if (get_option(options, "eliminate_identity_tensors", false)) {
            cout << "TP: trying to eliminate identities" << endl;

            // detect identity tensors that can be removed
            Array<bool> remove(index_sets.Size());
            remove = false;
            for (size_t i: Range(cfs)) {
                if (cfs[i]->GetDescription() == identity_descr)
                    if (substitute_id_index(index_sets, new_index_sets[i][0], new_index_sets[i][1],
                                            i, remove, true))
                        remove[i] = true;
            }

            // remove corresponding index sets
            new_index_sets.DeleteAll();
            for (size_t i: Range(index_sets))
                if (!remove[i])
                    new_index_sets.Append(index_sets[i]);

            // remove corresponding cfs
            new_cfs.DeleteAll();
            for (size_t i: Range(cfs))
                if (!remove[i])
                    new_cfs.Append(cfs[i]);
            index_sets = new_index_sets;
            cfs = new_cfs;
        }

        return {move(index_sets), move(cfs)};
    }

    shared_ptr<CoefficientFunction> EinsumCF(const string &index_signature,
                                             const Array<shared_ptr<CoefficientFunction>> &cfs,
                                             const map<string, bool> &options) {
        return make_shared<EinsumCoefficientFunction>(index_signature, cfs, options);
    }

    namespace tensor_internal {
        bool is_tensor_product(shared_ptr<CoefficientFunction> cf) {
            return dynamic_pointer_cast<EinsumCoefficientFunction>(cf) != nullptr;
        }

        string tensor_product_signature(shared_ptr<CoefficientFunction> cf) {
            if (auto tp = dynamic_pointer_cast<EinsumCoefficientFunction>(cf); tp)
                return tp->IndexSignature();
            throw NG_EXCEPTION("given cf is not a EinsumCF");
        }
    }

} // namespace ngfem
