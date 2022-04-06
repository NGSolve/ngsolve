/*********************************************************************/
/* File:   tensorcoefficient.cpp                                     */
/* Author: Matthias Rambausek                                        */
/* Date:   01. Apr. 2022                                             */
/*********************************************************************/

#include <cmath>
#include "../ngstd/ngstd.hpp"
#include "../ngstd/python_ngstd.hpp"

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
                    throw NG_EXCEPTION("Did not find any unused symbol in A-Z and a-z. Consider using 'optimize_path' and disabling 'einsum_expansion' at some point.");
                if (find(begin(existing), end(existing), new_char) != end(existing))
                    ++new_char;
                else
                    existing += new_char;
            }

            return string{existing.end() - nnew, existing.end()};
        }

        string form_index_signature(const vector<string>& parts)
        {
          stringstream signature{};
          // note: last of parts refers to result
          for (auto i : Range(parts.size() - 1))
            signature << (i == 0 ? "" : ",") << parts[i];

          signature << "->" << parts.back();
          return signature.str();
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

        Array<int> index_dimensions(const MultiIndex &mi)
        {
          Array<int> res_dims(mi.Size());
          for (size_t i: Range(res_dims.Size()))
            res_dims[i] = mi[i].dim;
          return res_dims;
        }

        string sanitize_signature(string signature)
        {
          if (signature.find("->") == string::npos)
          {
            map<char, size_t> symbol_count{};
            for (auto c: signature)
            {
              if (c == ',' || c == '-' || c == '>' || c == '\0')
                continue;
              if (symbol_count.count(c) == 0)
                symbol_count[c] = 0;
              symbol_count[c]++;
            }
            for (auto item : symbol_count)
              if (item.second != 2)
                throw NG_EXCEPTION("Signature does not contain '->'. ATM, this "
                                   "is only allowed if each index appears "
                                   "exactly twice.");
            signature += "->";
          }
          return move(signature);
        }

        vector<string> split_signature(string signature)
        {
          signature = sanitize_signature(signature);
          vector<string> parts{};
          stringstream new_part{};
          for (auto item: signature)
            {
              if (item == ',' || item == '-')
                {
                  parts.push_back(new_part.str());
                  new_part = stringstream{};
                }
              else if (item == '>' || item == '\0')
                continue;
              else
                new_part << item;
            }
          parts.push_back(new_part.str());
          return move(parts);
        }

        optional<string> substitute_id_index(string signature, pair<char, char> id_indices,
                                 size_t id_pos, const FlatArray<bool> marked_for_removal,
                                 bool recurse)
        {
            bool subs{false};
            auto parts = split_signature(signature);
            auto &result_part = parts.back();

            if (result_part.find(id_indices.first) != string::npos) {
                if (recurse)
                  return substitute_id_index(signature,
                                             {id_indices.second, id_indices.first},
                                             id_pos,
                                             marked_for_removal);
                else
                  return {};
            }

            for (size_t i2: Range(parts.size() - 1)) {
                // do not touch current id tensor pos
                if (i2 == id_pos || marked_for_removal[i2])
                    continue;

                auto &idx2 = parts[i2];
                for (char& c : idx2) {
                    // substitute
                    if (c == id_indices.first) {
                        c = id_indices.second;
                        subs = true;
                    }
                }
            }

            if (subs)
              return form_index_signature(parts);
            else if (recurse)
              return substitute_id_index(signature,
                                         {id_indices.second, id_indices.first},
                                         id_pos, marked_for_removal);
            else
              return {};
        }

        string
        expand_einsum_part(string target_signature_part, const string& nested_signature, const string& used_symbols)
        {
          auto parts_nested = split_signature(nested_signature);
          auto result_nested = parts_nested.back();

          if (result_nested.size() != target_signature_part.size())
            throw NG_EXCEPTION("Incompatible multi-indices");

          map<char, char> symbol_map{{',', ','}};
          for (auto i : Range(target_signature_part.size()))
            if (result_nested[i] != target_signature_part[i])
              symbol_map[result_nested[i]] = target_signature_part[i];
            else
              symbol_map[result_nested[i]] = result_nested[i];

          stringstream relabel{};
          for (const auto& part : parts_nested)
            for (auto item : part)
              {
                if (symbol_map.find(item) != symbol_map.end())
                  continue;

                if (used_symbols.find(item) != string::npos)
                  relabel << item;
                else
                  symbol_map[item] = item;
              }

          string relabel_str = relabel.str();
          string new_symbols = new_index_symbols(used_symbols, relabel_str.size());

          for (auto i : Range(relabel_str.size()))
            symbol_map[relabel_str[i]] = new_symbols[i];

          stringstream new_part{};
          for (auto c : nested_signature)
            {
              if (c == '-')
                break;

              new_part << symbol_map[c];
            }

          string new_part_str = new_part.str();
          return new_part_str;
        }

        pair<string, Array<shared_ptr<CoefficientFunction>>>
        flatten_einsum(string signature,
                       const Array<shared_ptr<CoefficientFunction>>& cfs,
                       const map<string, bool> &options) {

          const auto identity_descr = IdentityCF(1)->GetDescription();

          Array<shared_ptr<CoefficientFunction>> new_cfs;
          new_cfs.SetAllocSize(cfs.Size());

          cout << "TP: flatten einsum CF (no recursion)" << endl;
          auto parts = split_signature(signature);

          string used_symbols = signature;

          for (auto i : Range(cfs))
            if (auto cfi = dynamic_pointer_cast<EinsumCoefficientFunction>(cfs[i]); cfi)
              {
                auto nested_signature = cfi->ExpandedIndexSignature();
                parts[i] = expand_einsum_part(parts[i], nested_signature, used_symbols);
                used_symbols += parts[i];
                auto nested_inputs = cfi->ExpandedInputCoefficientFunctions();
                new_cfs.Append(nested_inputs);
              }
            else
              new_cfs.Append(cfs[i]);

          return {move(form_index_signature(parts)), move(new_cfs)};
        }

        pair<string, Array<shared_ptr<CoefficientFunction>>>
        expand_higher_order_identities(string signature,
                       const Array<shared_ptr<CoefficientFunction>>& cfs,
                       [[maybe_unused]] const map<string, bool> &options) {

          const auto identity_descr = IdentityCF(1)->GetDescription();

          Array<shared_ptr<CoefficientFunction>> new_cfs;
          new_cfs.SetAllocSize(cfs.Size());

          cout << "TP: expand higher-order identities" << endl;
          auto parts = split_signature(signature);

          for (auto i : Range(cfs))
            if (cfs[i]->GetDescription() == identity_descr)
            {
              auto dims = cfs[i]->Dimensions();
              stringstream new_part{};
              for (size_t j: Range(dims.Size() / 2)) {
                new_cfs.Append(IdentityCF(dims[j]));
                new_part << (j > 0 ? "," : "")
                         << parts[i][j] << parts[i][j + dims.Size() / 2];
              }
              parts[i] = new_part.str();
            }
            else
              new_cfs.Append(cfs[i]);

          return {move(form_index_signature(parts)), move(new_cfs)};
        }


        pair<string, Array<shared_ptr<CoefficientFunction>>>
        optimize_identities(string signature,
                            const Array<shared_ptr<CoefficientFunction>>& cfs,
                            [[maybe_unused]] const map<string, bool>& options)
        {
          const auto identity_descr = IdentityCF(1)->GetDescription();

          Array<shared_ptr<CoefficientFunction>> new_cfs;
          new_cfs.SetAllocSize(cfs.Size());
          auto parts = split_signature(signature);
          map<size_t, shared_ptr<CoefficientFunction>> cf_subs{};

          // detect identity tensors that can be removed
          Array<bool> remove(cfs.Size());
          remove = false;
          for (size_t i: Range(cfs))
          {
            if (cfs[i]->GetDescription() == identity_descr &&
                cfs[i]->Dimensions().Size() == 2)
              if (parts[i][0] == parts[i][1] &&
                  parts.back().find(parts[i][0]) == string::npos)
              {
                parts[i] = parts[i][0];
                cf_subs[i] = ConstantCF(cfs[i]->Dimensions()[0]);
                cf_subs[i]->SetDimensions(Array<int>{1});
              }
              else if (auto new_signature =
                      substitute_id_index(signature,
                                          {parts[i][0], parts[i][1]}, i,
                                          remove, true); new_signature)
              {
                signature = new_signature.value();
                parts = split_signature(signature);
                remove[i] = true;
              }
          }

          decltype(parts) new_parts{};
          // remove corresponding index sets
          for (auto i: Range(remove))
            if (!remove[i])
              new_parts.push_back(parts[i]);

          // result...
          new_parts.push_back(parts.back());

          // remove corresponding cfs
          for (size_t i: Range(cfs))
            if (cf_subs.count(i))
              new_cfs.Append(cf_subs[i]);
            else if (!remove[i])
              new_cfs.Append(cfs[i]);

          return {move(form_index_signature(new_parts)), move(new_cfs)};
        }


        shared_ptr<CoefficientFunction>
        optimize_legacy(const string& signature,
                      const Array<shared_ptr<CoefficientFunction>>& cfs,
                      [[maybe_unused]] const map<string, bool> &options) {

          const auto index_sets = compute_multi_indices(signature, cfs);
          for (size_t i: Range(cfs))
            if (cfs[i]->IsZeroCF())
            {
              auto dims = index_dimensions(index_sets[cfs.Size()]);
              return ZeroCF(dims);
            }

          const auto identity_descr = IdentityCF(1)->GetDescription();
          cout << "TP: trying to detect some 'legacy' operations" << endl;

          const bool optimize_identities =
              get_option(options, "optimize_identities", true);

          // Trace
          if (cfs.Size() == 1 && index_sets[0].Size() == 2) {
            if (index_sets[0][0].symbol == index_sets[0][1].symbol && index_sets[1].Size() == 0)
            {
              cout << "EinsumCF: detected trace" << endl;
              return TraceCF(cfs[0]);
            }
          }

          // Transpose
          if (cfs.Size() == 1)
            return optimize_transpose(signature, cfs, options);

          // Mat * Vec
          if (cfs.Size() == 2 && index_sets[0].Size() == 2 && index_sets[1].Size() == 1 &&
              index_sets[2].Size() == 1 && index_sets[0][0].symbol != index_sets[0][1].symbol) {
            if (index_sets[1][0].symbol == index_sets[0][1].symbol)
            {
              // NOTE: no other way to detect identity!
              if (cfs[0]->GetDescription() == identity_descr &&
                  optimize_identities)
              {
                cout << "EinsumCF: detected I * vec" << endl;
                return cfs[1];
              }
              cout << "EinsumCF: detected Mat * Vec" << endl;
              return cfs[0] * cfs[1];
            }
            else if (index_sets[1][0].symbol == index_sets[0][0].symbol)
            {
              // NOTE: no other way to detect identity!
              if (cfs[0]->GetDescription() == identity_descr &&
                  optimize_identities)
              {
                cout << "EinsumCF: detected I * vec" << endl;
                return cfs[1];
              }
              cout << "EinsumCF: detected Mat.trans * Vec" << endl;
              return TransposeCF(cfs[0]) * cfs[1];
            }
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

              if (cfs[i0]->GetDescription() == identity_descr &&
                  optimize_identities)
              {
                cout << "EinsumCF: detected I * Mat" << endl;
                return cfs[i1];
              }

              if (cfs[i1]->GetDescription() == identity_descr &&
                  optimize_identities)
              {
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
            if (cfs[0]->GetDescription() == identity_descr &&
                optimize_identities) {
              cout << "EinsumCF: detected I * Mat" << endl;
              return cfs[1];
            }

            if (cfs[1]->GetDescription() == identity_descr &&
                optimize_identities) {
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
            if (cfs[0]->GetDescription() == identity_descr &&
                optimize_identities) {
              cout << "EinsumCF: detected I * Mat.trans" << endl;
              return TransposeCF(cfs[1]);
            }

            if (cfs[1]->GetDescription() == identity_descr &&
                optimize_identities) {
              cout << "EinsumCF: detected Mat * I" << endl;
              return cfs[0];
            }

            cout << "EinsumCF: detected Mat * Mat.trans" << endl;
            return cfs[0] * TransposeCF(cfs[1]);
          }
          return {};
        }

        shared_ptr<CoefficientFunction>
        optimize_transpose(const string& signature,
                           const Array<shared_ptr<CoefficientFunction>>& cfs,
                           [[maybe_unused]] const map<string, bool> &aoptions)
        {
          const auto parts = split_signature(signature);
          if (set<char>{parts[0].begin(), parts[0].end()}.size() != parts[0].size() ||
              set<char>{parts[1].begin(), parts[1].end()}.size() != parts[1].size())
            return {};
          else if (parts[0] == parts[1])
            return cfs[0];
          else if (parts[0].size() == 2)
          {
            cout << "EinsumCF: detected transpose" << endl;
            return TransposeCF(cfs[0]);
          }
          else
          {
            cout << "EinsumCF: detected tensor transpose" << endl;
            Array<int> ordering;
            ordering.SetSize(cfs[0]->Dimensions().Size());
            for (auto i : Range(ordering))
              ordering[i] = parts[1].find(parts[0][i]);
            return MakeTensorTransposeCoefficientFunction(cfs[0],
                                                          move(ordering));
          }
        }

        shared_ptr<CoefficientFunction>
        optimize_path(const string &signature,
                      const Array<shared_ptr<CoefficientFunction>>& input_cfs,
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
            tp_inputs.Append(EinsumCF(new_signature, new_inputs, options));
          }
          return tp_inputs.Last();
        }

        LeviCivitaCoefficientFunction::LeviCivitaCoefficientFunction(int adim) : BASE(1), dim{adim}
        {
          Array<int> dims(dim);
          char symbol = 'a';
          for (size_t i : Range(dim)) {
            dims[i] = dim;
            mi.Append(Index{symbol++, i, static_cast<size_t>(dim)});
          }

          SetDimensions(dims.Part(0));
        }


        void LeviCivitaCoefficientFunction::GenerateCode(Code &code, FlatArray<int> inputs, int index,
                                    bool skip_zeroes) const
        {
          LocalHeap lh(100000, "Levi-Cevita code gen");

          auto dims = Dimensions();

          auto nzvec =
              nonzero_pattern(const_cast<LeviCivitaCoefficientFunction *>(this)
                                  ->shared_from_this());

          for (size_t I : Range(Dimension())) {

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

        void LeviCivitaCoefficientFunction::NonZeroPattern(const class ProxyUserData &ud,
                                      FlatVector<AutoDiffDiff<1, bool>> values) const
        {
          for (size_t I : Range(Dimension())) {
            const auto I_array = split(I, mi);

            if (is_even_iota_permutation(I_array.begin(), I_array.end()))
              values(I) = true;
            else if (is_odd_iota_permutation(I_array.begin(), I_array.end()))
              values(I) = true;
            else
              values(I) = false;
          }
        }

        void LeviCivitaCoefficientFunction::NonZeroPattern(const class ProxyUserData &ud,
                                      FlatArray<FlatVector<AutoDiffDiff<1, bool>>> input,
                                      FlatVector<AutoDiffDiff<1, bool>> values) const
        {
          NonZeroPattern(ud, values);
        }


        double
        LeviCivitaCoefficientFunction::Evaluate(
            const BaseMappedIntegrationPoint &ip) const
        {
          if (Dimension() == 1)
            return BASE::Evaluate(ip);
          throw Exception(
              "LeviCivitaCF scalar evaluate called for non-scalar result");
        }


        shared_ptr<CoefficientFunction>
        LeviCivitaCoefficientFunction::Diff(
            const CoefficientFunction *var,
            shared_ptr<CoefficientFunction> dir) const
        {
          if (this == var)
            return dir;

          return ZeroCF(Dimensions());
        }

        shared_ptr<CoefficientFunction>
        LeviCivitaCoefficientFunction::DiffJacobi(
            const CoefficientFunction *var) const
        {
          if (this == var)
            return IdentityCF(Dimensions());

          Array<int> dims{Dimensions()};
          dims.Append(Dimensions());
          return ZeroCF(dims);
        }


        EinsumCoefficientFunction::EinsumCoefficientFunction(
            const string &aindex_signature,
            const Array<shared_ptr<CoefficientFunction>> &acfs,
            const map<string, bool> &aoptions)
            :
              BASE(1,
                   std::find_if(acfs.begin(), acfs.end(),
                                [](const auto &item)
                                {
                                  return item->IsComplex();
                                }) != acfs.end()),
              original_inputs{acfs},
              original_index_signature{aindex_signature},
              max_mem{0},
              node{},
              options{aoptions} {

          if (get_option(options, "expand_einsum", true))
          {
            try
            {
              tie(expanded_index_signature, expanded_inputs) =
                  flatten_einsum(original_index_signature, original_inputs, options);
            }
            catch (const Exception& e) {
              cout << "Caught exception during Einsum flattening:\n"
                   << e.What() << endl;
              tie(expanded_index_signature, expanded_inputs) =
                  tie(original_index_signature, original_inputs);
              options["expand_einsum"] = false;
            }
          }
          else
            tie(expanded_index_signature, expanded_inputs) =
                tie(original_index_signature, original_inputs);

          if (get_option(options, "optimize_path", false))
          {
            if (get_option(options, "optimize_identities", false))
            {
              tie(index_signature, cfs) = expand_higher_order_identities(
                  expanded_index_signature, expanded_inputs, options);
              tie(index_signature, cfs) = optimize_identities(
                  index_signature, cfs, options);
              node = optimize_path(index_signature, cfs, options);
            }
            else
              node = optimize_path(expanded_index_signature, expanded_inputs, options);
          }
          else if (get_option(options, "optimize_identities", false))
          {
            tie(index_signature, cfs) = expand_higher_order_identities(original_index_signature, original_inputs, options);
            tie(index_signature, cfs) = optimize_identities(index_signature, cfs, options);
          }
          else
            tie (index_signature, cfs) = {original_index_signature, original_inputs};

          if (!node && cfs.Size() < 3 && get_option(options, "use_legacy_ops", false))
            node = optimize_legacy(index_signature, cfs, options);


          if (node)
            SetDimensions(node->Dimensions());
          else
          {
            // compute index mappings and nonzero patterns
            const auto index_sets = compute_multi_indices(index_signature, cfs);
            const auto& RI = index_sets[cfs.Size()];
            const auto& FI = index_sets[cfs.Size() + 1];

            if (RI.Size() > 0)
            {
              const auto dims = index_dimensions(RI);
              SetDimensions(dims);
            }

            for (size_t i: Range(cfs))
              max_mem += index_sets[i].TotalDim();

            index_maps = build_index_maps(index_sets, {});

            nz_inputs.SetSize(cfs.Size());
            for (size_t i: Range(cfs))
              nz_inputs[i] = nonzero_pattern(cfs[i]);

            nz_all.SetSize(index_maps.Height());
            nz_all = true;
            nz_result = nonzero_pattern(this);
            const auto cres = cfs.Size();

            for (size_t I: Range(index_maps.Height()))
            {
              const auto I_map = index_maps.Row(I);
              if (!nz_result(I_map(cres)))
              {
                nz_all(I) = false;
                continue;
              }
              for (size_t i: Range(cfs))
                if (!nz_inputs[i](I_map(i)))
                {
                  nz_all(I) = false;
                  continue;
                }
            }

            if (get_option(options, "sparse_evaluation", true))
              sparse_index_maps = build_index_maps(index_sets, nz_all);

          }
        }

        Matrix<int> EinsumCoefficientFunction::build_index_maps(
                const Array<MultiIndex>& index_sets, const optional<Vector<bool>>& nz_pattern)
        {
          Matrix<int> imaps;

          const auto& RI = index_sets[cfs.Size()];
          const auto& FI = index_sets[cfs.Size() + 1];
          const auto cres = cfs.Size();

          if (!nz_pattern)
          {
              imaps.SetSize(FI.TotalDim(), cfs.Size() + 1);
              for (size_t I: Range(imaps.Height()))
              {
                  const auto I_array = split(I, FI);
                  imaps(I, cres) = join(I_array, RI);
                  for (size_t i: Range(cfs))
                      imaps(I, i) = join(I_array, index_sets[i]);
              }
          }
          else
          {
              auto nnz = count(nz_pattern->Data(), nz_pattern->Data() + nz_pattern->Size(), true);
              imaps.SetSize(nnz, cfs.Size() + 1);
              size_t nzi = 0;
              for (size_t I: Range(FI.TotalDim()))
              {
                  if (!(*nz_pattern)[I])
                      continue;

                  const auto I_array = split(I, FI);
                  imaps(nzi, cres) = join(I_array, RI);
                  for (size_t i: Range(cfs))
                      imaps(nzi, i) = join(I_array, index_sets[i]);
                  ++nzi;
              }
          }
          return imaps;
        }

        string EinsumCoefficientFunction::GetDescription() const
        {
          stringstream descr{};
          descr << "EinsumCF " << original_index_signature;

          if (node)
            descr << " with optimized node " << node->GetDescription();

          return descr.str();
        }

        void EinsumCoefficientFunction::GenerateCode(
            Code &code, FlatArray<int> inputs, int index,
            bool skip_zeroes) const
        {
          if (node)
          {
            node->GenerateCode(code, inputs, index);
            return;
          }

          LocalHeap lh(100000, "TP code gen");

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
            const auto dims_i = cfs[i]->Dimensions();
            input_dims[i].Assign(dims_i.Size(), lh);
            input_dims[i] = dims_i;
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

          const FlatArray<int> result_dims = Dimensions();

          FlatArray<bool> declared(index_maps.Height(), lh);
          declared = false;

          const auto cres = cfs.Size();
          for (size_t I: Range(declared)) {
            const auto I_map = index_maps.Row(I);
            const auto res_idx = I_map(cres);

            if (skip_zeroes && !result_nz_vec[res_idx])
              continue;

            CodeExpr s;
            if (!result_nz_vec[res_idx] && !declared[res_idx])
              s = Var(0.0);
            else
              for (size_t i: Range(inputs)) {
                const auto I_local = I_map(i);
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

        void EinsumCoefficientFunction::NonZeroPattern(
            const class ProxyUserData &ud,
            FlatVector<AutoDiffDiff<1, bool>> values) const
        {
          if (node)
          {
            node->NonZeroPattern(ud, values);
            return;
          }

          Array<Vector<AutoDiffDiff<1, bool>>> vecs(cfs.Size());
          for (int i: Range(cfs)) {
            vecs[i].SetSize(cfs[i]->Dimension());
            cfs[i]->NonZeroPattern(ud, vecs[i]);
          }

          values = false;
          const auto cres = cfs.Size();
          for (size_t I: Range(index_maps.Height())) {
            if (!nz_all(I))
              continue;
            const auto& I_map = index_maps.Row(I);
            AutoDiffDiff<1, bool> tmp(true);
            for (size_t i: Range(vecs))
              tmp *= vecs[i](I_map(i));
            values(I_map(cres)) += tmp;
          }
        }

        void EinsumCoefficientFunction::NonZeroPattern(
            const class ProxyUserData &ud,
            FlatArray<FlatVector<AutoDiffDiff<1, bool>>> input,
            FlatVector<AutoDiffDiff<1, bool>> values) const
        {
          if (node)
          {
            node->NonZeroPattern(ud, input, values);
            return;
          }

          values = false;
          const auto cres = cfs.Size();
          for (size_t I: Range(index_maps.Height())) {
            if (!nz_all(I))
              continue;
            const auto& I_map = index_maps.Row(I);
            AutoDiffDiff<1, bool> tmp(true);
            for (size_t i: Range(input))
              tmp *= input[i](I_map(i));
            values(I_map(cres)) += tmp;
          }
        }

        shared_ptr<CoefficientFunction>
        EinsumCoefficientFunction::Diff(
            const CoefficientFunction *var,
            shared_ptr<CoefficientFunction> dir) const
        {
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
          // TODO: great potential for optimization when equivalent objects are
          //  identified in Compile
          return dres;
        }

        shared_ptr<CoefficientFunction>
        EinsumCoefficientFunction::DiffJacobi(
            const CoefficientFunction *var) const
        {
          if (this == var)
            return IdentityCF(Dimensions());

          Array<int> res_dims;
          res_dims.Append(Dimensions());
          res_dims.Append(var->Dimensions());

          auto dres = ZeroCF(res_dims);
          string new_symbols{};
          try
          {
            new_symbols = new_index_symbols(original_index_signature,
                                            var->Dimensions().Size());
          }
          catch (const Exception& e)
          {
            if (!options.at("optimize_path"))
            {
              cout << "Caught exception during DiffJacobi:\n"
                   << e.What()
                   << "\n"
                   << "Trying again with a broken-down EinsumCF." << endl;
              auto opts = options;
              opts["optimize_path"] = true;
              opts["expand_einsum"] = false;
              return Optimize(options)->DiffJacobi(var);
            }
            throw e;
          }

          auto parts = split_signature(original_index_signature);

          for (size_t i: Range(cfs.Size())) {
            auto new_inputs{original_inputs};
            new_inputs[i] = cfs[i]->DiffJacobi(var);
            if (new_inputs[i]->IsZeroCF())
              continue;
            auto new_parts{parts};
            new_parts[i] = parts[i] + new_symbols;
            new_parts.back() += new_symbols;
            dres = dres + EinsumCF(form_index_signature(new_parts), new_inputs, options);
          }

          // TODO: great potential for optimization when equivalent objects are
          //  identified in Compile
          return dres;
        }

        shared_ptr<EinsumCoefficientFunction>
        EinsumCoefficientFunction::Optimize(const map<string, bool> &aoptions) const
        {
          return make_shared<EinsumCoefficientFunction>(
              original_index_signature, original_inputs, aoptions
              );
        }
    }

    using namespace tensor_internal;


    shared_ptr<CoefficientFunction> LeviCivitaCF(int dimension) {
        return make_shared<LeviCivitaCoefficientFunction>(dimension);
    }

    shared_ptr<CoefficientFunction> EinsumCF(const string &index_signature,
                                             const Array<shared_ptr<CoefficientFunction>> &cfs,
                                             const map<string, bool> &options) {
        return make_shared<EinsumCoefficientFunction>(index_signature, cfs, options);
    }

} // namespace ngfem
