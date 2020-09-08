#ifndef NGSOLVE_PYTHON_LINALG_HPP
#define NGSOLVE_PYTHON_LINALG_HPP

#include <la.hpp>

// Workaround to ensure same lifetime for C++ and Python objects
// see https://github.com/pybind/pybind11/issues/1546
namespace pybind11::detail {
  using ngla::BaseMatrix;
  template<>
  struct type_caster<std::shared_ptr<BaseMatrix>>
  {
    PYBIND11_TYPE_CASTER (std::shared_ptr<BaseMatrix>, _("BaseMatrix"));

    using BaseCaster = copyable_holder_caster<BaseMatrix, std::shared_ptr<BaseMatrix>>;

    bool load (pybind11::handle src, bool b)
    {
      BaseCaster bc;
      bool success = bc.load (src, b);
      if (!success)
      {
        return false;
      }

      auto py_obj = reinterpret_borrow<object>(src);
      auto base_ptr = static_cast<std::shared_ptr<BaseMatrix>> (bc);

      // Construct a shared_ptr to the py::object
      auto py_obj_ptr = std::shared_ptr<object>{
          new object{py_obj},
          [](object * py_object_ptr) {
              // It's possible that when the shared_ptr dies we won't have the
              // gil (if the last holder is in a non-Python thread), so we
              // make sure to acquire it in the deleter.
              gil_scoped_acquire gil;
              delete py_object_ptr;
         }
      };

      value = std::shared_ptr<BaseMatrix> (py_obj_ptr, base_ptr.get ());
      return true;
    }

    static handle cast (std::shared_ptr<BaseMatrix> base,
                        return_value_policy rvp,
                        handle h)
    {
      return BaseCaster::cast (base, rvp, h);
    }
  };

  template <>
  struct is_holder_type<BaseMatrix, std::shared_ptr<BaseMatrix>> : std::true_type {};
}

#endif // NGSOLVE_PYTHON_LINALG_HPP

