#ifndef DEFS_HPP
#define DEFS_HPP

#include "openPMD/openPMD.hpp"

#include <jlcxx/array.hpp>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/module.hpp>
#include <jlcxx/stl.hpp>
#include <jlcxx/tuple.hpp>
#include <jlcxx/type_conversion.hpp>

#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using namespace openPMD;

template <std::size_t I> struct sized_uint;
template <> struct sized_uint<1> { using type = std::uint8_t; };
template <> struct sized_uint<2> { using type = std::uint16_t; };
template <> struct sized_uint<4> { using type = std::uint32_t; };
template <> struct sized_uint<8> { using type = std::uint64_t; };
template <std::size_t I> using sized_uint_t = typename sized_uint<I>::type;

#warning "TODO7"
// template <typename T> using array7 = std::array<T, 7>;

/*
 * Generate code fo all openPMD types. Use is e.g. as follows:
 *   #define USE_TYPE(NAME, ENUM, TYPE)                    \
 *     type.method("get_" NAME, &Attribute::get<TYPE>);
 *   FORALL_OPENPMD_TYPES(USE_TYPE)
 *   #undef USE_TYPE
 */

// We disable `long double` since Julia does not support this type
#define FORALL_OPENPMD_TYPES(MACRO)                                            \
  MACRO("CHAR", Datatype::CHAR, char)                                          \
  MACRO("UCHAR", Datatype::UCHAR, unsigned char)                               \
  MACRO("SHORT", Datatype::SHORT, short)                                       \
  MACRO("INT", Datatype::INT, int)                                             \
  MACRO("LONG", Datatype::LONG, long)                                          \
  MACRO("LONGLONG", Datatype::LONGLONG, long long)                             \
  MACRO("USHORT", Datatype::USHORT, unsigned short)                            \
  MACRO("UINT", Datatype::UINT, unsigned int)                                  \
  MACRO("ULONG", Datatype::ULONG, unsigned long)                               \
  MACRO("ULONGLONG", Datatype::ULONGLONG, unsigned long long)                  \
  MACRO("FLOAT", Datatype::FLOAT, float)                                       \
  MACRO("DOUBLE", Datatype::DOUBLE, double)                                    \
  /* MACRO("LONG_DOUBLE", Datatype::LONG_DOUBLE, long double) */               \
  MACRO("CFLOAT", Datatype::CFLOAT, std::complex<float>)                       \
  MACRO("CDOUBLE", Datatype::CDOUBLE, std::complex<double>)                    \
  /* MACRO("CLONG_DOUBLE", Datatype::CLONG_DOUBLE, std::complex<long           \
   * double>) */                                                               \
  MACRO("STRING", Datatype::STRING, std::string)                               \
  MACRO("VEC_CHAR", Datatype::VEC_CHAR, std::vector<char>)                     \
  MACRO("VEC_UCHAR", Datatype::VEC_UCHAR, std::vector<unsigned char>)          \
  MACRO("VEC_SHORT", Datatype::VEC_SHORT, std::vector<short>)                  \
  MACRO("VEC_INT", Datatype::VEC_INT, std::vector<int>)                        \
  MACRO("VEC_LONG", Datatype::VEC_LONG, std::vector<long>)                     \
  MACRO("VEC_LONGLONG", Datatype::VEC_LONGLONG, std::vector<long long>)        \
  MACRO("VEC_USHORT", Datatype::VEC_USHORT, std::vector<unsigned short>)       \
  MACRO("VEC_UINT", Datatype::VEC_UINT, std::vector<unsigned int>)             \
  MACRO("VEC_ULONG", Datatype::VEC_ULONG, std::vector<unsigned long>)          \
  MACRO("VEC_ULONGLONG", Datatype::VEC_ULONGLONG,                              \
        std::vector<unsigned long long>)                                       \
  MACRO("VEC_FLOAT", Datatype::VEC_FLOAT, std::vector<float>)                  \
  MACRO("VEC_DOUBLE", Datatype::VEC_DOUBLE, std::vector<double>)               \
  /* MACRO("VEC_LONG_DOUBLE", Datatype::VEC_LONG_DOUBLE, std::vector<long      \
   * double>) */                                                               \
  MACRO("VEC_CFLOAT", Datatype::VEC_CFLOAT, std::vector<std::complex<float>>)  \
  MACRO("VEC_CDOUBLE", Datatype::VEC_CDOUBLE,                                  \
        std::vector<std::complex<double>>)                                     \
  /* MACRO("VEC_CLONG_DOUBLE", Datatype::VEC_CLONG_DOUBLE,                     \
   * std::vector<std::complex<long double>>) */                                \
  MACRO("VEC_STRING", Datatype::VEC_STRING, std::vector<std::string>)          \
        /*TODO7 MACRO("ARR_DBL_7", Datatype::ARR_DBL_7, array7<double>) */    \
  MACRO("BOOL", Datatype::BOOL, bool)

#define FORALL_SCALAR_OPENPMD_TYPES(MACRO)                                     \
  MACRO("CHAR", Datatype::CHAR, char)                                          \
  MACRO("UCHAR", Datatype::UCHAR, unsigned char)                               \
  MACRO("SHORT", Datatype::SHORT, short)                                       \
  MACRO("INT", Datatype::INT, int)                                             \
  MACRO("LONG", Datatype::LONG, long)                                          \
  MACRO("LONGLONG", Datatype::LONGLONG, long long)                             \
  MACRO("USHORT", Datatype::USHORT, unsigned short)                            \
  MACRO("UINT", Datatype::UINT, unsigned int)                                  \
  MACRO("ULONG", Datatype::ULONG, unsigned long)                               \
  MACRO("ULONGLONG", Datatype::ULONGLONG, unsigned long long)                  \
  MACRO("FLOAT", Datatype::FLOAT, float)                                       \
  MACRO("DOUBLE", Datatype::DOUBLE, double)                                    \
  /* MACRO("LONG_DOUBLE", Datatype::LONG_DOUBLE, long double) */               \
  MACRO("CFLOAT", Datatype::CFLOAT, std::complex<float>)                       \
  MACRO("CDOUBLE", Datatype::CDOUBLE, std::complex<double>)                    \
  /* MACRO("CLONG_DOUBLE", Datatype::CLONG_DOUBLE, std::complex<long           \
   * double>) */                                                               \
  MACRO("STRING", Datatype::STRING, std::string)                               \
  /*TODO7 MACRO("ARR_DBL_7", Datatype::ARR_DBL_7, array7<double>) */    \
  MACRO("BOOL", Datatype::BOOL, bool)

// This C++ version is a bit more tedious to use than the macro version above
template <typename F, typename... Args>
void forall_openPMD_types(const F &f, Args &&...args) {
  f("CHAR", Datatype::CHAR, char{}, std::forward<Args>(args)...);
  f("UCHAR", Datatype::UCHAR, (unsigned char){}, std::forward<Args>(args)...);
  f("SHORT", Datatype::SHORT, short{}, std::forward<Args>(args)...);
  f("INT", Datatype::INT, int{}, std::forward<Args>(args)...);
  f("LONG", Datatype::LONG, long{}, std::forward<Args>(args)...);
  f("LONGLONG", Datatype::LONGLONG, (long long){}, std::forward<Args>(args)...);
  f("USHORT", Datatype::USHORT, (unsigned short){},
    std::forward<Args>(args)...);
  f("UINT", Datatype::UINT, (unsigned int){}, std::forward<Args>(args)...);
  f("ULONG", Datatype::ULONG, (unsigned long){}, std::forward<Args>(args)...);
  f("ULONGLONG", Datatype::ULONGLONG, (unsigned long long){},
    std::forward<Args>(args)...);
  f("FLOAT", Datatype::FLOAT, float{}, std::forward<Args>(args)...);
  f("DOUBLE", Datatype::DOUBLE, double{}, std::forward<Args>(args)...);
  // f("LONG_DOUBLE", Datatype::LONG_DOUBLE, (long double){},
  // std::forward<Args>(args)...);
  f("CFLOAT", Datatype::CFLOAT, std::complex<float>{},
    std::forward<Args>(args)...);
  f("CDOUBLE", Datatype::CDOUBLE, std::complex<double>{},
    std::forward<Args>(args)...);
  // f("CLONG_DOUBLE", Datatype::CLONG_DOUBLE, std::complex<long double>{},
  // std::forward<Args>(args)...);
  f("STRING", Datatype::STRING, std::string{}, std::forward<Args>(args)...);
  f("VEC_CHAR", Datatype::VEC_CHAR, std::vector<char>{},
    std::forward<Args>(args)...);
  f("VEC_UCHAR", Datatype::VEC_UCHAR, std::vector<unsigned char>{},
    std::forward<Args>(args)...);
  f("VEC_SHORT", Datatype::VEC_SHORT, std::vector<short>{},
    std::forward<Args>(args)...);
  f("VEC_INT", Datatype::VEC_INT, std::vector<int>{},
    std::forward<Args>(args)...);
  f("VEC_LONG", Datatype::VEC_LONG, std::vector<long>{},
    std::forward<Args>(args)...);
  f("VEC_LONGLONG", Datatype::VEC_LONGLONG, std::vector<long long>{},
    std::forward<Args>(args)...);
  f("VEC_USHORT", Datatype::VEC_USHORT, std::vector<unsigned short>{},
    std::forward<Args>(args)...);
  f("VEC_UINT", Datatype::VEC_UINT, std::vector<unsigned int>{},
    std::forward<Args>(args)...);
  f("VEC_ULONG", Datatype::VEC_ULONG, std::vector<unsigned long>{},
    std::forward<Args>(args)...);
  f("VEC_ULONGLONG", Datatype::VEC_ULONGLONG, std::vector<unsigned long long>{},
    std::forward<Args>(args)...);
  f("VEC_FLOAT", Datatype::VEC_FLOAT, std::vector<float>{},
    std::forward<Args>(args)...);
  f("VEC_DOUBLE", Datatype::VEC_DOUBLE, std::vector<double>{},
    std::forward<Args>(args)...);
  // f("VEC_LONG_DOUBLE", Datatype::VEC_LONG_DOUBLE, std::vector<long double>{},
  // std::forward<Args>(args)...);
  f("VEC_CFLOAT", Datatype::VEC_CFLOAT, std::vector<std::complex<float>>{},
    std::forward<Args>(args)...);
  f("VEC_CDOUBLE", Datatype::VEC_CDOUBLE, std::vector<std::complex<double>>{},
    std::forward<Args>(args)...);
  // f("VEC_CLONG_DOUBLE", Datatype::VEC_CLONG_DOUBLE,
  // std::vector<std::complex<long double>>{}, std::forward<Args>(args)...);
  f("VEC_STRING", Datatype::VEC_STRING, std::vector<std::string>{},
    std::forward<Args>(args)...);
#warning "TODO7"
  // f("ARR_DBL_7", Datatype::ARR_DBL_7, array7<double>{},
  //   std::forward<Args>(args)...);
  f("BOOL", Datatype::BOOL, bool{}, std::forward<Args>(args)...);
}

namespace {
template <typename T, typename U>
std::vector<std::pair<T, U>> map_to_vector_pair(const std::map<T, U> &m) {
  std::vector<std::pair<std::string, bool>> vp;
  vp.reserve(m.size());
  for (const auto &p : m)
    vp.push_back(p);
  return vp;
}

template <typename T, typename U>
std::vector<std::tuple<T, U>> map_to_vector_tuple(const std::map<T, U> &m) {
  std::vector<std::tuple<std::string, bool>> vp;
  vp.reserve(m.size());
  for (const auto &p : m)
    vp.emplace_back(p.first, p.second);
  return vp;
}

template <typename T> std::shared_ptr<T> create_aliasing_shared_ptr(T *ptr) {
  auto null_deleter = [](T *) {};
  return std::shared_ptr<T>(ptr, null_deleter);
}

template <typename T>
std::shared_ptr<T> capture_vector_as_buffer(std::vector<T> &vec) {
  if constexpr (std::is_same_v<T, bool>) {
    // We cannot handle std::vector<bool> because it is special
    std::abort();
  } else {
    auto deleter = [](T *) { /* do not delete anything */ };
    std::shared_ptr<T> ptr(vec.data(), std::move(deleter));
    return ptr;
  }
}

template <typename T> std::shared_ptr<T> capture_vector(std::vector<T> vec) {
  if constexpr (std::is_same_v<T, bool>) {
    // Copy the vector, because std::vector<bool> is special
    T *dataptr = new T[vec.size()];
    std::shared_ptr<T> ptr(dataptr, std::default_delete<T[]>());
    std::copy(vec.begin(), vec.end(), dataptr);
    return ptr;
  } else {
    // Capture the vector
    T *dataptr = vec.data();
    auto deleter = [vec = std::move(vec)](T *) {
      // We moved the vector into the anonymous function, and thus it will be
      // destructed when the anonymous function is destructed. There is no need
      // to call a destructor manually.
    };
    std::shared_ptr<T> ptr(dataptr, std::move(deleter));
    return ptr;
  }
}

#warning "TODO7"
// template <typename T, std::size_t N>
// void add_array_type(jlcxx::Module &mod, const std::string &name) {
//   mod.add_type<std::array<T, N>>(name)
//       .template constructor<>()
//       .template constructor<const std::array<T, N> &>()
//       .method("size1", &std::array<T, N>::size)
//       .method("getindex1",
//               [](const std::array<T, N> &a, std::size_t n) { return a[n]; });
//   jlcxx::stl::apply_stl<std::array<T, N>>(mod);
// }

// template <typename T, std::size_t N>
// void map_array_type(jlcxx::Module &mod, const std::string &name) {
//   mod.map_type<std::array<T, N>>(name);
//   mod.method("size1", [](const std::array<T, N> &a) { return a.size(); });
//   mod.method("getindex1",
//              [](const std::array<T, N> &a, std::size_t n) { return a[n]; });
//   jlcxx::stl::apply_stl<std::array<T, N>>(mod);
// }

template <typename T, typename U>
void add_pair_type(jlcxx::Module &mod, const std::string &name) {
  mod.add_type<std::pair<T, U>>(name)
      .template constructor<>()
      .template constructor<const std::pair<T, U> &>()
      .method("first", [](const std::pair<T, U> &p) { return p.first; })
      .method("second", [](const std::pair<T, U> &p) { return p.second; });
  jlcxx::stl::apply_stl<std::pair<T, U>>(mod);
}

} // namespace

#warning "TODO7"
// namespace jlcxx {
// template <> struct IsMirroredType<std::array<double, 7>> : std::false_type {};
// } // namespace jlcxx

// We use one function per header file
void define_julia_Access(jlcxx::Module &mod);
void define_julia_Attributable(jlcxx::Module &mod);
void define_julia_Attribute(jlcxx::Module &mod);
void define_julia_BaseRecordComponent(jlcxx::Module &mod);
void define_julia_ChunkInfo(jlcxx::Module &mod);
template <typename Eltype, typename Keytype = std::string>
void define_julia_Container(jlcxx::Module &mod);
void define_julia_Dataset(jlcxx::Module &mod);
void define_julia_Datatype(jlcxx::Module &mod);
void define_julia_Format(jlcxx::Module &mod);
void define_julia_Iteration(jlcxx::Module &mod);
void define_julia_Mesh(jlcxx::Module &mod);
void define_julia_MeshRecordComponent(jlcxx::Module &mod);
void define_julia_RecordComponent(jlcxx::Module &mod);
void define_julia_RecordComponent_load_chunk(
    jlcxx::Module &mod, jlcxx::TypeWrapper<RecordComponent> &type);
void define_julia_RecordComponent_make_constant(
    jlcxx::Module &mod, jlcxx::TypeWrapper<RecordComponent> &type);
void define_julia_RecordComponent_store_chunk(
    jlcxx::Module &mod, jlcxx::TypeWrapper<RecordComponent> &type);
void define_julia_Series(jlcxx::Module &mod);
void define_julia_UnitDimension(jlcxx::Module &mod);
void define_julia_ReadIterations(jlcxx::Module &mod);
void define_julia_WriteIterations(jlcxx::Module &mod);
void define_julia_shared_ptr(jlcxx::Module &mod);
void define_julia_version(jlcxx::Module &mod);

#endif // #ifndef DEFS_HPP
