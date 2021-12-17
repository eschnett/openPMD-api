#include "defs.hpp"

#include "Container.hpp"

#include <cstdint>

#warning "TODO"
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  std::cerr<<"define_julia_module.0\n";
#warning "TODO7"
  std::cerr<<"add_array_type\n";
  add_array_type<double, 7>(mod, "array_double_7");
  std::cerr<<"add_pair_type\n";
  add_pair_type<std::string, bool>(mod, "pair_string_bool");

  std::cerr<<"define_julia_shared_ptr\n";
  define_julia_shared_ptr(mod);

  // The order of these calls matters. Julia types need to be defined before
  // they are used.

  // Stand-alone classes
  std::cerr<<"define_julia_Access\n";
  define_julia_Access(mod);
  define_julia_ChunkInfo(mod);
  define_julia_Datatype(mod);
  define_julia_Format(mod);
  define_julia_UnitDimension(mod);
  // All classes below need at least Datatype

  std::cerr<<"define_julia_Attribute\n";
  define_julia_Attribute(mod);
  define_julia_Attributable(mod);
  define_julia_Dataset(mod);

  std::cerr<<"define_julia_BaseRecordComponent\n";
  define_julia_BaseRecordComponent(mod); // needs: Attributable
  define_julia_RecordComponent(mod);     // needs: BaseRecordComponent
  define_julia_MeshRecordComponent(mod); // needs: RecordComponent

  std::cerr<<"define_julia_Container<MeshRecordComponent>\n";
  define_julia_Container<MeshRecordComponent>(
      mod); // needs: Attributable, MeshRecordComponent

  std::cerr<<"define_julia_Mesh\n";
  define_julia_Mesh(mod); // needs: Container<MeshRecordComponent>

  std::cerr<<"define_julia_Container<Mesh>\n";
  define_julia_Container<Mesh>(mod); // needs: Attributable

  std::cerr<<"define_julia_Iteration\n";
  define_julia_Iteration(mod); // needs: Attributable, Container<Mesh>

  std::cerr<<"define_julia_Container<Iteration,uint64_t>\n";
  define_julia_Container<Iteration, std::uint64_t>(mod); // needs: Attributable

  std::cerr<<"define_julia_WriteIterations\n";
  define_julia_WriteIterations(mod); // needs: Iteration

  // The main class
  std::cerr<<"define_julia_Series\n";
  define_julia_Series(mod);

  // Handle metadata
  std::cerr<<"define_julia_version\n";
  define_julia_version(mod);
  std::cerr<<"define_julia_module.9\n";
}
