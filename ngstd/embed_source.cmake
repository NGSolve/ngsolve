# usage: cmake -DINPUT=<file> -DOUTPUT=<cpp> -DSYMBOL=<name> -P embed_source.cmake
# embeds INPUT verbatim as "extern const std::string SYMBOL"

if(NOT INPUT OR NOT OUTPUT OR NOT SYMBOL)
  message(FATAL_ERROR "usage: cmake -DINPUT=<file> -DOUTPUT=<cpp> -DSYMBOL=<name> -P embed_source.cmake")
endif()

file(READ "${INPUT}" content)
get_filename_component(inputname "${INPUT}" NAME)

# MSVC truncates a string literal above 16380 bytes (C2026), and the limit
# applies to adjacent literals after concatenation - so join chunks at run time
set(chunksize 15000)
string(LENGTH "${content}" len)

# SYMBOL may be qualified (ns::name): the definition is placed in that namespace
string(FIND "${SYMBOL}" "::" nspos REVERSE)
if(nspos GREATER -1)
  string(SUBSTRING "${SYMBOL}" 0 ${nspos} ns)
  math(EXPR namepos "${nspos} + 2")
  string(SUBSTRING "${SYMBOL}" ${namepos} -1 name)
  set(nsopen "namespace ${ns} {\n")
  set(nsclose "}\n")
else()
  set(name "${SYMBOL}")
  set(nsopen "")
  set(nsclose "")
endif()

set(body "// generated from ${inputname} - do not edit\n#include <string>\n${nsopen}extern const std::string ${name};\nconst std::string ${name} = std::string()\n")
set(pos 0)
while(pos LESS len)
  string(SUBSTRING "${content}" ${pos} ${chunksize} chunk)
  string(APPEND body "  + R\"ngs_embed(${chunk})ngs_embed\"\n")
  math(EXPR pos "${pos} + ${chunksize}")
endwhile()
string(APPEND body "  ;\n${nsclose}")

file(WRITE "${OUTPUT}" "${body}")
