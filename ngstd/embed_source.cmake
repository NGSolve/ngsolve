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

set(body "// generated from ${inputname} - do not edit\n#include <string>\nextern const std::string ${SYMBOL};\nconst std::string ${SYMBOL} = std::string()\n")
set(pos 0)
while(pos LESS len)
  string(SUBSTRING "${content}" ${pos} ${chunksize} chunk)
  string(APPEND body "  + R\"ngs_embed(${chunk})ngs_embed\"\n")
  math(EXPR pos "${pos} + ${chunksize}")
endwhile()
string(APPEND body "  ;\n")

file(WRITE "${OUTPUT}" "${body}")
