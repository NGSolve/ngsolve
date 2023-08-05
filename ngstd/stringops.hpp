#ifndef _STRINGOPS_HPP
#define _STRINGOPS_HPP

namespace ngstd
{

  // NGS_DLL_HEADER bool StringFitsPattern(const string & str, const string & pattern);
  NGS_DLL_HEADER bool StringFitsPattern(string_view str, string_view pattern);

}

#endif // _STRINGOPS_HPP
