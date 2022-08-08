#pragma once

#ifdef _STATIC_CPPLIB
   #define _DECLSPEC
#elif defined LIBEDMDLL_EXPORTS
   #define    __declspec(dllexport)
#else
   #define    __declspec(dllimport)
#endif