#if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#define EXPORTFUNC __declspec(dllexport)
#else
#define EXPORTFUNC __attribute__((visibility("default")))
#endif
