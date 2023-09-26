#ifdef DGP_INLINE
#undef DGP_INLINE
#endif

#ifndef DGP_STATIC_LIBRARY
#  define DGP_INLINE inline
#else
#  define DGP_INLINE
#endif