#if NDIM == 1
#define KJI_DO(lo,hi) i=lo,hi
#define KJI_DO_VEC(hi) do i=1,hi(1)
#define CLOSE_DO
#define IJK i
#define DTIMES(TXT) TXT
#elif NDIM == 2
#define KJI_DO(lo,hi) j=lo,hi; do i=lo,hi
#define KJI_DO_VEC(hi) j=1,hi(2); do i=1,hi(1)
#define CLOSE_DO end do
#define IJK i, j
#define DTIMES(TXT) TXT, TXT
#elif NDIM == 3
#define KJI_DO(lo,hi) j=lo,hi; do i=lo,hi; do k=lo,hi
#define KJI_DO_VEC(hi) k=1,hi(3); do j=1,hi(2); do i=1,hi(1)
#define CLOSE_DO end do; end do
#define IJK i, j, k
#define DTIMES(TXT) TXT, TXT, TXT
#endif
