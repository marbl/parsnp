#ifndef _threadstorage_h_
#define _threadstorage_h_

// aed 9/7/7: MUSCLE v3.6 made prolific use of file and function-scope static variables.
// When running multi-threaded, the shared nature of static variables
// results in nasty race conditions.
//
// The TLS template defines thread-local storage for a particular variable
//

#ifdef _OPENMP
#define MAX_THREAD_COUNT	64
#define OMP_GET_THREAD_NUM	omp_get_thread_num()
#include <omp.h>
#else
#define MAX_THREAD_COUNT	1
#define OMP_GET_THREAD_NUM	0

#endif


#define NELEMS(o)	sizeof(o)/sizeof(o[0])

template<typename T>
class TLS
{
public:
	TLS(){};

	TLS( T t_val ){
		for(int i = 0; i < MAX_THREAD_COUNT; i++)
			t[i] = t_val;
	}

	T& get()
	{
	  //	  printf("Max thread count: %d \n", MAX_THREAD_COUNT);
		return t[OMP_GET_THREAD_NUM];
	}
private:
	TLS( const TLS& tls );	// disallow copying
	T t[MAX_THREAD_COUNT];
};


/**
 * A thread-local storage class specifically to handle string constant initializers.
 * The above TLS code gives the compiler fits for string constants.
 */
template<typename T>
class TLSstr
{
public:
	TLSstr(){};

	TLSstr( T t_val ){
		for(int i = 0; i < MAX_THREAD_COUNT; i++)
			memcpy(t[i], t_val, sizeof(t_val));
	}

	T& get()
	{
		return t[OMP_GET_THREAD_NUM];
	}
private:
	TLSstr( const TLSstr& tls );	// disallow copying
	T t[MAX_THREAD_COUNT];
};


#endif	// _threadstorage_h_

