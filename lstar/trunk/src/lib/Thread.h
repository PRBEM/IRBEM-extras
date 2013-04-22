//
//  Thread.h
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/19/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#ifndef UBJDevelopment_Thread_h
#define UBJDevelopment_Thread_h

#include "Definitions.h"
#include <pthread.h>
#include <stdexcept>
#include <cerrno>
#include <cstring> // for strerror
#include <cstdio>
#include <vector>

namespace UBK {

    class Locker;

    //!
    //! Mutex lock class.
    //! Wrapper to pthread's mutex lock.
    //! @sa Locker
    //!
    class Mutex {
        pthread_mutex_t __l;

        Mutex(Mutex const&);
        Mutex& operator=(Mutex const&);

    protected:
        virtual void lock() {
            int err;
            if ( (err=pthread_mutex_lock(&__l)) ) {
                throw std::runtime_error(strerror(err));
            }
        };
        virtual bool tryLock() {
            int err = pthread_mutex_trylock(&__l);
            if (EINVAL == err) {
                throw std::runtime_error(strerror(err));
            }
            return EBUSY != err;
        };
        virtual void unlock() {
            int err;
            if ( (err=pthread_mutex_unlock(&__l)) ) {
                throw std::runtime_error(strerror(err));
            }
        };

    public:
        virtual ~Mutex() {
            int err;
            if ( (err=pthread_mutex_destroy(&__l)) ) {
                fprintf(stderr, "Mutex destroy failed. Reason: %s\n", strerror(err));
                fprintf(stderr, "Terminating.");
                std::terminate();
            }
        };
        Mutex() {
            int err;
            if ( (err=pthread_mutex_init(&__l, NULL)) ) {
                throw std::runtime_error(strerror(err));
            }
        };

        friend class Locker;
    };

    //!
    //! Locker class.
    //! Takes the reference to Mutex object, locks and unlock on deallocation.
    //! @sa Mutex
    //!
    class Locker {
        Mutex& __l;

    public:
        virtual ~Locker() {
            __l.unlock();
        };
        Locker(Mutex& _l) : __l(_l) {
            __l.lock();
        };
    };

    //!
    //! Thread class.
    //! Creates a thread and assigns a task.
    //! @note Functor should implement operator()().
    //! If the passed task throw exceptions, the worker thread calls std::terminate.
    //! @sa ThreadFor
    //!
    template<class _Nullary>
    class Thread {
        BOOL _joinable;
        pthread_t __t;
        _Nullary __f;

        Thread(Thread<_Nullary> const&);
        Thread<_Nullary>& operator=(Thread<_Nullary> const&);

        static void *worker(void *_ctx) {//throw() {
            try {
                Thread<_Nullary>* thiS = reinterpret_cast< Thread<_Nullary>* >(_ctx);
                thiS->__f();
            } catch (std::exception &e) {
                printf("Exception on worker. Reason: %s\n", e.what());
                std::terminate();
            } catch (...) {
                printf("Unknown exception on worker.\n");
                std::terminate();
            }
            pthread_exit(NULL);
        };

    public:
        ~Thread() {
            if (this->joinable()) {
                this->join();
            }
        };
        Thread(_Nullary _f) : __f(_f), _joinable(NO) {
            int err;
            if ( (err=pthread_create(&__t, NULL, Thread::worker, this)) ) {
                throw std::runtime_error(strerror(err));
            }
            _joinable = YES;
        };

        void join() {
            if (!this->joinable()) {
                throw std::runtime_error("Thread is not joinable.");
            }
            int err;
            if ( (err=pthread_join(__t, NULL)) ) {
                throw std::runtime_error(strerror(err));
            }
            _joinable = NO;
        };
        void detach() {
            if (!this->joinable()) {
                throw std::runtime_error("Thread is not joinable.");
            }
            int err;
            if ( (err=pthread_detach(__t)) ) {
                throw std::runtime_error(strerror(err));
            }
            _joinable = NO;
        };
        bool joinable() const {
            return _joinable;
        };
    };

    //!
    //! Loop parallelization.
    //! Indexed-loop (for loop) parallelization.
    //! @note Functor should implement operator()(long);
    //! @sa Thread
    //!
    template<class _Unary>
    class ThreadFor {
        Mutex __key;
        long __wid;
        _Unary __f;

        ThreadFor(ThreadFor<_Unary> const&);
        ThreadFor<_Unary>& operator=(ThreadFor<_Unary> const&);

        long wid() {
            Locker l(__key);
            return --__wid;
        };
        void operator()() {
            long idx;
            while ( (idx=this->wid())>=0 ) {
                __f(idx);
            }
        };

    public:
        typedef ThreadFor<_Unary> self_type;
        typedef Thread<self_type&> thread_type;

        ~ThreadFor() {};
        ThreadFor(unsigned long _n_threads, unsigned long _n_jobs, _Unary _f) :  __wid(_n_jobs),__f(_f) {
            using namespace std;

            if (!_n_threads) {
                throw invalid_argument("_n_threads == 0.");
            }
            if (!_n_jobs) {
                return;
            } else if (_n_threads > _n_jobs) {
                _n_threads = _n_jobs;
            }

            vector< thread_type* > t(_n_threads, NULL);
            for (long idx=0; idx<_n_threads; idx++) {
                t[idx] = new thread_type(*this);
            }
            for (long idx=0; idx<_n_threads; idx++) {
                delete t[idx];
            }
        };

        friend class Thread<self_type&>;
    };

}

#endif
