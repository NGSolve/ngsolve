#ifndef PAJE_INTERFACE_HPP_INCL__
#define PAJE_INTERFACE_HPP_INCL__

#include <limits>
#include <omp.h>
#include "array.hpp"

namespace ngstd
{
  extern class PajeTrace *trace;
  class PajeTrace
    {
    public:
      double start_time;
      int nthreads;

      static size_t max_tracefile_size;


      // Approximate number of events to trace. Tracing will
      // be stopped if any thread reaches this number of events
      int max_num_events_per_thread;
      bool tracing_enabled;

      static void SetMaxTracefileSize( size_t max_size ) {
          max_tracefile_size = max_size;
      }

      std::string tracefile_name;

      struct Job
        {
          int job_id;
          const std::type_info *type;
          double start_time;
          double stop_time;
        };

      struct Task
        {
          int thread_id;

          int id;
          int id_type;

          int additional_value;

          double start_time;
          double stop_time;

          static constexpr int ID_NONE = -1;
          static constexpr int ID_JOB = 1;
          static constexpr int ID_TIMER = 2;
        };

      struct TimerEvent
        {
          int timer_id;
          double time;
          bool is_start;
          int thread_id;

          bool operator < (const TimerEvent & other) const { return time < other.time; }
        };

      struct ThreadLink
        {
          int thread_id;
          int key;
          double time;
          bool is_start;
          bool operator < (const ThreadLink & other) const { return time < other.time; }
        };

      Array<Array<Task> > tasks;
      Array<Job> jobs;
      Array<TimerEvent> timer_events;
      Array<Array<ThreadLink> > links;

      double GetTime()
        {
          return omp_get_wtime() - start_time;
        }


    public:
      void StopTracing()
        {
          if(tracing_enabled)
            {
              cout << "Maximum number of traces reached, tracing is stopped now. To increase the tracefile size, set in the pde file:" << endl;
              cout << "flags tracer = -max_size=size_in_megabytes" << endl;
            }
          tracing_enabled = false;
        }

      PajeTrace(int anthreads, std::string aname = "" )
        {
          start_time = omp_get_wtime();


          nthreads = anthreads;
          tracefile_name = aname;

          int bytes_per_event=33;
          max_num_events_per_thread = min2( (size_t)std::numeric_limits<int>::max, max_tracefile_size/bytes_per_event/(2*nthreads+1)*10/7);
          cout << "Tracefile size = " << max_tracefile_size/1024/1024 << "MB." << endl;
          cout << "Tracing " << max_num_events_per_thread << " events per thread." << endl;

          tasks.SetSize(nthreads);
          int reserve_size = min2(1000000, max_num_events_per_thread);
          for(auto & t : tasks)
            t.SetAllocSize(reserve_size);

          links.SetSize(nthreads);
          for(auto & l : links)
            l.SetAllocSize(reserve_size);

          jobs.SetAllocSize(reserve_size);
          timer_events.SetAllocSize(reserve_size);

          tracing_enabled = true;
        }

      ~PajeTrace()
        {
          if(tracefile_name.size()>0)
            Write(tracefile_name);
        }

      void StartTimer(int timer_id)
        {
          if(timer_events.Size() == max_num_events_per_thread)
            StopTracing();
          if(tracing_enabled)
            timer_events.Append(TimerEvent{timer_id, GetTime(), true});
        }

      void StopTimer(int timer_id)
        {
          if(timer_events.Size() == max_num_events_per_thread)
            StopTracing();
          if(tracing_enabled)
            timer_events.Append(TimerEvent{timer_id, GetTime(), false});
        }

      int StartTask(int thread_id, int id, int id_type = Task::ID_NONE, int additional_value = -1)
        {
          if(tasks[thread_id].Size() == max_num_events_per_thread)
            StopTracing();
          if(tracing_enabled)
            {
              int task_num = tasks[thread_id].Size();
              tasks[thread_id].Append( Task{thread_id, id, id_type, additional_value, GetTime(), 0.0} );
              return task_num;
            }
          return -1;
        }

      void StopTask(int thread_id, int task_num)
        {
          if(task_num>=0)
            tasks[thread_id][task_num].stop_time = GetTime();
        }

      void StartJob(int job_id, const std::type_info & type)
        {
          if(jobs.Size() == max_num_events_per_thread)
            StopTracing();
          if(tracing_enabled)
            jobs.Append( Job{job_id, &type, GetTime(), 0.0 } );
        }

      void StopJob()
        {
          if(tracing_enabled)
            jobs.Last().stop_time = GetTime();
        }

      void StartLink(int thread_id, int key)
        {
          if(links[thread_id].Size() == max_num_events_per_thread)
            StopTracing();
          if(tracing_enabled)
            links[thread_id].Append( ThreadLink{thread_id, key, GetTime(), true} );
        }

      void StopLink(int thread_id, int key)
        {
          if(links[thread_id].Size() == max_num_events_per_thread)
            StopTracing();
          if(tracing_enabled)
            links[thread_id].Append( ThreadLink{thread_id, key, GetTime(), false} );
        }

      void Write( std::string filename );

    };
}

#endif // PAJE_INTERFACE_HPP_INCL__
