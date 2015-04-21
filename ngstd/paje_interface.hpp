#ifndef PAJE_INTERFACE_HPP_INCL__
#define PAJE_INTERFACE_HPP_INCL__

#include <vector>
#include <omp.h>

namespace ngstd
{
  extern class PajeTrace *trace;
  class PajeTrace
    {
    public:
      double start_time;
      int nthreads;
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
          int task_id;
          int job_id;
          int thread_id;
          double start_time;
          double stop_time;
        };

      struct TimerEvent
        {
          int timer_id;
          double time;
          bool is_start;

          bool operator < (const TimerEvent & other) { return time < other.time; }
        };

      struct ThreadLink
        {
          int thread_id;
          int key;
          double time;
          bool is_start;
          bool operator < (const ThreadLink & other) { return time < other.time; }
        };

      std::vector<std::vector<Task> > tasks;
      std::vector<Job> jobs;
      std::vector<TimerEvent> timer_events;
      std::vector<std::vector<ThreadLink> > links;

      double GetTime()
        {
          return omp_get_wtime() - start_time;
        }


    public:
      PajeTrace(int anthreads, std::string aname = "" )
        {
          start_time = omp_get_wtime();

          nthreads = anthreads;
          tracefile_name = aname;

          tasks.resize(nthreads);
          for(auto & t : tasks)
            t.reserve(1000000);

          links.resize(nthreads);
          for(auto & l : links)
            l.reserve(1000000);

          jobs.reserve(1000000);
          timer_events.reserve(1000000);
        }

      ~PajeTrace()
        {
          if(tracefile_name.size()>0)
            Write(tracefile_name);
        }

      void StartTimer(int timer_id)
        {
          timer_events.push_back(TimerEvent{timer_id, GetTime(), true});
        }

      void StopTimer(int timer_id)
        {
          timer_events.push_back(TimerEvent{timer_id, GetTime(), false});
        }

      void StartTask(int thread_id, int task_id, int job_id = -1)
        {
          tasks[thread_id].push_back( Task{task_id, job_id, thread_id, GetTime(), 0.0} );
        }

      void StopTask(int thread_id)
        {
          tasks[thread_id].back().stop_time = GetTime();
        }

      void StartJob(int job_id, const std::type_info & type)
        {
          jobs.push_back( Job{job_id, &type, GetTime(), 0.0 } );
        }

      void StopJob()
        {
          jobs.back().stop_time = GetTime();
        }

      void StartLink(int thread_id, int key)
        {
          links[thread_id].push_back( ThreadLink{thread_id, key, GetTime(), true} );
        }

      void StopLink(int thread_id, int key)
        {
          links[thread_id].push_back( ThreadLink{thread_id, key, GetTime(), false} );
        }

      void Write( std::string filename );

    };
}

#endif // PAJE_INTERFACE_HPP_INCL__
