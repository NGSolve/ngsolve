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

      int StartTask(int thread_id, int id, int id_type = Task::ID_NONE, int additional_value = -1)
        {
          int task_num = tasks[thread_id].size();
          tasks[thread_id].push_back( Task{thread_id, id, id_type, additional_value, GetTime(), 0.0} );
          return task_num;
        }

      void StopTask(int thread_id, int task_num)
        {
          tasks[thread_id][task_num].stop_time = GetTime();
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
