#ifndef PAJE_INTERFACE_HPP_INCL__
#define PAJE_INTERFACE_HPP_INCL__
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cxxabi.h>
#include <omp.h>

namespace ngstd
{
  class PajeTrace
    {
      stringstream trace;
      double start_time;
      int nthreads;

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

      vector<vector<Task> > tasks;
      vector<Job> jobs;

    public:
      void Init(int anthreads)
        {
          start_time = omp_get_wtime();

          nthreads = anthreads;

          tasks.resize(nthreads);
          for(auto & t : tasks)
            t.reserve(100000);
          jobs.reserve(100000);
        }

      void StartTask(int thread_id, int task_id, int job_id)
        {
//           Task t;
//           t.start_time = omp_get_wtime() - start_time;
//           t.task_id = task_id;
//           t.job_id = job_id;
//           t.thread_id = thread_id;
//           tasks[thread_id].push_back(t);
          tasks[thread_id].push_back( Task{task_id, job_id, thread_id, omp_get_wtime()-start_time, 0.0} );
        }

      void StopTask(int thread_id)
        {
          tasks[thread_id].back().stop_time = omp_get_wtime() - start_time;
        }

      void StartJob(int job_id, const std::type_info & type)
        {
          jobs.push_back( Job{job_id, &type, omp_get_wtime()-start_time, 0.0 } );
        }

      void StopJob()
        {
          jobs.back().stop_time = omp_get_wtime() - start_time;
        }

      void Write(const char *filename)
        {
          trace << header << endl;
          trace << PajeDefineContainerType << '\t'
            << "TM" << '\t'
            << 0 << '\t'
            << "\"Task Manager\"" << endl;

          trace << PajeDefineContainerType << '\t'
            << "T" << '\t'
            << "TM" << '\t'
            << "\"Thread\"" << endl;

          trace << PajeDefineStateType << '\t'
            << "J" << '\t'
            << "TM" << '\t'
            << "\"Job\"" << endl;

          trace << PajeDefineStateType << '\t'
            << "TS" << '\t'
            << "T" << '\t'
            << "\"Task\"" << endl;

          trace << PajeCreateContainer << '\t'
            << std::setprecision(15) << 1000*(omp_get_wtime()-start_time) << '\t'
            << "tm" << '\t'
            << "TM\t"
            << "0\t"
            << "\"The task manager\"" << endl;

          for (int i=0; i<nthreads; i++)
            {
              trace << PajeCreateContainer << '\t'
                << std::setprecision(15) << 1000*(omp_get_wtime()-start_time) << '\t'
                << 't' << i << '\t'
                << "T\t"
                << "tm\t"
                << "\"Thread " << i << '"' << endl;
            }

          int job_counter = 0;
          map<const std::type_info *, int> job_map;
          map<const std::type_info *, string> job_names;

          for(Job & j : jobs)
            if(job_map.find(j.type) == job_map.end())
              {
                job_map[j.type] = job_counter++;
                int status;
                char * realname = abi::__cxa_demangle(j.type->name(), 0, 0, &status);
                string name = realname;
                name.erase( name.find('}')+1);
                job_names[j.type] = name;
                free(realname);
              }


          int njob_types = job_map.size();
          for( auto & job : job_map )
            {
              string name = job_names[job.first];
              unsigned char u;
              for (char c : name)
                u+=c;
              double x = 1.0*u/255;
              double d = 1.0/6.0;
              double r,g,b;
              if(x<d)
                r=1, g=6*x,b=0;
              else if (x<2*d)
                r=1.0-6*(x-d),g=1,b=0;
              else if (x<3*d)
                r=0, g=1,b=6*(x-2*d);
              else if (x<4*d)
                r=0, g=1-6*(x-3*d),b=1;
              else if (x<5*d)
                r=6*(x-4*d), g=0,b=1;
              else
                r=1, g=0,b=1-5*(x-d);

              trace << PajeDefineEntityValue << '\t'
                << 'J' << job.second << '\t'
                << 'J' << '\t'
                << '"' << job_names[job.first] << '"' << '\t'
                << '"' << r << ' ' << g << ' ' << b << '"' << '\t' << endl;

              trace << PajeDefineEntityValue << '\t'
                << "TS" << job.second << '\t'
                << "TS" << '\t'
                << '"' << job_names[job.first] << '"' << '\t'
                << '"' << r << ' ' << g << ' ' << b << '"' << '\t' << endl;
            }


          for(Job & j : jobs)
            {
              trace << PajePushState << '\t'
                << std::setprecision(15) << 1000*j.start_time << '\t'
                << 'J' << '\t'
                << "tm" << '\t'
                << 'J' << job_map[j.type] << endl;

              trace << PajePopState << '\t'
                << std::setprecision(15) << 1000*j.stop_time << '\t'
                << 'J' << '\t'
                << "tm" << '\t'
                << endl;
            }

          for(auto & vtasks : tasks)
            for (Task & t : vtasks) {
                stringstream value;
                value << "\"taskid: " << t.task_id << ','
                  << jobs[t.job_id-1].job_id << '"';
                trace << PajePushState << '\t'
                  << std::setprecision(15) << 1000*t.start_time << '\t'
                  << "TS" << '\t'
                  << 't' << t.thread_id << '\t'
                  << "TS" << job_map[jobs[t.job_id-1].type] << endl;

                trace << PajePopState << '\t'
                  << std::setprecision(15) << 1000*t.stop_time << '\t'
                  << "TS" << '\t'
                  << 't' << t.thread_id << '\t'
                  << endl;
            }
          ofstream out(filename);
          out << trace.str() << endl;
          out.close();
        }

    private:
      enum
        {
          PajeDefineContainerType = 0,
          PajeDefineVariableType = 1,
          PajeDefineStateType = 2,
          PajeDefineEventType = 3,
          PajeDefineLinkType = 4,
          PajeDefineEntityValue = 5,
          PajeCreateContainer = 6,
          PajeDestroyContainer = 7,
          PajeSetVariable = 8,
          PajeAddVariable = 9,
          PajeSubVariable = 10,
          PajeSetState = 11,
          PajePushState = 12,
          PajePopState = 13,
          PajeResetState = 14,
          PajeStartLink = 15,
          PajeEndLink = 16,
          PajeNewEvent = 17
        };

      static constexpr const char *header =
        "%EventDef PajeDefineContainerType 0 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineVariableType 1 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%       Color color \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineStateType 2 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineEventType 3 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%       Color color \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineLinkType 4 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       StartContainerType string \n"
        "%       EndContainerType string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDefineEntityValue 5 \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Name string \n"
        "%       Color color \n"
        "%EndEventDef \n"
        "%EventDef PajeCreateContainer 6 \n"
        "%       Time date \n"
        "%       Alias string \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeDestroyContainer 7 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Name string \n"
        "%EndEventDef \n"
        "%EventDef PajeSetVariable 8 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value double \n"
        "%EndEventDef\n"
        "%EventDef PajeAddVariable 9 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value double \n"
        "%EndEventDef\n"
        "%EventDef PajeSubVariable 10 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value double \n"
        "%EndEventDef\n"
        "%EventDef PajeSetState 11 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%EndEventDef\n"
        "%EventDef PajePushState 12 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%EndEventDef\n"
        "%EventDef PajePopState 13 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%EndEventDef\n"
        "%EventDef PajeResetState 14 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%EndEventDef\n"
        "%EventDef PajeStartLink 15 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%       StartContainer string \n"
        "%       Key string \n"
        "%EndEventDef\n"
        "%EventDef PajeEndLink 16 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%       EndContainer string \n"
        "%       Key string \n"
        "%EndEventDef\n"
        "%EventDef PajeNewEvent 17 \n"
        "%       Time date \n"
        "%       Type string \n"
        "%       Container string \n"
        "%       Value string \n"
        "%EndEventDef\n";
    };
}

#endif // PAJE_INTERFACE_HPP_INCL__
