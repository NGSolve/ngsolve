#include <fstream>
#include <cxxabi.h>
#include <map>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sstream>
#include "paje_interface.hpp"
#include <ngstd.hpp>

namespace ngstd {
    using std::cout;
    using std::endl;
    using std::string;
    using std::stringstream;

    PajeTrace *trace;
    void PajeTrace::Write(const char *filename)
      {
        cout << "Write traces..." << endl;
        cout << "Number of Jobs: " << jobs.size() << endl;

        std::ofstream trace_stream(filename);
        auto DefineContainer = [&] ( char const * alias, char const * parent, char const * name ) 
          {
            trace_stream 
              << PajeDefineContainerType << '\t' 
              << alias << '\t' 
              << parent << '\t' 
              << '"' << name << '"' 
              << '\n';
          };

        auto DefineState = [&] ( char const * alias, char const * parent, char const * name ) 
          {
            trace_stream 
              << PajeDefineStateType << '\t' 
              << alias << '\t' 
              << parent << '\t' 
              << '"' << name << '"' 
              << '\n';
          };

        auto CreateContainer = [&] ( char const * alias, char const * type, char const * parent, char const * name ) 
          { 
            trace_stream 
              << PajeCreateContainer << '\t' 
              << std::setprecision(15) << 1000*(GetTime()) << '\t'
              << alias << '\t' 
              << type << '\t' 
              << parent << '\t' 
              << '"' << name << '"' 
              << '\n';
          };
        auto DefineEntityValue = [&] ( char const * alias, char const * type, char const * name, double r, double g, double b) 
          {
            trace_stream 
              << PajeDefineEntityValue << '\t' 
              << alias << '\t' 
              << type << '\t' 
              << '"' << name << '"' << '\t'
              << '"' << r << ' ' << g << ' ' << b << '"' 
              << '\n';
          };
        auto PushState = [&] ( double time, char const * type, char const * container, char const * value, int container_id = -1 , int value_id = -1)
          {
            trace_stream << PajePushState;
            trace_stream << '\t' << std::setprecision(15) << 1000*time;
            trace_stream << '\t' << type;
            trace_stream << '\t' << container;
            if(container_id>=0)
              trace_stream << container_id;
            if(value_id==-1)
              trace_stream<< '\t' << '"' << value << '"';
            else
              trace_stream << '\t' << value;
            if(value_id>=0)
              trace_stream << value_id;
              trace_stream << '\n';
          };

        auto PopState = [&] ( double time, char const * type, char const * container , int container_id = -1)
          {
            trace_stream << PajePopState;
            trace_stream << '\t' << std::setprecision(15) << 1000*time;
            trace_stream << '\t' << type;
            trace_stream << '\t' << container;
            if(container_id>=0)
              trace_stream << container_id;
            trace_stream << '\n';
          };

        trace_stream << header << endl;
        DefineContainer( "TM", "0", "Task Manager" );
        DefineContainer( "T", "TM", "Thread");
        DefineContainer( "TIM", "TM", "Timers");

        DefineState( "J", "TM", "Job" );
        DefineState( "TS", "T", "Task" );
        DefineState( "TimS", "TIM", "Timer state" );

        CreateContainer("tm", "TM", "0", "The task manager" );
        CreateContainer("tim", "TIM", "tm", "Timer" );

        for (int i=0; i<nthreads; i++)
          {
            stringstream alias;
            alias << 't' << i;
            stringstream name;
            name << "Thread " << i;

            CreateContainer(alias.str().c_str(), "T", "tm", name.str().c_str() );
          }

        int job_counter = 0;
        std::map<const std::type_info *, int> job_map;
        std::map<const std::type_info *, string> job_names;

        for(Job & j : jobs)
          if(job_map.find(j.type) == job_map.end())
            {
              job_map[j.type] = job_counter++;
              int status;
              char * realname = abi::__cxa_demangle(j.type->name(), 0, 0, &status);
              string name = realname;
              // name.erase( name.find('}')+1);
              job_names[j.type] = name;
              free(realname);
            }


        int njob_types = job_map.size();
        for( auto & job : job_map )
          {
            string name = job_names[job.first];
            unsigned char u = 0;

            for (char c : name)
              {
                u+=c;
                u+=128;
              }

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

            stringstream alias;
            alias << 'J' << job.second;
            DefineEntityValue( alias.str().c_str(), "J", job_names[job.first].c_str(), r, g, b );

            stringstream task_alias;
            task_alias << "TS" << job.second;
            DefineEntityValue( task_alias.str().c_str(), "TS", job_names[job.first].c_str(), r, g, b );

          }

        for(Job & j : jobs)
          {
            stringstream job_value;
            job_value << 'J' << job_map[j.type];
            PushState( j.start_time, "J", "tm", job_value.str().c_str() );
            PopState( j.stop_time, "J", "tm" );
          }

        std::sort (timer_events.begin(), timer_events.end());
        for(auto & event : timer_events)
          {
            if(event.is_start)
              PushState( event.time, "TimS", "tim", NgProfiler::GetName(event.timer_id).c_str() );
            else
              PopState( event.time, "TimS", "tim" );
          }

        for(auto & vtasks : tasks)
          {
            for (Task & t : vtasks) {
                int value_id = job_map[jobs[t.job_id-1].type];
                PushState( t.start_time, "TS", "t", "TS", t.thread_id, value_id );
                PopState( t.stop_time, "TS", "t", t.thread_id );
            }
          }
        trace_stream.close();
      }
}
