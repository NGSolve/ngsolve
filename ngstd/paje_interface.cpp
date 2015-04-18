#include <fstream>
#include <cxxabi.h>
#include <map>
#include <set>
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
        auto DefineContainer = [&] ( int alias, int parent, char const * name ) 
          {
            trace_stream 
              << PajeDefineContainerType << '\t' 
              << alias << '\t' 
              << parent << '\t' 
              << '"' << name << '"' 
              << '\n';
          };

        auto DefineState = [&] ( int alias, int parent, char const * name ) 
          {
            trace_stream 
              << PajeDefineStateType << '\t' 
              << alias << '\t' 
              << parent << '\t' 
              << '"' << name << '"' 
              << '\n';
          };

        auto CreateContainer = [&] ( int alias, int type, int parent, char const * name ) 
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
        auto DefineEntityValue = [&] ( int alias, int type, char const * name, double r, double g, double b) 
          {
            trace_stream 
              << PajeDefineEntityValue << '\t' 
              << alias << '\t' 
              << type << '\t' 
              << '"' << name << '"' << '\t'
              << '"' << r << ' ' << g << ' ' << b << '"' 
              << '\n';
          };
        auto PushState = [&] ( double time, int type, int container, int value )
          {
            trace_stream << PajePushState;
            trace_stream << '\t' << std::setprecision(15) << 1000*time;
            trace_stream << '\t' << type;
            trace_stream << '\t' << container;
//             if(container_id>=0)
//               trace_stream << container_id;
//             if(value_id==-1)
//               trace_stream<< '\t' << '"' << value << '"';
//             else
//               trace_stream << '\t' << value;
//             if(value_id>=0)
              trace_stream << '\t' << value;
              trace_stream << '\n';
          };

        auto PopState = [&] ( double time, int type, int container )
          {
            trace_stream << PajePopState;
            trace_stream << '\t' << std::setprecision(15) << 1000*time;
            trace_stream << '\t' << type;
            trace_stream << '\t' << container;
            trace_stream << '\n';
          };

        trace_stream << header << endl;
        const int container_type_task_manager = 1;
        const int container_type_thread = 2;
        const int container_type_timer = 3;

        DefineContainer( container_type_task_manager, 0, "Task Manager" );
        DefineContainer( container_type_thread, container_type_task_manager, "Thread");
        DefineContainer( container_type_timer, container_type_task_manager, "Timers");

        const int state_type_job = 1;
        const int state_type_task = 2;
        const int state_type_timer = 3;

        DefineState( state_type_job, container_type_task_manager, "Job" );
        DefineState( state_type_task, container_type_thread, "Task" );
        DefineState( state_type_timer, container_type_timer, "Timer state" );

        const int container_task_manager = -1;
        const int container_timer = -2;
        CreateContainer(container_task_manager, container_type_task_manager, 0, "The task manager" );
        CreateContainer(container_timer, container_type_timer, container_task_manager, "Timer" );

        for (int i=0; i<nthreads; i++)
          {
            stringstream name;
            name << "Thread " << i;
            CreateContainer(i, container_type_thread, container_task_manager, name.str().c_str() );
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

        auto Hue2RGB = [] ( double x, double &r, double &g, double &b )
          {
            double d = 1.0/6.0;
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
          };

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
            double r,g,b;
            Hue2RGB(x,r,g,b);
            DefineEntityValue( job.second, state_type_job, job_names[job.first].c_str(), r, g, b );
            DefineEntityValue( job.second, state_type_task, job_names[job.first].c_str(), r, g, b );
          }


        for(Job & j : jobs)
          {
            PushState( j.start_time, state_type_job, container_task_manager, job_map[j.type] );
            PopState( j.stop_time, state_type_job, container_task_manager );
          }

        std::sort (timer_events.begin(), timer_events.end());
        std::set<int> timer_ids;
        for(auto & event : timer_events)
          timer_ids.insert(event.timer_id);
        for(auto id : timer_ids)
          DefineEntityValue( id, state_type_timer, NgProfiler::GetName(id).c_str(), 0.0, 1.0, 1.0 );

        for(auto & event : timer_events)
          {
            if(event.is_start)
              PushState( event.time, state_type_timer, container_timer, event.timer_id );
            else
              PopState( event.time, state_type_timer, container_timer );
          }

        for(auto & vtasks : tasks)
          {
            for (Task & t : vtasks) {
                int value_id = job_map[jobs[t.job_id-1].type];
                PushState( t.start_time, state_type_task, t.thread_id, value_id );
                PopState( t.stop_time, state_type_task, t.thread_id );
            }
          }
        trace_stream.close();
      }
}
