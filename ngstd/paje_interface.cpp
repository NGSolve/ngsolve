#include <cxxabi.h>
#include <map>
#include <set>
#include <algorithm>

#include <ngstd.hpp>
#include "paje_interface.hpp"

namespace ngstd
{
  size_t PajeTrace::max_tracefile_size;

  using std::string;
  class PajeFile
    {
    public:
      static void Hue2RGB ( double x, double &r, double &g, double &b )
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

      int alias_counter;

//       std::ofstream trace_stream;
      FILE * ctrace_stream;

//       void Write( double x )
//         {
//           trace_stream << std::setprecision(15) << '\t' << x;
//         }
// 
//       void WriteTime( double x )
//         {
//           trace_stream << std::setprecision(15) << '\t' << 1000*x;
//         }
// 
//       void WriteAlias( int n )
//         {
//           trace_stream << '\t' << 'a' << n;
//         }
// 
//       void Write( int n )
//         {
//           trace_stream << '\t' << n;
//         }
// 
//       void Write( string s )
//         {
//           trace_stream << '\t' << '"' << s << '"';
//         }
// 
//       void WriteColor( double r, double g, double b )
//         {
//           trace_stream << '\t' << '"' << r << ' ' << g << ' ' << b << '"';
//         }
// 
//       void WriteColor( double hue )
//         {
//           double r,g,b;
//           Hue2RGB( hue, r, g, b );
//           trace_stream << '\t' << '"' << r << ' ' << g << ' ' << b << '"';
//         }

      struct PajeEvent
        {
          double time;
          std::function<void()> func;
          bool operator < (const PajeEvent & other) const { return time < other.time; }
        };

      std::vector<PajeEvent> events;

    public:
      enum PType
        {
          CONTAINER=0,
          VARIABLE=1,
          STATE=2,
          EVENT=3,
          LINK=4
        };

      PajeFile( string filename )
//         : trace_stream(filename, ios::binary)
        {
//           trace_stream << header;
          ctrace_stream = fopen ("ngs0.trace","w");
          fprintf(ctrace_stream, "%s", header );
          alias_counter = 0;
        }
      int DefineContainerType ( int parent_type, string name )
        {
          int alias = ++alias_counter;
//           trace_stream << PajeDefineContainerType;
//           WriteAlias(alias);
//           if(parent_type!=0) WriteAlias(parent_type);
//           else Write(parent_type);
//           Write(name);
//           trace_stream << '\n';
          if(parent_type!=0)
            fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\n", PajeDefineContainerType, alias, parent_type, name.c_str() );
          else
            fprintf( ctrace_stream, "%d\ta%d\t%d\t\"%s\"\n", PajeDefineContainerType, alias, parent_type, name.c_str() );
          return alias;
        }

      int DefineVariableType ( int container_type, string name )
        {
          int alias = ++alias_counter;
//           trace_stream << PajeDefineVariableType;
//           WriteAlias(alias);
//           WriteAlias(container_type);
//           Write(name);
//           WriteColor( 1.0, 1.0, 1.0 );
//           trace_stream << '\n';
          fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\t\"1.0 1.0 1.0\"\n", PajeDefineVariableType, alias, container_type, name.c_str() );
          return alias;
        }

      int DefineStateType ( int type, string name )
        {
          int alias = ++alias_counter;
//           trace_stream << PajeDefineStateType;
//           WriteAlias(alias);
//           WriteAlias((int)type);
//           Write(name);
//           trace_stream << '\n';
          fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\n", PajeDefineStateType, alias, type, name.c_str() );
          return alias;
        }

//       int DefineEventType ()
//         {
//           Write("event not implemented");
//         }

      int DefineLinkType (int parent_container_type, int start_container_type, int stop_container_type, string name)
        {
          int alias = ++alias_counter;
//           trace_stream << PajeDefineLinkType;
//           WriteAlias(alias);
//           WriteAlias(parent_container_type);
//           WriteAlias(start_container_type);
//           WriteAlias(stop_container_type);
//           Write(name);
//           trace_stream << '\n';
          fprintf( ctrace_stream, "%d\ta%d\ta%d\ta%d\ta%d\t\"%s\"\n", PajeDefineLinkType, alias, parent_container_type, start_container_type, stop_container_type, name.c_str() );
          return alias;
        }

      int DefineEntityValue (int type, string name, double hue = -1)
        {
          if(hue==-1)
            {
              std::hash<string> shash;
              size_t h = shash(name);
              h ^= h>>32;
              h = (uint32_t)h;
              hue = h*1.0/std::numeric_limits<uint32_t>::max();
            }

          int alias = ++alias_counter;
//           trace_stream << PajeDefineEntityValue;
//           WriteAlias(alias);
//           WriteAlias(type);
//           Write(name);
//           WriteColor(hue);
//           trace_stream << '\n';
          double r,g,b;
          Hue2RGB( hue, r, g, b );
          fprintf( ctrace_stream, "%d\ta%d\ta%d\t\"%s\"\t\"%.15g %.15g %.15g\"\n", PajeDefineEntityValue, alias, type, name.c_str(), r,g,b );
          return alias;
        }

      int CreateContainer ( int type, int parent, string name )
        {
          int alias = ++alias_counter;
//           trace_stream << PajeCreateContainer;
//           Write(0.0);
//           WriteAlias(alias);
//           WriteAlias(type);
//           if(parent!=0) WriteAlias(parent);
//           else Write(parent);
//           Write(name);
//           trace_stream << '\n';
          if(parent!=0)
            fprintf( ctrace_stream, "%d\t0\ta%d\ta%d\ta%d\t\"%s\"\n", PajeCreateContainer, alias, type, parent, name.c_str() );
          else
            fprintf( ctrace_stream, "%d\t0\ta%d\ta%d\t%d\t\"%s\"\n", PajeCreateContainer, alias, type, parent, name.c_str() );
          return alias;
        }
      void DestroyContainer ()
        {}

      void SetVariable (double time, int type, int container, double value )
        {
          events.push_back( PajeEvent {time, [=]() {
//           trace_stream << PajeSetVariable;
//           WriteTime(time);
//           WriteAlias(type);
//           WriteAlias(container);
//           Write(value);
//           trace_stream << '\n';
          fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\t%.15g\n", PajeSetVariable, time, type, container, value );
                           } } );
        }

      void AddVariable (double time, int type, int container, double value )
        {
          events.push_back( PajeEvent {time, [=]() {
//           trace_stream << PajeAddVariable;
//           WriteTime(time);
//           WriteAlias(type);
//           WriteAlias(container);
//           Write(value);
//           trace_stream << '\n';
          fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\t%.15g\n", PajeAddVariable, time, type, container, value );
                           } } );
        }

      void SubVariable (double time, int type, int container, double value )
        {
          events.push_back( PajeEvent {time, [=]() {
//           trace_stream << PajeSubVariable;
//           WriteTime(time);
//           WriteAlias(type);
//           WriteAlias(container);
//           Write(value);
//           trace_stream << '\n';
          fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\t%.15g\n", PajeSubVariable, time, type, container, value );
                           } } );
        }

      void SetState ()
        {}

      void PushState ( double time, int type, int container, int value, int id = 0, bool value_is_alias = true )
        {
          events.push_back( PajeEvent {time, [=]() {
//                            trace_stream << PajePushState;
//                            WriteTime(time);
//                            WriteAlias(type);
//                            WriteAlias(container);
//                            if(value_is_alias) WriteAlias(value);
//                            else Write(value);
//                            Write(id);
//                            trace_stream << '\n';
                           if(value_is_alias)
                             fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\ta%d\t%d\n", PajePushState, time, type, container, value, id);
                           else
                             fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\t%d\t%d\n", PajePushState, time, type, container, value, id);
                           } } );
        }

      void PopState ( double time, int type, int container )
        {
          events.push_back( PajeEvent {time, [=]() {
//                            trace_stream << PajePopState;
//                            WriteTime(time);
//                            WriteAlias(type);
//                            WriteAlias(container);
//                            trace_stream << '\n';
                           fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\n", PajePopState, time, type, container );
                           } } );
        }

      void ResetState ()
        {}

      void StartLink ( double time, int type, int container, int value, int start_container, int key )
        {
          events.push_back( PajeEvent {time, [=]() {
//           trace_stream << PajeStartLink;
//           WriteTime(time);
//           WriteAlias(type);
//           WriteAlias(container);
//           Write(value);
//           WriteAlias(start_container);
//           Write(key);
//           trace_stream << '\n';
              fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\t%d\ta%d\t%d\n", PajeStartLink, time, type, container, value, start_container, key );
          } } );
        }

      void EndLink ( double time, int type, int container, int value, int end_container, int key )
        {
          events.push_back( PajeEvent {time, [=]() {
//           trace_stream << PajeEndLink;
//           WriteTime(time);
//           WriteAlias(type);
//           WriteAlias(container);
//           Write(value);
//           WriteAlias(end_container);
//           Write(key);
//           trace_stream << '\n';
              fprintf( ctrace_stream, "%d\t%.15g\ta%d\ta%d\t%d\ta%d\t%d\n", PajeEndLink, time, type, container, value, end_container, key );
          } } );
        }

      void NewEvent ()
        {}

      void WriteEvents()
        {
          cout << "Wite trace file..." << endl;
          double t = omp_get_wtime();
          cout << __LINE__ << endl;
          std::sort (events.begin(), events.end());
          cout << __LINE__ << endl;
          for( auto & event : events )
            event.func();
          cout << __LINE__ << endl;
          cout << "took " << omp_get_wtime() - t << " seconds." << endl;
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
        "%       Id string \n"
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

  PajeTrace *trace;

  void PajeTrace::Write( string filename )
    {
      std::cout << "Write traces..." << std::endl;
      std::cout << "Number of Jobs: " << jobs.Size() << std::endl;
      if(!tracing_enabled)
        {
          cout << "Tracing stopped during computation due to tracefile size limit of " << max_tracefile_size/1024/1024 << " megabytes." << endl;
          cout << "To increase the limit, set in the pde file:" << endl;
          cout << "flags tracer = -max_size=size_in_megabytes" << endl << endl;
        }

      PajeFile paje(filename);

      const int container_type_task_manager = paje.DefineContainerType( 0, "Task Manager" );
      const int container_type_node = paje.DefineContainerType( container_type_task_manager, "Node");
      const int container_type_thread = paje.DefineContainerType( container_type_task_manager, "Thread");
      const int container_type_timer = container_type_thread; //paje.DefineContainerType( container_type_task_manager, "Timers");
      const int container_type_jobs = paje.DefineContainerType( container_type_task_manager, "Jobs");

      const int state_type_job = paje.DefineStateType( container_type_jobs, "Job" );
      const int state_type_task = paje.DefineStateType( container_type_thread, "Task" );
      const int state_type_timer = paje.DefineStateType( container_type_timer, "Timer state" );

      const int variable_type_active_threads = paje.DefineVariableType( container_type_jobs, "Active threads" );

      const int container_task_manager = paje.CreateContainer( container_type_task_manager, 0, "The task manager" );
      const int container_jobs = paje.CreateContainer( container_type_jobs, container_task_manager, "Jobs" );
      paje.SetVariable( 0.0, variable_type_active_threads, container_jobs, 0.0 );
      const int container_node0 = paje.CreateContainer( container_type_node, container_task_manager, "Node 0" );

      std::vector <int> thread_aliases;
      for (int i=0; i<nthreads; i++)
        {
          char name[20];
          sprintf(name, "Thread %d", i);
          thread_aliases.push_back( paje.CreateContainer( container_type_thread, container_node0, name ) );
        }

      int job_counter = 0;
      std::map<const std::type_info *, int> job_map;
      std::map<const std::type_info *, int> job_task_map;

      for(Job & j : jobs)
        if(job_map.find(j.type) == job_map.end())
          {
            int status;
            char * realname = abi::__cxa_demangle(j.type->name(), 0, 0, &status);
            string name = realname;
            job_map[j.type] = paje.DefineEntityValue( state_type_job, name, -1 );
            job_task_map[j.type] = paje.DefineEntityValue( state_type_task, name, -1 );
            // name.erase( name.find('}')+1);
            free(realname);
          }

      for(Job & j : jobs)
        {
          paje.PushState( j.start_time, state_type_job, container_jobs, job_map[j.type] );
          paje.PopState( j.stop_time, state_type_job, container_jobs );
        }

//       QuickSort(timer_events);
      std::set<int> timer_ids;
      std::map<int,int> timer_aliases;

      for(auto & event : timer_events)
        timer_ids.insert(event.timer_id);

      for(auto & vtasks : tasks)
        for (Task & t : vtasks)
          if(t.id_type==Task::ID_TIMER)
            timer_ids.insert(t.id);

      for(auto id : timer_ids)
        timer_aliases[id] = paje.DefineEntityValue( state_type_timer, NgProfiler::GetName(id).c_str(), -1 );

      int timerdepth = 0;
      int maxdepth = 0;
      for(auto & event : timer_events)
        {
          if(event.is_start)
            {
              timerdepth++;
              maxdepth = timerdepth>maxdepth ? timerdepth : maxdepth;
            }
          else
              timerdepth--;
        }

      std::vector<int> timer_container_aliases;
      timer_container_aliases.resize(maxdepth);
      for(int i=0; i<maxdepth; i++)
        {
          char name[30];
          sprintf(name, "Timer level %d", i);
          timer_container_aliases[i] =  paje.CreateContainer( container_type_timer, container_task_manager, name );
        }

      timerdepth = 0;
      for(auto & event : timer_events)
        {
          if(event.is_start)
            paje.PushState( event.time, state_type_timer, timer_container_aliases[timerdepth++], timer_aliases[event.timer_id] );
          else
            paje.PopState( event.time, state_type_timer, timer_container_aliases[--timerdepth] );
        }

      for(auto & vtasks : tasks)
        {
          for (Task & t : vtasks) {
              int value_id = t.id;

              switch(t.id_type)
                {
                case Task::ID_JOB:
                  value_id = job_task_map[jobs[t.id-1].type];
                  paje.AddVariable( t.start_time, variable_type_active_threads, container_jobs, 1.0 );
                  paje.SubVariable( t.stop_time, variable_type_active_threads, container_jobs, 1.0 );
                  paje.PushState( t.start_time, state_type_task, thread_aliases[t.thread_id], value_id, t.additional_value, true );
                  paje.PopState( t.stop_time, state_type_task, thread_aliases[t.thread_id] );
                  break;
                case Task::ID_TIMER:
                  value_id = timer_aliases[t.id];
                  paje.PushState( t.start_time, state_type_timer, thread_aliases[t.thread_id], value_id, t.additional_value, true );
                  paje.PopState( t.stop_time, state_type_timer, thread_aliases[t.thread_id] );
                  break;
                default:
                  paje.PushState( t.start_time, state_type_task, thread_aliases[t.thread_id], value_id, t.additional_value, false );
                  paje.PopState( t.stop_time, state_type_task, thread_aliases[t.thread_id] );
                  break;
                }


          }

        }
      paje.WriteEvents();

      // Merge link event
      int nlinks = 0;
      for( auto & l : links)
        nlinks += l.Size();

      std::vector<ThreadLink> links_merged;
      links_merged.reserve(nlinks);
      std::vector<int> pos(nthreads);

      int nlinks_merged = 0;
      while(nlinks_merged < nlinks)
        {
          int minpos = 0;
          double mintime = 1e300;
          for (int t = 0; t<nthreads; t++)
            {
              if(pos[t] < links[t].Size() && links[t][pos[t]].time < mintime)
                {
                  minpos = t;
                  mintime = links[t][pos[t]].time;
                }
            }
          links_merged.push_back( links[minpos][pos[minpos]] );
          pos[minpos]++;
          nlinks_merged++;
        }

      std::vector<ThreadLink> started_links;

      int link_type = paje.DefineLinkType(container_type_node, container_type_thread, container_type_thread, "links");

      // match links
      for ( auto & l : links_merged )
        {
          if(l.is_start)
            {
              started_links.push_back(l);
            }
          else
            {
              int i = 0;
              while(i<started_links.size())
                {
                  while(i<started_links.size() && started_links[i].key == l.key)
                    {
                      ThreadLink & sl = started_links[i];
                      // Avoid links on same thread
                      if(sl.thread_id != l.thread_id)
                        {
                          paje.StartLink( sl.time, link_type, container_node0, l.key, thread_aliases[sl.thread_id], l.key);
                          paje.EndLink(    l.time, link_type, container_node0, l.key, thread_aliases[l.thread_id], l.key);
                        }
                      started_links.erase(started_links.begin()+i);
                    }
                  i++;
                }
            }
        }
    }
}
