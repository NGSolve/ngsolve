#ifndef _NGSOLVE_SOCKETS_HPP
#define _NGSOLVE_SOCKETS_HPP

// based on the socket class by Rob Tougher
// http://www.linuxgazette.com/issue74/tougher.html



#include <ngstd.hpp>


#include <iostream>
#include <typeinfo>
#include <sys/types.h>



#ifdef WIN32

#include <winsock.h>

#define socklen_t int
#define MSG_NOSIGNAL 0

#else

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <arpa/inet.h>

#define SOCKET int

#endif

#include <string>

#define SOCKETCOMMANDLENGTH 2


namespace ngstd
{






  class SocketException : public Exception 
  {
  public:
    SocketException (const string & s) : Exception (s) { ; }
    virtual ~SocketException () { ; }
  };





  
  // const int MAXHOSTNAME = 200;
  const int MAXCONNECTIONS = 5;
  // const int MAXRECV = 500;
  // const int MAXRECV = 10000;
  
  class Socket
  {
  private:
    
    SOCKET m_sock;
    sockaddr_in m_addr;
  
    mutable int latesterror;
    
  public:
    Socket();
    virtual ~Socket();

    // Server initialization
    void create();
    void bind (int port);
    void listen() const;
    bool accept ( Socket& ) const;

    // Client initialization
    void connect (const string & host, int port);

    // Data Transimission
    void send (const string & s) const;
    void send (string & s) const 
    {
      send ( (const string &)s);
    }

    void recv (std::string & str) const;

    template <typename T>
    void Tsend (const T & data) const
    {
      int status = ::send (m_sock, &data, sizeof(data), MSG_NOSIGNAL );
      if (status < 0)
        throw SocketException (string("problem sending ")
                               + typeid(T).name() + " " 
                               + ToString(data) + string("\n"));
    }

    template <typename T>
    void Trecv (T & data)
    {
      int status = ::recv (m_sock, &data, sizeof(data), MSG_WAITALL);
      if (status < 0)
        throw SocketException (string("problem receiving ")
                               + typeid(T).name() + string("\n"));
    }


    void set_non_blocking ( const bool );

    bool is_valid() const { return m_sock != -1; }

    virtual string GetLatestError(void) const;

  };





  class ServerSocket : public Socket
  {
  public:
    ServerSocket (int port);
    ServerSocket () { ; };
    virtual ~ServerSocket() { ; }

    const ServerSocket& operator << (const string & str) const;
    const ServerSocket& operator >> (string & str) const;
    
    void accept (ServerSocket & sock);
  };



  class ClientSocket : public Socket
  {

  public:
    ClientSocket (int port, const string & host="localhost" );
    virtual ~ClientSocket(){};

    const ClientSocket& operator << ( const std::string& ) const;
    const ClientSocket& operator >> ( std::string& ) const;

  };




  

}

#endif // _NGSOLVE_SOCKETS_HPP

