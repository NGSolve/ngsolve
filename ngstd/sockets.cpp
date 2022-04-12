// based on the socket class by Rob Tougher
// http://www.linuxgazette.com/issue74/tougher.html


#ifdef WIN32
#include <winsock2.h>
#endif

#include "sockets.hpp"

/*
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
*/
#include <errno.h>
#include <fcntl.h>

namespace ngstd
{

#ifdef WIN32
  static bool sockets_warmed_up = false;
#endif



  Socket :: Socket() 
    : m_sock ( -1 )
  {
    memset ( &m_addr, 0, sizeof ( m_addr ) );
    
#ifdef WIN32
    if(!sockets_warmed_up)
      {
	WSADATA wsa;
	if(WSAStartup(MAKEWORD(1,1),&wsa))
	  {
	    cerr << "WSAStartup() failed, " << GetLastError() << endl;
	  }
      }
#endif
  }

  Socket::~Socket()
  {
    if (is_valid())
#ifdef WIN32		
      ::closesocket (m_sock);
#else
    ::close (m_sock);
#endif
  }



  void Socket::create()
  {
    m_sock = socket (AF_INET, SOCK_STREAM, 0);
    
    if (m_sock == -1)
      throw SocketException ("call to socket failed");

    // TIME_WAIT - argh
    int on = 1;
    if (setsockopt (m_sock, SOL_SOCKET, SO_REUSEADDR, 
                    ( const char* ) &on, sizeof ( on ) ) == -1 )
      {
        throw SocketException (GetLatestError());
      }
  }



  void Socket::bind (int port)
  {
    if ( !is_valid() ) throw SocketException ("not a valid socket");

    m_addr.sin_family = AF_INET;
    m_addr.sin_addr.s_addr = INADDR_ANY;
    m_addr.sin_port = htons ( port );

    int bind_return = ::bind ( m_sock,
			       ( struct sockaddr * ) &m_addr,
			       sizeof ( m_addr ) );

    if (bind_return == -1)
      throw SocketException (GetLatestError());
  }


  void Socket::listen() const
  {
    if ( !is_valid() ) throw SocketException ("not a valid socket");

    int listen_return = ::listen ( m_sock, MAXCONNECTIONS );
    if ( listen_return == -1 )
      throw Exception (GetLatestError());
  }


  bool Socket::accept (Socket & new_socket) const
  {
    socklen_t addr_length = sizeof ( m_addr );
    new_socket.m_sock = ::accept (m_sock, ( sockaddr * )&m_addr, &addr_length );

#ifdef WIN32
    if ( new_socket.m_sock == INVALID_SOCKET )
      return false;
    else
      return true;
#else
    if ( new_socket.m_sock <= 0 )
      {
        if (errno == EINTR) // interrupted system call (e.g. exit)
          return true;

        throw SocketException (GetLatestError());
      }
    else
      return true;
#endif
  }


  void Socket::send ( const std::string & s ) const
  {
    int length = s.length();
    int status = ::send (m_sock, &length, sizeof(length), MSG_NOSIGNAL );    
    if (status <= 0)
      throw SocketException (GetLatestError());

    status = ::send (m_sock, s.c_str(), length+1, 0 /* MSG_NOSIGNAL */ );    
    if (status != length+1)
      {
        cout << "length = " << length << ", status = " << status << endl;
      }
    if (status <= 0)
      throw SocketException (GetLatestError());      
  }


  void Socket::recv (string & s) const
  {
    // new version by Joachim
    int length;
    int status = ::recv (m_sock, &length, sizeof(length), 0 );

    if (status == 0)
      {
        s = "";
        return;
      }

    if (status < 0)
      {
        s = "";
        throw SocketException (GetLatestError());
      }
    
    char * hstr = new char[length+1];

    status = ::recv (m_sock, hstr, length+1, MSG_WAITALL );    
    if (status != length+1)
      {
        s = "";
        delete [] hstr;
        cout << "receive, status = " << status << endl;
        throw SocketException (GetLatestError());
      }

    s = hstr;
    delete [] hstr;
  }



  void Socket::connect ( const string & host, int port )
  {
    if (!is_valid()) throw SocketException ("not a valid socket");

    m_addr.sin_family = AF_INET;
    m_addr.sin_port = htons ( port );

#ifdef WIN32
    struct sockaddr_storage addr;
    int len = sizeof(addr);

    char * dummych; 
    dummych = new char[host.size()+1]; 
    strcpy(dummych,host.c_str());

    int status = WSAStringToAddress(dummych, AF_INET, NULL,(LPSOCKADDR)&addr, &len);
    delete [] dummych;
    latesterror = WSAGetLastError();
    memcpy(&(m_addr.sin_addr),&((struct sockaddr_in *) &addr)->sin_addr,4);

    if ( status != 0)
      throw SocketException (GetLatestError());      

#else
    int status = inet_pton ( AF_INET, host.c_str(), &m_addr.sin_addr );
    if (errno == EAFNOSUPPORT )
      throw Exception ("EAFNOSUPPORT");
#endif
    status = ::connect ( m_sock, ( sockaddr * ) &m_addr, sizeof ( m_addr ) );
    if (status != 0)
      throw SocketException (GetLatestError());
  }


  void Socket::set_non_blocking ( const bool b )
  {
#ifdef WIN32
    cerr << "Socket::set_non_blocking not yet implemented for Windows" << endl;
    exit(10);
#else

    int opts = fcntl (m_sock, F_GETFL);

    if (opts < 0) return;

    if ( b )
      opts = (opts | O_NONBLOCK);
    else
      opts = (opts & ~O_NONBLOCK);

    fcntl (m_sock, F_SETFL,opts);
#endif
  }


  string Socket::GetLatestError(void) const
  {
#ifdef WIN32
    latesterror = WSAGetLastError();

    if(latesterror == WSAEADDRINUSE)
      return "Address already in use";
    if(latesterror == WSAECONNABORTED)
      return "Software caused connection abort";
    if(latesterror == WSAECONNREFUSED)
      return "Connection refused";
    if(latesterror == WSAECONNRESET)
      return "Connection reset by peer";
    if(latesterror == WSAEDESTADDRREQ)
      return "Destination address required";
    if(latesterror == WSAEHOSTUNREACH)
      return "No route to host";
    if(latesterror == WSAEMFILE)
      return "Too many open files";
    if(latesterror == WSAENETDOWN)
      return "Network is down";
    if(latesterror == WSAENETRESET)
      return "Network dropped connection";
    if(latesterror == WSAENOBUFS)
      return "No buffer space available";
    if(latesterror == WSAENETUNREACH)
      return "Network is unreachable";
    if(latesterror == WSAETIMEDOUT)
      return "Connection timed out";
    if(latesterror == WSAHOST_NOT_FOUND)
      return "Host not found";
    if(latesterror == WSASYSNOTREADY)
      return "Network sub-system is unavailable";
    if(latesterror == WSANOTINITIALISED)
      return "WSAStartup() not performed";
    if(latesterror == WSANO_DATA)
      return "Valid name, no data of that type";
    if(latesterror == WSANO_RECOVERY)
      return "Non-recoverable query error";
    if(latesterror == WSATRY_AGAIN)
      return "Non-authoritative host found";
    if(latesterror == WSAVERNOTSUPPORTED)
      return "Wrong WinSock DLL version";
#else
    latesterror = errno;

    if(latesterror == EACCES)
      return "The requested address is protected, and the current user has inadequate permission to access it.";
    if(latesterror == EADDRINUSE)
      return "The specified address is already in use.";
    if(latesterror == EADDRNOTAVAIL)
      return "The specified address is invalid or not available from the local machine, or for AF_CCITT sockets which use wild card addressing, the specified address space overlays the address space of an existing bind.";
    if(latesterror == EAFNOSUPPORT)
      return "The specified address is not a valid address for the address family of this socket.";
    if(latesterror == EBADF)
      return "no valid file descriptor";
    if(latesterror == EDESTADDRREQ)
      return "No addr parameter was specified.";
    if(latesterror == EFAULT)
      return "addr is not a valid pointer.";
    if(latesterror == EINVAL)
      return "The socket is already bound to an address, the socket has been shut down, addrlen is a bad value, or an attempt was made to bind() an AF_UNIX socket to an NFS-mounted (remote) name.";
    if(latesterror == ENETDOWN)
      return "The x25ifname field name specifies an interface that was shut down, or never initialized, or whose Level 2 protocol indicates that the link is not working: Wires might be broken, the interface hoods on the modem are broken, the modem failed, the phone connection failed (this error can be returned by AF_CCITT only), noise interfered with the line for a long period of time.";
    if(latesterror == ENETUNREACH)
      return "The X.25 Level 2 protocol is down. The X.25 link is not working: Wires might be broken, or connections are loose on the interface hoods at the modem, the modem failed, or noise interfered with the line for an extremely long period of time.";
    if(latesterror == ENOBUFS)
      return "No buffer space is available. The bind() cannot complete.";
    if(latesterror == ENOMEM)
      return "No memory is available. The bind() cannot complete.";
    if(latesterror == ENODEV)
      return "The x25ifname field name specifies a nonexistent interface. (This error can be returned by AF_CCITT only.)";
    if(latesterror == ENOTSOCK)
      return "s is a valid file descriptor, but it is not a socket.";
    if(latesterror == EOPNOTSUPP)
      return "The socket referenced by s does not support address binding.";
    if(latesterror == EISCONN)
      return "The connection is already bound. (AF_VME_LINK.)";
    if(latesterror == EOPNOTSUPP)
      return "The socket is not of a type that supports the operation.";
    if(latesterror == ECONNREFUSED)
      return "No one listening on the remote address.";
    if(latesterror == ETIMEDOUT)
      return "Timeout while attempting connection. The server may be too busy to accept new connections. Note that for IP sockets the timeout may be very long when syncookies are enabled on the server.";
    if(latesterror == ENETUNREACH)
      return "Network is unreachable.";
    if(latesterror == EINPROGRESS)
      return "The socket is non-blocking and the connection cannot be completed immediately. It is possible to select(2) or poll(2) for completion by selecting the socket for writing. After select indicates writability, use getsockopt(2) to read the SO_ERROR option at level SOL_SOCKET to determine whether connect completed successfully (SO_ERROR is zero) or unsuccessfully (SO_ERROR is one of the usual error codes listed here, explaining the reason for the failure).";
    if(latesterror == EALREADY)
      return "The socket is non-blocking and a previous connection attempt has not yet been completed.";
    if(latesterror == EAGAIN)
      return "No more free local ports or insufficient entries in the routing cache. For PF_INET see the net.ipv4.ip_local_port_range sysctl in ip(7) on how to increase the number of local ports.";
    if(latesterror == EPERM)
      return "The user tried to connect to a broadcast address without having the socket broadcast flag enabled or the connection request failed because of a local firewall rule.";
#endif

    
    // from errno-base.h:
    if (errno == EINTR)
      return "Interrupted system call";

    return "Unknown error.";
  }


  ServerSocket::ServerSocket (int port)
  {
    try
      {
        Socket::create();
        Socket::bind (port);
        Socket::listen();
      }
    catch (SocketException& e)
      {
        e.Append ("\nduring server socket creation");
        throw e;
      }
  }

  const ServerSocket& ServerSocket::operator << ( const std::string& s ) const
  {
    Socket::send (s);
    return *this;
  }

  const ServerSocket& ServerSocket::operator >> (string & s) const
  {
    Socket::recv (s);
    return *this;
  }


  void ServerSocket::accept ( ServerSocket& sock )
  {
    if ( ! Socket::accept ( sock ) )
      {
	throw SocketException ( "Could not accept socket." );
      }
  }





  ClientSocket::ClientSocket (int port, const string & host)
  {
    try
      {
        Socket::create();
        Socket::connect (host, port);
      }
    catch (SocketException& e)
      {
        e.Append ("\nduring client socket initialization");
        throw e;
      }
  }


  const ClientSocket& ClientSocket::operator << ( const std::string& s ) const
  {
    Socket::send ( s );
    return *this;
  }

  const ClientSocket& ClientSocket::operator >> ( std::string& s ) const
  {
    Socket::recv (s);
    return *this;
  }
}



