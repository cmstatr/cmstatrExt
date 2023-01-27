#ifndef RF_ERROR_H
#define RF_ERROR_H

#include <cstdio>
#include <string>
#include <stdexcept>
#include <memory>

template<typename ... Args>
inline void Rf_error( const std::string& format, Args ... args )
{
  int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
  if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
  auto size = static_cast<size_t>( size_s );
  std::unique_ptr<char[]> buf( new char[ size ] );
  std::snprintf( buf.get(), size, format.c_str(), args ... );
  throw std::runtime_error( std::string( buf.get(), buf.get() + size - 1 )); // We don't want the '\0' inside
}

#endif
