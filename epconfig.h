#ifndef EPCONFIG_H_
#define EPCONFIG_H_

#include <iostream>
#include <map>
#include <string>

namespace epconfig
{

  //---------------------------------------------------------------------------
  // The configuration::data is a simple map string (key, value) pairs.
  // The file is stored as a simple listing of those pairs, one per line.
  // The key is separated from the value by an equal sign '='.
  // Commentary begins with the first non-space character on the line a hash or
  // semi-colon ('#' or ';').
  //
  // Example:
  //   # This is an example
  //   source.directory = C:\Documents and Settings\Jennifer\My Documents\
  //   file.types = *.jpg;*.gif;*.png;*.pix;*.tif;*.bmp
  //
  // Notice that the configuration file format does not permit values to span
  // more than one line, commentary at the end of a line, or [section]s.
  //   
  struct data : std::map <std::string, std::string>
    {
    // Here is a little convenience method...
    bool iskey( const std::string& s ) const
      {
      return count( s ) != 0;
      }
    };

  //---------------------------------------------------------------------------
  // The extraction operator reads configuration::data until EOF.
  // Invalid data is ignored.
  //
  std::istream& operator >> ( std::istream& ins, data& d );

  //---------------------------------------------------------------------------
  // The insertion operator writes all configuration::data to stream.
  //
  std::ostream& operator << ( std::ostream& outs, const data& d );
}

#endif /* EPCONFIG_H_ */
