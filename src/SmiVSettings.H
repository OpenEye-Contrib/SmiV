//
// file SmiVSettings.H
// David Cosgrove
// 10th July 2009
//
// This is a public class that holds all the settings that the user used to
// start smiv with, to make them easier to pass round the place.

#ifndef DAC_SMIV_SETTINGS
#define DAC_SMIV_SETTINGS

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

// **************************************************************************

class SmiVSettings {

public :

  SmiVSettings( int argc , char **argv );

  bool operator!() const;

  void print_usage( std::ostream &os ) const;

  const std::string &mol_file() { return mol_file_; }
  const std::string &smarts_file() { return smarts_file_; }
  const std::string &mdl_file() { return mdl_file_; }
  const std::string &data_file() { return data_file_; }
  const std::string &usage_text() { return usage_text_; }

private :

  std::string mol_file_;
  std::string smarts_file_;
  std::string mdl_file_; // substructure query in MDL mol file format
  std::string data_file_;
  std::string usage_text_;

  void build_program_options( boost::program_options::options_description &desc );

};

#endif
