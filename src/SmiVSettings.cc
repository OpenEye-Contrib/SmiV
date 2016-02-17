//
// file SmiVSettings.cc
// David Cosgrove
// 10th July 2009

#include "SmiVSettings.H"

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <iostream>

using namespace std;
namespace po = boost::program_options;

// *****************************************************************************
SmiVSettings::SmiVSettings( int argc , char **argv ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 0 );
  }

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();

}

// ***************************************************************************
bool SmiVSettings::operator!() const {

  return false;

}

// ***************************************************************************
void SmiVSettings::print_usage( ostream &os ) const {

  os << usage_text_ << endl;

}

// **************************************************************************
void SmiVSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text" )
    ( "molecule-file,M" , po::value<string>( &mol_file_ ) ,
      "Input molecule filename" )
    ( "mdl-query-file,Q" , po::value<string>( &mdl_file_ ) ,
      "Input MDL substructure query file")
    ( "smarts-file,S" , po::value<string>( &smarts_file_ ) ,
      "Input SMARTS filename" )
    ( "data-file,D" , po::value<string>( &data_file_ ) ,
      "Arbitrary data filename" );

}

