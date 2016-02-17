//
// file SmiVRecord.cc
// David Cosgrove
// AstraZeneca
// 4th January 2010
//

#include "SmiVRecord.H"

using namespace std;
using namespace OEChem;

namespace DACLIB {
  void apply_daylight_aromatic_model( OEMolBase &mol );
}

// ****************************************************************************
SmiVRecord::SmiVRecord() {

}

// ****************************************************************************
SmiVRecord::SmiVRecord( const string &smi , const string &smi_name ) :
    in_smi_( smi ) , smi_name_( smi_name ) {

}

// ****************************************************************************
SmiVRecord::SmiVRecord( const OEMolBase &mol ) {

  OECreateSmiString( in_smi_ , mol , OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo );
  OECreateIsoSmiString( can_smi_ , mol );
  smi_name_ = mol.GetTitle();

}

// ****************************************************************************
SmiVRecord::~SmiVRecord() {

}

// ****************************************************************************
void SmiVRecord::create_can_smi() {

  OEMol mol;
  OEParseSmiles( mol , in_smi_ );
  DACLIB::apply_daylight_aromatic_model( mol );
  OECreateIsoSmiString( can_smi_ , mol );

}
