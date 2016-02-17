//
// File SmiVRecord.H
// David Cosgrove
// AstraZeneca
// 4th January 2010
//
// This class holds various bits of interest for a SmivPanel, such as SMILES strings
// for the molecule.

#ifndef DAC_SMIV_RECORD
#define DAC_SMIV_RECORD

#include <string>

#include <oechem.h>

// ****************************************************************************

class SmiVRecord {

public :

  SmiVRecord();
  SmiVRecord( const std::string &smi , const std::string &smi_name );
  SmiVRecord( const OEChem::OEMolBase &mol );
  ~SmiVRecord();

  std::string in_smi() { return in_smi_; }
  std::string can_smi() { return can_smi_; }
  std::string smi_name() { return smi_name_; }
  void set_smi_name( const std::string &new_name ) {
    smi_name_ = new_name;
  }

  void create_can_smi(); // from in_smi_, via an OEMol

protected :

  std::string in_smi_;
  std::string smi_name_;
  std::string can_smi_;

};

#endif // DAC_SMIV_RECORD
