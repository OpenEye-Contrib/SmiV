//
// file QTSmartsIntPickDialog.cc
// Dave Cosgrove
// AstraZeneca
// 12th January 2010
//

#include "DACOEMolAtomIndex.H"
#include "QTMolDisplay2D.H"
#include "QTSmartsIntPickDialog.H"

#include <QHBoxLayout>

#include <oechem.h>

#include <boost/scoped_ptr.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {
  // in apply_daylight_arom_model_to_oemol.cc
  void apply_daylight_aromatic_model( OEMolBase &mol );
}

// ****************************************************************************
QTSmartsIntPickDialog::QTSmartsIntPickDialog( const string &smi ,
                                              bool with_name_input ,
                                              QWidget *parent ,
                                              Qt::WindowFlags f ) :
QTSmartsEditDialog( with_name_input , parent , f ) {

  add_int_pick_display( smi );
  connect( smarts_box_ , SIGNAL( returnPressed() ) ,
           this , SLOT( slot_select_mol_disp_atoms() ) );

}

// ****************************************************************************
void QTSmartsIntPickDialog::slot_int_pick_atom_selected() {

  vector<OEAtomBase *> sel_atoms = int_pick_disp_->sel_atoms();

  QString smarts = build_smarts_from_disp_mol_sel_atoms( sel_atoms );
  set_smarts( smarts );

}

// ****************************************************************************
void QTSmartsIntPickDialog::slot_select_mol_disp_atoms() {

  QString q_smarts = get_full_smarts();
  if( q_smarts == "Bad SMARTS") {
    return;
  }
  string smarts = q_smarts.toLocal8Bit().data();

  int_pick_disp_->select_atoms_by_smarts( smarts );

}

// ****************************************************************************
void QTSmartsIntPickDialog::add_int_pick_display( const string &smi ) {

  int_pick_disp_ = new DACLIB::QTMolDisplay2D;
  mol_disp_hbox_->insertWidget( 0 , int_pick_disp_ );

  disp_mol_.reset( OENewMolBase( OEMolBaseType::OEDefault ) );
  OEParseSmiles( *disp_mol_ , smi );
  DACLIB::apply_daylight_aromatic_model( *disp_mol_ );

  int_pick_disp_->set_display_molecule( disp_mol_.get() );

  connect( int_pick_disp_ , SIGNAL( atom_selected( unsigned int ) ) ,
           this , SLOT( slot_int_pick_atom_selected() ) );

}

// ****************************************************************************
QString QTSmartsIntPickDialog::build_smarts_from_disp_mol_sel_atoms( vector<OEAtomBase *> &sel_atoms ) {

  scoped_ptr<OEMolBase> tmp_mol( OENewMolBase( *disp_mol_ , OEMolBaseType::OEDefault ) );
  vector<char> in_mol( DACLIB::max_atom_index( *tmp_mol ) , 0 );
  for( int i = 0 , is = sel_atoms.size() ; i < is ; ++i ) {
    in_mol[DACLIB::atom_index( *sel_atoms[i] )] = 1;
  }

  for( OEIter<OEAtomBase> atom = tmp_mol->GetAtoms() ; atom ; ++atom ) {
    if( !in_mol[DACLIB::atom_index( *atom )] ) {
      tmp_mol->DeleteAtom( atom );
    }
  }

  string smi;
  OECreateSmiString( smi , *tmp_mol );

  return QString( smi.c_str() );

}
