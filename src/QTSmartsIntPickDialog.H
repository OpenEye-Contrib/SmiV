//
// file QTSmartsIntPickDialog.H
// Dave Cosgrove
// AstraZeneca
// 12th January 2010
//
// This class derives from QTSmartsEditDialog, and takes adds a QTMolDisplay2D
// from which the user can pick atoms to define the SMARTS.

#ifndef DAC_QT_SMARTS_INT_PICK_DIALOG
#define DAC_QT_SMARTS_INT_PICK_DIALOG

#include "QTSmartsEditDialog.H"

#include <string>

#include <boost/scoped_ptr.hpp>

#include <oechem.h>

class QString;

// ****************************************************************************

class QTSmartsIntPickDialog : public QTSmartsEditDialog {

  Q_OBJECT

public :

  QTSmartsIntPickDialog( const std::string &smi , bool with_name_input = true ,
                         QWidget *parent = 0 , Qt::WindowFlags f = 0 );

public slots :

  void slot_int_pick_atom_selected();
  void slot_select_mol_disp_atoms();

private :

  DACLIB::QTMolDisplay2D *int_pick_disp_;
  boost::scoped_ptr<OEChem::OEMolBase> disp_mol_;

  void add_int_pick_display( const std::string &smi );

  QString build_smarts_from_disp_mol_sel_atoms( std::vector<OEChem::OEAtomBase *> &sel_atoms );

};

#endif // DAC_QT_SMARTS_INT_PICK_DIALOG
