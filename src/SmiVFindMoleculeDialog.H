//
// file SmivFindMoleculeDialog.H
// David Cosgrove
// AstraZeneca
// 15th January 2010
//
// Dialog for finding a molecule in SmiV by name. Bit more general than that,
// probably, but that's what it was written for.
//

#ifndef DAC_SMIV_FIND_MOLECULE_DIALOG
#define DAC_SMIV_FIND_MOLECULE_DIALOG

#include <QDialog>

class QComboBox;
class QLineEdit;

// ****************************************************************************

class SmiVFindMoleculeDialog : public QDialog {

  Q_OBJECT

public :

  SmiVFindMoleculeDialog( QWidget *parent = 0 , Qt::WindowFlags f = 0 );

private :

  QComboBox *search_mode_;
  QLineEdit *mol_name_;

  void build_widget();
  QWidget *build_action_box();

  void emit_search_details();

private slots :

  void slot_ok_clicked();
  void slot_apply_clicked();

signals :

  void molecule_search_name( QString mol_name , int search_mode );

};

#endif // DAC_SMIV_FIND_MOLECULE_DIALOG
