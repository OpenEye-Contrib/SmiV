//
// file QTSmilesEditDialog.H
// Dave Cosgrove
// AstraZeneca
// 20th January 2010
//
// This dialog allows the input of a SMILES string with interactive
// drawing of corresponding molecule.

#ifndef DAC_QT_SMILES_EDIT_DIALOG
#define DAC_QT_SMILES_EDIT_DIALOG

#include <QDialog>

class QLabel;
class QLineEdit;
class QPushButton;
class QString;
class QVBoxLayout;

namespace DACLIB {
  class QTMolDisplay2D;
}

// *****************************************************************************

namespace DACLIB {

  class QTSmilesEditDialog : public QDialog {

  Q_OBJECT

  public :

    QTSmilesEditDialog( QWidget *parent = 0 , Qt::WindowFlags f = 0 );

    void set_smiles( const QString &new_smi );
    QString get_smiles() const;

    void set_name( const QString &new_name );
    QString get_name() const;

  protected :

    QVBoxLayout *vbox_;
    QLineEdit *smiles_box_ , *name_box_;
    QPushButton *apply_button_;

    DACLIB::QTMolDisplay2D *mol_disp_;

    bool data_changed_;

  private :

    virtual void build_widget();
    virtual QWidget *build_action_box();

  private slots :

    void slot_ok_clicked();
    void slot_apply_clicked();
    void slot_update_mol_display();
    void slot_atom_selected( unsigned int oe_ind );
    void slot_smiles_cursor_moved( int old_pos , int new_pos );
    void slot_data_changed();

  signals :

    void atom_selected( unsigned int oe_ind );
    void new_smiles( QString new_smi , QString new_name );
    void new_smiles_other();

  };

} // EO namespace DACLIB

#endif
