//
// file QTSmartsEditDialog.H
// Dave Cosgrove
// 31st October 2007
//
// This class puts up a modal dialog that allows the user to change a SMARTS
// definition and see a picture of it.

#ifndef DAC_QT_SMARTS_EDIT_DIALOG
#define DAC_QT_SMARTS_EDIT_DIALOG

#include <string>
#include <vector>

#include <QDialog>
#include <QLineEdit>

using namespace std;

// ************************************************************************

class QString;
class QLabel;
class QHBoxLayout;
class QVBoxLayout;

namespace DACLIB {
  class QTMolDisplay2D;
}

class QTSmartsEditDialog : public QDialog {

Q_OBJECT

public :

  QTSmartsEditDialog( bool with_name_input ,
                      QWidget *p = 0 , Qt::WindowFlags f = 0 );

  void set_smarts( const QString &new_smarts );
  QString get_smarts() const {
    return smarts_box_->text();
  }
  void set_smarts_name( const QString &new_name );
  QString get_smarts_name() const;
  
  void set_smarts_sub_defns( const std::vector<std::pair<std::string,std::string> > &new_subs ) {
    sub_defs_ = new_subs;
  }

protected :

  QVBoxLayout *vbox_; // vertical stack of elements in the widget
  QHBoxLayout *mol_disp_hbox_; /* holds smarts_disp_ and anything else that
                                  sub-classes might need */
  QLineEdit *smarts_box_ , *smarts_name_box_;
  DACLIB::QTMolDisplay2D *smarts_disp_;

  virtual QString get_full_smarts();
  virtual bool ok_to_close();

private :

  QLabel *status_label_;

  bool with_name_input_; // whether input of a SMARTS name is present
  std::vector<std::pair<std::string,std::string> > sub_defs_; // for parsing vector bindings

  virtual void build_widget();
  QWidget *build_action_box();

private slots :

  void slot_display_smarts();
  void slot_ok_clicked();

signals :

  void smarts_atom_selected( unsigned int oe_ind );

};

#endif
