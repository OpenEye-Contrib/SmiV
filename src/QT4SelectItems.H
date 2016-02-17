//
// file QT4SelectItems.H
// Dave Cosgrove
// AstraZeneca
// 27th June 2006
//
// Puts up a QListWidget in a box allowing selections from a list of strings.
// Similar to QTSelectItems.H but in Qt4 speak.
// This must always be used as a modal dialog, since the address of the vector
// holding the selections is kept, and this must be in scope when it is filled.
// It is filled when either the ok or cancel buttons are pressed.

#ifndef DAC_QT4_SELECT_ITEMS__
#define DAC_QT4_SELECT_ITEMS__

#include <string>
#include <vector>

#include <QDialog>

using namespace std;

// ***************************************************************************

class QListWidget;
class QListWidgetItem;
class QVBoxLayout;

class QT4SelectItems : public QDialog {

  Q_OBJECT

protected slots :

  virtual void okClicked();
  virtual void cancelClicked();
  void slot_list_double_clicked( QListWidgetItem *item ); // select item and out

protected :

  vector<char>  *selected_items_; // pointer to what's passed in

  QWidget       *action_box_;
  QVBoxLayout   *vlayout_;
  QListWidget   *list_widget_;

  void build_action_box();

public :

  // selected_items needs to be big enough to hold the result, and can hold
  // pre-selections if required. I.e. any item in item_labels that has a
  // corresponding 1 in selected_items will be set as selected when the dialog
  // pops up.
  QT4SelectItems( const string &label , vector<QString> &item_labels ,
                  vector<char> &selected_items , bool radio_box = false ,
                  QWidget *par = 0 );
  ~QT4SelectItems() {}

};

#endif
