//
// File SmiVPanel.H
// David Cosgrove
// AstraZeneca
// 4th January 2010
//
// This class holds the display and slider for a set of SmivRecords.

#ifndef DAC_SMIV_PANEL
#define DAC_SMIV_PANEL

#include <string>
#include <vector>

#include <QWidget>

#include <boost/shared_ptr.hpp>

// ********************************************************************************

class SmiVRecord;
class QCheckBox;
class QKeyEvent;
class QLabel;
class QLineEdit;
class QMouseEvent;
class QSlider;
class QString;

namespace DACLIB {
  class QTMolDisplay2D;
}

namespace OEChem {
  class OESubSearch;
}

typedef boost::shared_ptr<SmiVRecord> pSmiVRec;

// ********************************************************************************

class SmiVPanel : public QWidget {

  Q_OBJECT

public :

  SmiVPanel( QWidget *parent = 0 , Qt::WindowFlags f = 0 );

  void add_data( const std::vector<pSmiVRec> &new_recs );
  void set_title( const QString &new_title );
  void set_subsearches( const std::vector<std::pair<boost::shared_ptr<OEChem::OESubSearch>,std::string> > &ss );

  // search from current position down for the named molecule, using the given
  // search mode - 0 for Exact Match, 1 for Starts With, 2 for Contains,
  // and display molecule if found.
  // returns whether sucessful or not.
  bool show_molecule( std::string mol_name , int search_mode );

  boost::shared_ptr<SmiVRecord> current_smiv_rec() const;
  std::vector<pSmiVRec> smiv_recs() const { return smiv_recs_; }

  bool is_selected() const;
  void set_selected( bool new_val );

  void write_smiles_to_stream( std::ostream &os ) const;

  void go_to_first_mol();
  void go_to_last_mol();

protected :

  DACLIB::QTMolDisplay2D *mol_disp_;
  QSlider *mol_slider_;
  QLineEdit *in_smi_;
  QLineEdit *can_smi_;
  QLabel *msg_ , *title_; // arbitrary messages
  QCheckBox *sel_box_;

  std::vector<pSmiVRec> smiv_recs_;
  std::vector<std::pair<pSmiVRec,int> > dropped_recs_; // the record and its original sequence number
  std::vector<std::pair<boost::shared_ptr<OEChem::OESubSearch>,std::string> > sub_searches_; // used for colouring molecules
  bool selected_;

  void keyPressEvent( QKeyEvent *event );

  void build_widget();
  void colour_atoms(); // using sub_searches_
  void drop_current_mol();
  void undo_last_drop();
  // move the mol_slider_ by the given step, if possible
  void change_current_mol( int step );

private slots :

  void slot_mol_slider_changed();
  void slot_selection_box_changed();

signals :

  void selection_box_changed( QWidget *wid );
  void new_display_mol( QString mol_name );

};

#endif // DAC_SMIV_PANEL
