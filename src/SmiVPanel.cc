//
// File SmiVPanel.cc
// David Cosgrove
// AstraZeneca
// 5th January 2010
//

#include "SmiVPanel.H"
#include "SmiVRecord.H"

#include "QTMolDisplay2D.H"

#include <QCheckBox>
#include <QKeyEvent>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMouseEvent>
#include <QSlider>

#include <oechem.h>

#include <boost/scoped_ptr.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;

// ****************************************************************************
SmiVPanel::SmiVPanel( QWidget *parent , Qt::WindowFlags f ) :
QWidget( parent , f ) , selected_( false ) {

  build_widget();

}

// ****************************************************************************
void SmiVPanel::add_data( const vector<pSmiVRec> &new_recs ) {

  if( new_recs.empty() ) {
    smiv_recs_.clear();
    mol_slider_->setDisabled( true );
    mol_disp_->clear_display_molecule();
    in_smi_->setText( "" );
    can_smi_->setText( "" );
    msg_->setText( "" );
  } else {
    smiv_recs_ = new_recs;
    mol_slider_->setEnabled( true );
    mol_slider_->setMinimum( 0 );
    mol_slider_->setMaximum( int( smiv_recs_.size() - 1 ) );
    slot_mol_slider_changed();
  }

  sub_searches_.clear();
  dropped_recs_.clear();

}

// ****************************************************************************
void SmiVPanel::set_subsearches( const vector<pair<boost::shared_ptr<OESubSearch>,string> > &ss ) {

  sub_searches_ = ss;
  colour_atoms();

}

// ****************************************************************************
// search from current position + 1 down for the named molecule, using the given
// search mode - 0 for Exact Match, 1 for Starts With, 2 for Contains
bool SmiVPanel::show_molecule( string mol_name , int search_mode ) {

  int start_num = 0 == search_mode ? 0 : mol_slider_->value() + 1;
  bool found_mol( false );
  for( int i = start_num , is = smiv_recs_.size() ; i < is ; ++i ) {
    switch( search_mode ) {
    case 0 :
      if( smiv_recs_[i]->smi_name() == mol_name ) {
        found_mol = true;
      }
      break;
    case 1 :
      if( smiv_recs_[i]->smi_name().substr( 0 , mol_name.length() ) == mol_name ) {
        found_mol = true;
      }
      break;
    case 2 :
      if( string::npos != smiv_recs_[i]->smi_name().find( mol_name ) ) {
        found_mol = true;
      }
      break;
    }
    if( found_mol ) {
      mol_slider_->setValue( i );
      return true;
    }
  }

  return false;

}

// ****************************************************************************
pSmiVRec SmiVPanel::current_smiv_rec() const {

  if( smiv_recs_.empty() ) {
    return pSmiVRec();
  } else {
    return smiv_recs_[mol_slider_->value()];
  }

}

// ****************************************************************************
bool SmiVPanel::is_selected() const {

  return sel_box_->isChecked();

}

// ****************************************************************************
void SmiVPanel::set_selected( bool new_val ) {

  sel_box_->setChecked( new_val );

}

// ****************************************************************************
void SmiVPanel::set_title( const QString &new_title ) {

  title_->setText( new_title );

}

// ****************************************************************************
void SmiVPanel::keyPressEvent( QKeyEvent *event ) {

  if( event->modifiers() == Qt::ControlModifier ) {
    switch( event->key() ) {
    case Qt::Key_D :
      drop_current_mol();
      break;
    case Qt::Key_Z :
      undo_last_drop();
      break;
    case Qt::Key_N :
      change_current_mol( +1 );
      break;
    case Qt::Key_P :
      change_current_mol( -1 );
      break;
    case Qt::Key_E :
      go_to_last_mol();
      break;
    case Qt::Key_A :
      go_to_first_mol();
      break;
    default :
      break;
    }
  }

}

// ****************************************************************************
void SmiVPanel::build_widget() {

  QVBoxLayout *vbox = new QVBoxLayout;

  mol_disp_ = new DACLIB::QTMolDisplay2D;
  vbox->addWidget( mol_disp_ , 1 ); // add stretch, so it's this that grows, not the QLabels

  QHBoxLayout *hbox1 = new QHBoxLayout;
  hbox1->addWidget( new QLabel( "Input SMILES  ") );
  in_smi_ = new QLineEdit;
  in_smi_->setReadOnly( true );
  hbox1->addWidget( in_smi_ );
  vbox->addLayout( hbox1 );

  hbox1 = new QHBoxLayout;
  hbox1->addWidget( new QLabel( "Canon. SMILES ") );
  can_smi_ = new QLineEdit;
  can_smi_->setReadOnly( true );
  hbox1->addWidget( can_smi_ );
  vbox->addLayout( hbox1 );

  msg_ = new QLabel;
  vbox->addWidget( msg_ );

  hbox1 = new QHBoxLayout;
  hbox1->addLayout( vbox );

  mol_slider_ = new QSlider;
  mol_slider_->setInvertedAppearance( true );
  mol_slider_->setInvertedControls( true );
  mol_slider_->setSingleStep( 1 );
  mol_slider_->setPageStep( 1 );
  mol_slider_->setEnabled( false );
  hbox1->addWidget( mol_slider_ );

  connect( mol_slider_ , SIGNAL( valueChanged( int ) ) ,
           this , SLOT( slot_mol_slider_changed() ) );

  vbox = new QVBoxLayout;
  sel_box_ = new QCheckBox;
  title_ = new QLabel;
  title_->setAlignment( Qt::AlignHCenter );
  title_->setWordWrap( true );
  QHBoxLayout *hbox2 = new QHBoxLayout;
  hbox2->addWidget( sel_box_ );
  hbox2->addWidget( title_ , 1 );
  vbox->addLayout( hbox2 );

  connect( sel_box_ , SIGNAL( toggled( bool ) ) ,
           this , SLOT( slot_selection_box_changed() ) );

  vbox->addLayout( hbox1 );
  setLayout( vbox );

}

// ****************************************************************************
void SmiVPanel::colour_atoms() {

  if( sub_searches_.empty() ) {
    mol_disp_->clear_atom_colours();
    return;
  }

  // colouring by subsearch requires black and white molecule drawing
  mol_disp_->set_coloured_mol( false );
  mol_disp_->colour_atoms( sub_searches_ );

}

// ****************************************************************************
void SmiVPanel::drop_current_mol() {

  if( smiv_recs_.empty() ) {
    return;
  }

  int mol_num = mol_slider_->value();
  dropped_recs_.push_back( make_pair( smiv_recs_[mol_num] , mol_num ) );
  copy( smiv_recs_.begin() + mol_num + 1 , smiv_recs_.end() , smiv_recs_.begin() + mol_num );
  smiv_recs_.pop_back();
  mol_slider_->setMinimum( 0 );
  mol_slider_->setMaximum( int( smiv_recs_.size() - 1 ) );
  slot_mol_slider_changed();

}

// ****************************************************************************
void SmiVPanel::undo_last_drop() {

  if( dropped_recs_.empty() ) {
    return;
  }
  
  smiv_recs_.push_back( dropped_recs_.back().first );
  int ins_pos = dropped_recs_.back().second;
  mol_slider_->setMinimum( 0 );
  mol_slider_->setMaximum( int( smiv_recs_.size() - 1 ) );
  mol_slider_->setValue( ins_pos );
  if( smiv_recs_.size() > 1 ) {
    copy_backward( smiv_recs_.begin() + ins_pos , smiv_recs_.end() - 1 , smiv_recs_.end() );
    smiv_recs_[ins_pos] = dropped_recs_.back().first;
  }

  dropped_recs_.pop_back();
  slot_mol_slider_changed();

}

// ****************************************************************************
// move the mol_slider_ by the given step, if possible
void SmiVPanel::change_current_mol( int step ) {

  int new_val = mol_slider_->value() + step;
  if( new_val < 0 ) {
    new_val = 0;
  } else if( new_val > mol_slider_->maximum() ) {
    new_val = mol_slider_->maximum();
  }

  mol_slider_->setValue( new_val );

}

// ****************************************************************************
void SmiVPanel::write_smiles_to_stream( ostream &os ) const {

  for( int i = 0 , is = smiv_recs_.size() ; i < is ; ++i ) {
    os << smiv_recs_[i]->in_smi() << " " << smiv_recs_[i]->smi_name() << endl;
  }

}

// ****************************************************************************
void SmiVPanel::go_to_first_mol() {

  mol_slider_->setValue( 0 );

}

// ****************************************************************************
void SmiVPanel::go_to_last_mol() {

  mol_slider_->setValue( mol_slider_->maximum() );

}

// ****************************************************************************
void SmiVPanel::slot_mol_slider_changed() {

  if( smiv_recs_.empty() ) {
    mol_disp_->clear_display_molecule();
    in_smi_->setText( "" );
    can_smi_->setText( "" );
    msg_->setText( "No Molecules");
    return;
  }

  int mol_num = mol_slider_->value();
  in_smi_->setText( smiv_recs_[mol_num]->in_smi().c_str() );
  in_smi_->setCursorPosition( 0 );
  if( smiv_recs_[mol_num]->can_smi().empty() ) {
    smiv_recs_[mol_num]->create_can_smi();
  }
  can_smi_->setText( smiv_recs_[mol_num]->can_smi().c_str() );
  can_smi_->setCursorPosition( 0 );

  scoped_ptr<OEMolBase> mol( OENewMolBase( OEMolBaseType::OEDefault ) );
  OEParseSmiles( *mol , smiv_recs_[mol_num]->in_smi() );
  mol->SetTitle( smiv_recs_[mol_num]->smi_name() );
  mol_disp_->set_display_molecule( mol.get() );
  colour_atoms();

  QString msg = QString( "Displaying mol %1 of %2.").arg( mol_num + 1 ).arg( smiv_recs_.size() );
  msg_->setText( msg );

  emit new_display_mol( QString( smiv_recs_[mol_num]->smi_name().c_str() ) );

}

// ****************************************************************************
void SmiVPanel::slot_selection_box_changed() {

  emit( selection_box_changed( this ) );

}
