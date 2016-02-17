//
// file SmiVFindMoleculeDialog.cc
// David Cosgrove
// AstraZeneca
// 15th January 2010
//

#include "SmiVFindMoleculeDialog.H"

#include <QComboBox>
#include <QFrame>
#include <QLayout>
#include <QLineEdit>
#include <QPushButton>

// ****************************************************************************
SmiVFindMoleculeDialog::SmiVFindMoleculeDialog( QWidget *parent ,
                                                Qt::WindowFlags f ) :
QDialog( parent , f ) {

  build_widget();

}

// ****************************************************************************
void SmiVFindMoleculeDialog::build_widget() {

  mol_name_ = new QLineEdit;
  search_mode_ = new QComboBox;
  search_mode_->addItem( "Exact Match" );
  search_mode_->addItem( "Starts With" );
  search_mode_->addItem( "Contains" );

  QVBoxLayout *vbox = new QVBoxLayout;
  vbox->addWidget( mol_name_ );
  vbox->addWidget( search_mode_ );
  vbox->addWidget( build_action_box() );

  setLayout( vbox );

}

// ****************************************************************************
QWidget *SmiVFindMoleculeDialog::build_action_box() {

  QFrame *action_frame = new QFrame;
  action_frame->setFrameStyle( QFrame::Box );

  QHBoxLayout *hlayout = new QHBoxLayout;

  QPushButton *button = new QPushButton( "Ok" );
  hlayout->addWidget( button );
  connect( button , SIGNAL( clicked() ) , this , SLOT( slot_ok_clicked() ) );

  button = new QPushButton( "Apply" );
  hlayout->addWidget( button );
  connect( button , SIGNAL( clicked() ) ,
           this , SLOT( slot_apply_clicked() ) );

  button = new QPushButton( "Cancel" );
  hlayout->addWidget( button );
  button->setDefault( false );
  button->setAutoDefault( false );
  connect( button , SIGNAL( clicked() ) , this , SLOT( reject() ) );

  action_frame->setLayout( hlayout );

  return action_frame;

}

// ****************************************************************************
void SmiVFindMoleculeDialog::slot_ok_clicked() {

  emit_search_details();
  accept();

}

// ****************************************************************************
void SmiVFindMoleculeDialog::slot_apply_clicked() {

  emit_search_details();

}

// ****************************************************************************
void SmiVFindMoleculeDialog::emit_search_details() {

  emit( molecule_search_name( mol_name_->text() ,
                              search_mode_->currentIndex() ) );

}

