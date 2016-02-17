//
// file QTSmilesEditDialog.cc
// Dave Cosgrove
// AstraZeneca
// 20th January 2010
//

#include <iostream>

#include "QTMolDisplay2D.H"
#include "QTSmilesEditDialog.H"

#include <QFrame>
#include <QLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QString>

#include <oechem.h>
#include <oedepict.h>

using namespace std;
using namespace OEChem;
using namespace OEDepict;

namespace DACLIB {

// *****************************************************************************
QTSmilesEditDialog::QTSmilesEditDialog( QWidget *parent ,
					Qt::WindowFlags f ) :
  QDialog( parent , f ) {

  build_widget();

}

// *****************************************************************************
void QTSmilesEditDialog::set_smiles( const QString &new_smi ) {

  smiles_box_->setText( new_smi );
  slot_update_mol_display();

}

// *****************************************************************************
QString QTSmilesEditDialog::get_smiles() const {

  return smiles_box_->text();

}

// *****************************************************************************
void QTSmilesEditDialog::set_name( const QString &new_smi ) {

  name_box_->setText( new_smi );

}

// *****************************************************************************
QString QTSmilesEditDialog::get_name() const {

  return name_box_->text();

}

// *****************************************************************************
void QTSmilesEditDialog::build_widget() {

  vbox_ = new QVBoxLayout;

  mol_disp_ = new DACLIB::QTMolDisplay2D;
  vbox_->addWidget( mol_disp_ , 1 ); // so it stretches properly
  connect( mol_disp_ , SIGNAL( atom_selected( unsigned int ) ) ,
           this , SLOT( slot_atom_selected( unsigned int ) ) );

  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addWidget( new QLabel( "SMILES : " ) );
  smiles_box_ = new QLineEdit;
  hbox->addWidget( smiles_box_ );
  vbox_->addLayout( hbox );

  connect( smiles_box_ , SIGNAL( returnPressed() ) ,
	   this , SLOT( slot_update_mol_display() ) );
  connect( smiles_box_ , SIGNAL( textChanged( const QString & ) ) ,
	   this , SLOT( slot_update_mol_display() ) );
  connect( smiles_box_ , SIGNAL( cursorPositionChanged( int , int ) ) ,
           this , SLOT( slot_smiles_cursor_moved( int , int ) ) );
  connect( smiles_box_ , SIGNAL( textChanged( const QString & ) ) ,
           this , SLOT( slot_data_changed() ) );

  hbox = new QHBoxLayout;
  hbox->addWidget( new QLabel( "Name : " ) );
  name_box_ = new QLineEdit;
  hbox->addWidget( name_box_ );
  vbox_->addLayout( hbox );

  connect( name_box_ , SIGNAL( textChanged( const QString & ) ) ,
           this , SLOT( slot_data_changed() ) );
  vbox_->addWidget( build_action_box() );

  setLayout( vbox_ );

  connect( this , SIGNAL( atom_selected( unsigned int )) ,
           mol_disp_ , SLOT( slot_atom_selected( unsigned int ) ) );

}

// *****************************************************************************
QWidget *QTSmilesEditDialog::build_action_box() {

  QFrame *action_frame = new QFrame;
  action_frame->setFrameStyle( QFrame::Box );

  QHBoxLayout *hlayout = new QHBoxLayout;

  QPushButton *button = new QPushButton( "Ok" );
  hlayout->addWidget( button );
  connect( button , SIGNAL( clicked() ) , this , SLOT( slot_ok_clicked() ) );

  apply_button_ = new QPushButton( "Apply" );
  hlayout->addWidget( apply_button_ );
  connect( apply_button_ , SIGNAL( clicked() ) , this , SLOT( slot_apply_clicked() ) );

  button = new QPushButton( "Cancel" );
  hlayout->addWidget( button );
  button->setDefault( false );
  button->setAutoDefault( false );
  connect( button , SIGNAL( clicked() ) , this , SLOT( reject() ) );

  action_frame->setLayout( hlayout );

  return action_frame;

}

// *****************************************************************************
void QTSmilesEditDialog::slot_ok_clicked() {

  if( name_box_->text().isEmpty() ) {
    QMessageBox::information( this , "Missing Molecule Name" , "You must give a name for the molecule" );
    return;
  }

  if( data_changed_ ) {
    emit new_smiles( smiles_box_->text() , name_box_->text() );
  }

  if( smiles_box_->text().isEmpty() ) {
    reject();
  } else {
    accept();
  }

}

// *****************************************************************************
void QTSmilesEditDialog::slot_apply_clicked() {

  if( name_box_->text().isEmpty() ) {
    QMessageBox::information( this , "Missing Molecule Name" , "You must give a name for the molecule" );
    return;
  }

  if( smiles_box_->text().isEmpty() ) {
    QMessageBox::information( this , "Missing SMILES" , "You must give a SMILES for the molecule" );
    return;
  }

  emit new_smiles( smiles_box_->text() , name_box_->text() );

  data_changed_ = false;

}

// *****************************************************************************
void QTSmilesEditDialog::slot_update_mol_display() {

  string smi = smiles_box_->text().toLocal8Bit().data();
  if( smi.empty() ) {
    mol_disp_->clear_display_molecule();
    return;
  }

  OEMol mol;
  OEParseSmiles( mol , smi );
  if( !name_box_->text().isEmpty() ) {
    mol.SetTitle( name_box_->text().toLocal8Bit().data() );
  }

  mol_disp_->set_display_molecule( &(mol.SCMol()) );

}

// *****************************************************************************
void QTSmilesEditDialog::slot_atom_selected( unsigned int oe_ind ) {

  // I assume this spell came from picto originally. It looks a tad fragile, to me,
  // in that it relies on oe_ind, the atom->GetIdx() being contiguous and counting from
  // zero in the same order as the SMILES string. I think that's a reasonable assumption
  // in this case, but it probably isn't generally reliable.
  QString smiles = smiles_box_->text();
  for( int i = 0 , is = smiles.length() ; i < is ; ++i ) {
    unsigned int cp = OEDepictSmilesAtomCount( smiles.left( i ).toLocal8Bit().data() );
    if( cp == oe_ind ) {
      smiles_box_->setCursorPosition( i );
      //smiles_box_->setSelection( i , 1 );;
      smiles_box_->setFocus( Qt::OtherFocusReason );
    }
  }
}

// *****************************************************************************
void QTSmilesEditDialog::slot_smiles_cursor_moved( int old_pos , int new_pos ) {

  if( new_pos ) {
    --new_pos;
  }

  QString bit_of_string = smiles_box_->text().left( new_pos );
  emit atom_selected( OEDepictSmilesAtomCount( bit_of_string.toLocal8Bit().data( ) ) );

}

// *****************************************************************************
void QTSmilesEditDialog::slot_data_changed() {

  data_changed_ = true;

}

} // EO namespace DACLIB


