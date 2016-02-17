//
// file QTSmartsEditDialog.cc
// Dave Cosgrove
// 31st October 2007
//

#include <iostream>

#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QString>
#include <QVBoxLayout>

#include <oechem.h>

#include "QTMolDisplay2D.H"
#include "QTSmartsEditDialog.H"

using namespace std;
using namespace OEChem;
using namespace OEPlatform;
using namespace OESystem;

// **************************************************************************
QTSmartsEditDialog::QTSmartsEditDialog( bool with_name_input ,
                                        QWidget *p ,
                                        Qt::WindowFlags f )
  : QDialog( p , f ) , smarts_name_box_( 0 ) ,
    with_name_input_( with_name_input ) {

  build_widget();

}

// *************************************************************************
void QTSmartsEditDialog::set_smarts( const QString &new_smarts ) {

  smarts_box_->setText( new_smarts );
  slot_display_smarts();

}

// *************************************************************************
void QTSmartsEditDialog::set_smarts_name( const QString &new_name ) {

  if( with_name_input_ && smarts_name_box_ ) {
    smarts_name_box_->setText( new_name );
  }

}

// ***************************************************************************
QString QTSmartsEditDialog::get_smarts_name() const {

  if( with_name_input_ && smarts_name_box_ )
    return smarts_name_box_->text();
  else
    return QString( "" );

}

// *************************************************************************
QString QTSmartsEditDialog::get_full_smarts() {

  string smt = smarts_box_->text().toLocal8Bit().data();
  SmartsLexReplace( smt , sub_defs_ );

  return smt.c_str();

}

// *************************************************************************
bool QTSmartsEditDialog::ok_to_close() {

  // must have a SMARTS definition in the box.
  if( smarts_box_->text().isEmpty() ) {
    QMessageBox::warning( this , "SMARTS Error" , "Need a SMARTS definition." );
    return false;
  }

  // same with the name, if we have one
  if( with_name_input_ && smarts_name_box_->text().isEmpty() ) {
    QMessageBox::warning( this , "SMARTS Error" , "Need a SMARTS name." );
    return false;
  }

  return true;

}

// *************************************************************************
void QTSmartsEditDialog::build_widget() {

  vbox_ = new QVBoxLayout;
  mol_disp_hbox_ = new QHBoxLayout;

  smarts_disp_ = new DACLIB::QTMolDisplay2D;
  mol_disp_hbox_->addWidget( smarts_disp_ );
  vbox_->addLayout( mol_disp_hbox_ , 1 ); /* set stretch, so it's the mol display that
                                             grows then the dialog grows vertically. */

  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addWidget( new QLabel( "SMARTS : " ) );
  smarts_box_ = new QLineEdit();
  hbox->addWidget( smarts_box_ );
  vbox_->addLayout( hbox );

  if( with_name_input_ ) {
    QHBoxLayout *hbox1 = new QHBoxLayout;
    hbox1->addWidget( new QLabel( "Name : " ) );
    smarts_name_box_ = new QLineEdit();
    hbox1->addWidget( smarts_name_box_ );
    vbox_->addLayout( hbox1 );
  }

  status_label_ = new QLabel();
  vbox_->addWidget( status_label_ );

  QWidget *action_box = build_action_box();
  vbox_->addWidget( action_box );

  setLayout( vbox_ );

  connect( smarts_box_ , SIGNAL( returnPressed() ) ,
           this , SLOT( slot_display_smarts() ) );
  if( with_name_input_ ) {
    connect( smarts_name_box_ , SIGNAL( returnPressed() ) ,
             this , SLOT( slot_display_smarts() ) );
  }

}

// *************************************************************************
QWidget *QTSmartsEditDialog::build_action_box() {

  QFrame *action_frame = new QFrame;
  action_frame->setFrameStyle( QFrame::Box );

  QHBoxLayout *hlayout = new QHBoxLayout;

  QPushButton *button = new QPushButton( "Ok" );
  hlayout->addWidget( button );
  // don't have a default button, as we want return to force a read of the SMARTS
  // typed in for display etc.
  button->setDefault( false );
  button->setAutoDefault( false );
  connect( button , SIGNAL( clicked() ) , this , SLOT( slot_ok_clicked() ) );

  button = new QPushButton( "Cancel" );
  hlayout->addWidget( button );
  button->setDefault( false );
  button->setAutoDefault( false );
  connect( button , SIGNAL( clicked() ) , this , SLOT( reject() ) );

  action_frame->setLayout( hlayout );

  return action_frame;

}

// *********************************************************************
void QTSmartsEditDialog::slot_display_smarts() {

  // There has been a change in OEDepict v2.0 which breaks the code below, such that
  // it frequently throws a Fatal error due to a bond style error.  For now, I'm just
  // putting the SMARTS through as a SMILES string which gives a different, not quite
  // as adequate, drawing.
  QString q_smarts = get_full_smarts();
  string smarts = q_smarts.toLocal8Bit().data();

#ifdef NOTYET
  // it seems like the bug's been fixed
  if( smarts.empty() ) {
    status_label_->setText( "" );
    smarts_disp_->clear_display_molecule();
    return;
  }
  OEMolBase *mol = OENewMolBase();
  OEParseSmiles( *mol , smarts );
  smarts_disp_->set_display_molecule( mol );
#endif

  OEPlatform::oeosstream oeerrs;
  OEThrow.SetOutputStream( oeerrs );
  OESubSearch subs( smarts.c_str() );
  OEThrow.SetOutputStream( oeerr );
  string errstr = oeerrs.str();
  if( !errstr.empty() ) {
    QMessageBox::warning( this , "SMARTS Error" , QString( errstr.c_str() ) );
    status_label_->setText( "Bad SMARTS" );
    smarts_disp_->clear_display_molecule();
    return;
  }

  OEQMol qmol( subs.GetPattern() );
  if( !qmol ) {
    status_label_->setText( "Bad SMARTS" );
    smarts_disp_->clear_display_molecule();
    return;
  }

  status_label_->setText( "" );
  // seems to be needed to get atoms displayed correctly by atomic symbol.
  // Can't find it documented anywhere.
  OEDisassembleExpressions( qmol );
  smarts_disp_->set_display_molecule( &(qmol.SCMol()) );

}

// *********************************************************************
void QTSmartsEditDialog::slot_ok_clicked() {

  if( ok_to_close() ) {
    accept();
  } else {
    return;
  }

}
