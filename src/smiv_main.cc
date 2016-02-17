//
// file smiv_main.cc
// David Cosgrove
// AstraZeneca
// 9th February 2009
//
// This is the main function for smiv, a little program that displays a
// file of molecules one at a time, and lets the user edit the SMILES,
// and write a canonical SMILES file out.

#include <string>

#include <QApplication>
#include <QMessageBox>

#include "SmiV.H"

using namespace std;

namespace DACLIB {
  bool check_oechem_licence( string &err_msg );
}

// *************************************************************************
int main( int argc , char **argv ) {

  QApplication a( argc , argv );
  SmiV *smiv = new SmiV;
  smiv->setGeometry( 100 , 100 , 500 , 500 );
  smiv->show();

  string err_msg;
  if( !DACLIB::check_oechem_licence( err_msg ) ) {
    QString msg = QString( "OEChem Licence error :\n%1\n").arg( err_msg.c_str() );
    QMessageBox::warning( smiv , "OEChem Licence Error" , msg );
    exit( 1 );
  }

  smiv->parse_args( argc , argv );

  return a.exec();

}
