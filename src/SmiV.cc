//
// file SmiV.cc
// David Cosgrove
// AstraZeneca
// 9th February 2009.

#include "SmiV.H"
#include "SmivDataTable.H"
#include "SmiVFindMoleculeDialog.H"
#include "SmiVPanel.H"
#include "SmiVRecord.H"
#include "SmiVSettings.H"

#include "DACOEMolAtomIndex.H"
#include "SMARTSExceptions.H"
#include "QT4SelectItems.H"
#include "QTSmartsEditDialog.H"
#include "QTSmartsIntPickDialog.H"
#include "QTSmilesEditDialog.H"
#include "stddefs.H"

#include <QAction>
#include <QApplication>
#include <QFileDialog>
#include <QFileInfo>
#include <QHeaderView>
#include <QInputDialog>
#include <QLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QPushButton>
#include <QSlider>
#include <QSplitter>
#include <QStatusBar>
#include <QTableView>

#include <oechem.h>

#include <fstream>
#include <iostream>
#include <iterator>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

extern string BUILD_TIME; // put together in build_time.cc

namespace DACLIB {
  void apply_daylight_aromatic_model( OEMolBase &mol );
  void read_smarts_file( const string &smarts_file ,
                         vector<pair<string,string> > &input_smarts ,
                         vector<pair<string,string> > &smarts_sub_defn );
  OESubSearch *create_eosubsearch( const string &smarts ,
                                   const string &smarts_name ,
                                   bool reorder ,
                                   vector<pair<string,string> > &sub_defns );
}

// ****************************************************************************
SmiV::SmiV() : QMainWindow() , find_mol_dialog_( 0 ) , smiles_edit_dialog_( 0 ) {

  build_actions();
  build_menubar();
  build_widget();

  last_dir_ = QString( "." );

  // turn off reporting of errors
  OEThrow.SetOutputStream( OEPlatform::oenul );

}

// ****************************************************************************
SmiV::~SmiV() {

}

// ****************************************************************************
void SmiV::parse_args( int argc , char **argv ) {

  SmiVSettings ss( argc , argv );
  if( !ss ) {
    return; // it's not the end of the world
  }

  if( !ss.mol_file().empty() ) {
    read_mol_file( QString( ss.mol_file().c_str() ) );
  }
  if( !ss.smarts_file().empty() ) {
    read_smarts_file( QString( ss.smarts_file().c_str() ));
  }
  if( !ss.mdl_file().empty() ) {
    read_mdl_query_file( QString( ss.mdl_file().c_str() ) );
  }
  if( !ss.data_file().empty() ) {
    read_data_file( QString( ss.data_file().c_str() ) );
  }
  usage_text_ = ss.usage_text();

}

// ****************************************************************************
void SmiV::slot_quit() {

  exit( 0 );

}

// ****************************************************************************
void SmiV::slot_file_read_mols() {

  QString filename =
      QFileDialog::getOpenFileName( this , "Choose molecule file" , last_dir_ ,
                                    "Molecules (*.smi *.mol2 *.sdf *.oeb)" );
  if( filename.isEmpty() ) {
    return;
  }

  read_mol_file( filename );

}

// ****************************************************************************
void SmiV::slot_file_reread_mols() {

  if( !last_mol_file_.isEmpty() ) {
    smiv_recs_.clear();
    read_mol_file( last_mol_file_ );
  } else {
    QMessageBox::warning( this , "No Mol file" , "You have not yet read a molecule file to re-read." );
  }
}

// ****************************************************************************
void SmiV::slot_file_read_smarts() {

  QString filename =
      QFileDialog::getOpenFileName( this , "Choose SMARTS file" , last_dir_ ,
                                    "SMARTS (*.smt);;Any File (*)" );
  if( filename.isEmpty() )
    return;

  read_smarts_file( filename );

}

// ****************************************************************************
void SmiV::slot_file_reread_smarts() {

  if( !last_smarts_file_.isEmpty() ) {
    smarts_.clear();
    smarts_sub_defn_.clear();
    read_smarts_file( last_smarts_file_ );
  } else {
    QMessageBox::warning( this , "No SMARTS file" , "You have not yet read a SMARTS file to re-read." );
  }

}

// ****************************************************************************
void SmiV::slot_file_read_mdl_query() {

  QString filename =
      QFileDialog::getOpenFileName( this , "Choose MDL Query file" , last_dir_ ,
                                    "MDL (*.mol)" );
  if( filename.isEmpty() )
    return;

  read_mdl_query_file( filename );

}

// ****************************************************************************
void SmiV::slot_file_read_data() {

  QString filename =
      QFileDialog::getOpenFileName( this , "Choose data file" , last_dir_ ,
                                    "CSV (*.csv);;Any File(*)" );
  if( filename.isEmpty() )
    return;

  read_data_file( filename );

}


// ****************************************************************************
void SmiV::slot_write_smiles() {

  QString filename =
    QFileDialog::getSaveFileName( this , "SMILES file" ,
                                  last_dir_ , "SMILES file (*.smi)" );

  if( filename.isEmpty() )
    return;

  write_smiles_file( filename );

}

// ****************************************************************************
void SmiV::slot_clear_molecules() {

  smiv_recs_.clear();
  left_panel_->add_data( smiv_recs_ );
  right_panel_->add_data( smiv_recs_ );
  right_panel_->hide();

}

// ****************************************************************************
void SmiV::slot_write_smarts() {

  QString filename =
    QFileDialog::getSaveFileName( this , "SMARTS file" ,
                                  last_dir_ , "SMARTS file (*.smt)" );

  if( filename.isEmpty() )
    return;

  write_smarts_file( filename );

}

// ****************************************************************************
void SmiV::slot_clear_smarts() {

  smarts_.clear();
  smarts_sub_defn_.clear();

}

// ****************************************************************************
void SmiV::slot_find_mol() {

  if( !find_mol_dialog_ ) {
    find_mol_dialog_ = new SmiVFindMoleculeDialog( this );
    connect( find_mol_dialog_ , SIGNAL( molecule_search_name( QString , int ) ) ,
             this , SLOT( slot_molecule_search_name( QString , int ) ) );
  }

  find_mol_dialog_->show();

}

// *****************************************************************************
void SmiV::slot_smarts_match() {

  if( smarts_.empty() ) {
    QMessageBox::information( this , "SMARTS Match" , "No SMARTS defined." );
    return;
  }

  do_smarts_matching();

}

// *****************************************************************************
void SmiV::slot_smarts_edit() {

  if( smarts_.empty() ) {
    slot_smarts_input_int_pick();
    return;
  }

  // first, select SMARTS pattern to edit
  vector<char> sel_smarts;
  get_query_to_use( sel_smarts , smarts_ , true );
  if( 0 == count( sel_smarts.begin() , sel_smarts.end() , 1 ) ) {
    return;
  }
  vector<char>::iterator p = std::find( sel_smarts.begin() , sel_smarts.end() , 1 );

  string int_pick_smi = get_active_smiles();

  QTSmartsIntPickDialog sipd( int_pick_smi , this );
  sipd.set_smarts_sub_defns( smarts_sub_defn_ );
  sipd.set_smarts( smarts_[distance( sel_smarts.begin() , p )].second.c_str() );
  sipd.set_smarts_name(smarts_[distance( sel_smarts.begin() , p )].first.c_str() );
  if( QDialog::Accepted != sipd.exec() ) {
    return;
  }

  add_smarts_definition( sipd );

}

// *****************************************************************************
void SmiV::slot_smarts_input_int_pick() {

  string int_pick_smi = get_active_smiles();

  QTSmartsIntPickDialog sipd( int_pick_smi , this );
  sipd.set_smarts_sub_defns( smarts_sub_defn_ );
  if( QDialog::Accepted != sipd.exec() ) {
    return;
  }

  add_smarts_definition( sipd );

}

// *****************************************************************************
void SmiV::slot_mdl_query_match() {

  if( mdl_queries_.empty() ) {
    QMessageBox::information( this , "MDL Query Match" , "No MDL query defined." );
    return;
  }

  do_mdl_query_matching();

}

// *****************************************************************************
void SmiV::slot_show_about_box() {

  QString msg;
  msg = QString( "SmiV\n"
                 "Copyright AstraZeneca 2016\n"
                 "Built on %1\nusing OEChem version %2\nand Qt version %3\n" )
    .arg( BUILD_TIME.c_str() ).arg( OEChemGetRelease() ).arg( QT_VERSION_STR );
  msg += usage_text_.c_str();

  QMessageBox::information( this , "About SmiV" , msg ,
                            QMessageBox::Ok | QMessageBox::Default );

}

// ****************************************************************************
// search mode is 0 for Exact Match, 1 for Starts With, 2 for Contains
void SmiV::slot_molecule_search_name( QString search_name , int search_mode ) {

#ifdef NOTYET
  cout << "Looking for " << search_name.toStdString() << " in mode " << search_mode
      << endl;
#endif

  SmiVPanel *panel = get_active_panel();
  bool found_in_panel = panel->show_molecule( search_name.toLocal8Bit().data() , search_mode );
  if( !found_in_panel ) {
    panel = get_inactive_panel();
    if( !panel->isHidden() ) {
      panel->show_molecule( search_name.toLocal8Bit().data() , search_mode );
    }
  }

}

// ****************************************************************************
void SmiV::slot_panel_selection_changed( QWidget *panel ) {

  if( panel == left_panel_ ) {
    right_panel_->set_selected( !left_panel_->is_selected() );
  } else {
    left_panel_->set_selected( !right_panel_->is_selected() );
  }

}

// ****************************************************************************
void SmiV::slot_input_smiles() {

#ifdef NOTYET
  cout << "SmiV::slot_input_smiles" << endl;
#endif

  if( !smiles_edit_dialog_ ) {
    smiles_edit_dialog_ = new DACLIB::QTSmilesEditDialog( this );
    connect( smiles_edit_dialog_ , SIGNAL( new_smiles( QString , QString ) ) ,
             this , SLOT( slot_new_smiles( QString , QString ) ) );
  }

  if( QDialog::Accepted != smiles_edit_dialog_->exec() ) {
    return;
  }

}

// ****************************************************************************
void SmiV::slot_edit_smiles() {

  if( !smiles_edit_dialog_ ) {
    smiles_edit_dialog_ = new DACLIB::QTSmilesEditDialog( this );
    connect( smiles_edit_dialog_ , SIGNAL( new_smiles( QString , QString ) ) ,
             this , SLOT( slot_new_smiles( QString , QString ) ) );
  }

  SmiVPanel *sp = get_active_panel();
  pSmiVRec smiv_rec = sp->current_smiv_rec();
  if( smiv_rec ) {
    smiles_edit_dialog_->set_smiles( smiv_rec->in_smi().c_str() );
    smiles_edit_dialog_->set_name( smiv_rec->smi_name().c_str() );
  }

  if( QDialog::Accepted != smiles_edit_dialog_->exec() ) {
    return;
  }

}

// ****************************************************************************
void SmiV::slot_new_smiles( QString new_smiles , QString new_name ) {

  cout << "slot_new_smiles : " << new_smiles.toStdString() << " : " << new_name.toStdString() << endl;
  update_smiv_recs( new_smiles.toLocal8Bit().data() ,
                    new_name.toLocal8Bit().data() );

}

// ****************************************************************************
void SmiV::slot_full_list() {

  right_panel_->hide();
  left_panel_->add_data( smiv_recs_ );
  left_panel_->set_title( "All Molecules" );

}

// ****************************************************************************
void SmiV::slot_new_mol_list() {

  new_mol_list( QString( "" ) );

}

// ****************************************************************************
void SmiV::slot_save_mol_list() {

  // just writes the names
  QString filename = QFileDialog::getSaveFileName( this , "What filename?" , last_dir_);
  if( filename.isEmpty() ) {
    return;
  }
  ofstream ofs( filename.toLocal8Bit().data() );
  if( !ofs ) {
    QString msg = QString( "Couldn't open %1 for writing.").arg( filename );
    QMessageBox::warning( this , "Bad File" , msg );
    return;
  }

  SmiVPanel *sp = get_active_panel();
  vector<pSmiVRec> these_recs = sp->smiv_recs();
  transform( these_recs.begin() , these_recs.end() ,
             ostream_iterator<string>( ofs , "\n" ) , bind( &SmiVRecord::smi_name , _1 ) );

}

// ****************************************************************************
void SmiV::slot_show_mol_list() {

  QAction *action = qobject_cast<QAction *>(sender());
  if( action ) {
    show_mol_list( action->text() );
  }

}

// ****************************************************************************
void SmiV::slot_sort_data_table( int col_num ) {

  if( Qt::DescendingOrder == data_table_->last_sort_order() ) {
    data_table_->sort( col_num , Qt::AscendingOrder );
  } else {
    data_table_->sort( col_num , Qt::DescendingOrder );
  }

}

// ****************************************************************************
void SmiV::slot_data_table_cell_double_clicked( const QModelIndex &ind ) {

  int row_num( ind.row() );
  int col_num( ind.column() );

  if( row_num < 0 || row_num >= data_table_->rowCount() ||
      col_num < 0 || col_num >= data_table_->columnCount() ) {
    return; // not a useful cell
  }

  if( col_num == get_data_table_column_number( QString( "CoreSmiles" ) ) ) {
    show_core_smiles_molecule( row_num );
  } else if( col_num == get_data_table_column_number( QString( "CoreSmarts" ) ) ) {
    search_with_core_smarts( row_num );
  } else {
    show_molecule_from_data_table( row_num );
  }

}

// *****************************************************************************************
void SmiV::slot_data_table_show_row( QString row_name ) {

  // assume the compound name is in the first column of the table
  for( unsigned int i = 0 , is = data_table_->rowCount() ; i < is ; ++i ) {
    if( data_table_->data( i , 0 ).toString() == row_name ) {
      data_table_view_->scrollTo( data_table_->index( i , 0 ) );
      break;
    }
  }

}

// ****************************************************************************
void SmiV::build_actions() {

  build_file_actions();
  build_smarts_actions();
  build_mdl_query_actions();
  build_molecule_actions();
  build_help_actions();

}

// **************************************************************************
void SmiV::build_file_actions() {

  file_read_data_ = new QAction( "Read Data File" , this );
  file_read_data_->setStatusTip( "Read arbitrary data file into table" );
  connect( file_read_data_ , SIGNAL( triggered() ) , this , SLOT( slot_file_read_data() ) );

  file_quit_ = new QAction( "Quit" , this );
  file_quit_->setShortcut( QString( "Ctrl+Q" ) );
  file_quit_->setStatusTip( "Exit the application" );
  connect( file_quit_ , SIGNAL( triggered() ) , this , SLOT( slot_quit() ) );

}

// **************************************************************************
void SmiV::build_smarts_actions() {

  file_read_smarts_ = new QAction( "Read SMARTS File" , this );
  connect( file_read_smarts_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_file_read_smarts() ) );

  file_reread_smarts_ = new QAction( "Re-Read last SMARTS File" , this );
  connect( file_reread_smarts_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_file_reread_smarts() ) );

  smarts_match_ = new QAction( "Match" , this );
  connect( smarts_match_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_smarts_match() ) );
  smarts_match_->setShortcut( QString( "Ctrl+M" ) );

  smarts_input_edit_ = new QAction( "Edit Existing" , this );
  connect( smarts_input_edit_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_smarts_edit() ) );

  smarts_input_int_pick_ = new QAction( "Interactive Pick" , this );
  connect( smarts_input_int_pick_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_smarts_input_int_pick() ) );

  smarts_write_ = new QAction( "Write File" , this );
  connect( smarts_write_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_write_smarts() ) );

  clear_smarts_ = new QAction( "Clear SMARTS" , this );
  connect( clear_smarts_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_clear_smarts() ) );

}

// **************************************************************************
void SmiV::build_mdl_query_actions() {

  file_read_mdl_query_ = new QAction( "Read MDL Query File" , this );
  connect( file_read_mdl_query_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_file_read_mdl_query() ) );

  mdl_query_match_ = new QAction( "Match" , this );
  connect( mdl_query_match_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_mdl_query_match() ) );

}

// **************************************************************************
void SmiV::build_molecule_actions() {

  file_read_mol_ = new QAction( "Read Molecule File" , this );
  file_read_mol_->setShortcut( QString( "Ctrl+R" ) );
  connect( file_read_mol_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_file_read_mols() ) );
  file_reread_mol_ = new QAction( "Re-Read Last Molecule File" , this );
  connect( file_reread_mol_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_file_reread_mols() ) );

  file_write_smiles_ = new QAction( "Write SMILES File" , this );
  file_write_smiles_->setShortcut( QString( "Ctrl+W" ) );
  connect( file_write_smiles_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_write_smiles() ) );

  find_mol_ = new QAction( "Find Mol." , this );
  find_mol_->setShortcut( QString( "Ctrl+F" ) );
  connect( find_mol_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_find_mol() ) );

  input_smiles_ = new QAction( "Input SMILES" , this );
  connect( input_smiles_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_input_smiles() ) );

  edit_smiles_ = new QAction( "Edit SMILES" , this );
  connect( edit_smiles_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_edit_smiles() ) );

  new_list_ = new QAction( "Create New" , this );
  new_list_->setShortcut( QString( "Ctrl+L") );
  connect( new_list_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_new_mol_list() ) );

  save_list_ = new QAction( "Save" , this );
  connect( save_list_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_save_mol_list() ) );

  full_list_ = new QAction( "Full Set" , this );
  full_list_->setShortcut( QString( "Ctrl+A") );
  connect( full_list_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_full_list() ) );

  clear_mols_ = new QAction( "Clear Molecules" , this );
  connect( clear_mols_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_clear_molecules() ) );

}

// **************************************************************************
void SmiV::build_help_actions() {

  help_show_about_ = new QAction( "About SmiV" , this );
  connect( help_show_about_ , SIGNAL( triggered() ) ,
	   this , SLOT( slot_show_about_box() ) );

}

// ****************************************************************************
void SmiV::build_menubar() {

  QMenu *file_menu = menuBar()->addMenu( "File" );
  file_menu->addAction( file_read_mol_ );
  file_menu->addAction( file_reread_mol_ );
  file_menu->addAction( file_read_smarts_ );
  file_menu->addAction( file_reread_smarts_ );
  file_menu->addAction( file_read_mdl_query_ );
  file_menu->addAction( file_read_data_ );
  file_menu->addSeparator();
  file_menu->addAction( file_quit_ );

  QMenu *smarts_menu = menuBar()->addMenu( "SMARTS" );
  smarts_menu->addAction( smarts_match_ );
  smarts_menu->addAction( smarts_input_edit_ );
  smarts_menu->addAction( smarts_input_int_pick_ );
  smarts_menu->addAction( file_read_smarts_ );
  smarts_menu->addAction( file_reread_smarts_ );
  smarts_menu->addAction( smarts_write_ );
  smarts_menu->addAction( clear_smarts_ );

  QMenu *mdl_query_menu = menuBar()->addMenu( "MDL Query" );
  mdl_query_menu->addAction( mdl_query_match_ );

  QMenu *mol_menu = menuBar()->addMenu( "Molecules" );
  mol_menu->addAction( file_read_mol_ );
  mol_menu->addAction( file_reread_mol_ );
  mol_menu->addAction( file_write_smiles_ );
  mol_menu->addAction( find_mol_ );
  mol_menu->addAction( input_smiles_ );
  mol_menu->addAction( edit_smiles_ );
  mol_menu->addAction( clear_mols_ );
  mol_lists_menu_ = mol_menu->addMenu( "Lists");
  mol_lists_menu_->addAction( full_list_ );
  mol_list_separator_ = mol_lists_menu_->addSeparator(); // so we can insert ahead of it in menu
  mol_lists_menu_->addAction( new_list_ );
  mol_lists_menu_->addAction( save_list_ );

  QMenu *help_menu = menuBar()->addMenu( "Help" );
  help_menu->addAction( help_show_about_ );

}

// ****************************************************************************
void SmiV::build_widget() {

  QSplitter *splitter = new QSplitter;
  splitter->setOrientation( Qt::Vertical );

  QHBoxLayout *hbox = new QHBoxLayout;

  left_panel_ = new SmiVPanel;
  left_panel_->set_selected( true );
  hbox->addWidget( left_panel_ );

  connect( left_panel_ , SIGNAL( selection_box_changed( QWidget * ) ) ,
           this , SLOT( slot_panel_selection_changed( QWidget * ) ) );
  connect( left_panel_ , SIGNAL( new_display_mol( QString ) ) ,
           this , SLOT( slot_data_table_show_row( QString ) ) );

  right_panel_ = new SmiVPanel;
  hbox->addWidget( right_panel_ );
  right_panel_->hide();

  connect( right_panel_ , SIGNAL( selection_box_changed( QWidget * ) ) ,
           this , SLOT( slot_panel_selection_changed( QWidget * ) ) );
  connect( right_panel_ , SIGNAL( new_display_mol( QString ) ) ,
           this , SLOT( slot_data_table_show_row( QString ) ) );

  data_table_ = new SmivDataTable;
  data_table_view_ = new QTableView;
  data_table_view_->setModel( data_table_ );
  data_table_view_->hide();
  connect( data_table_view_->horizontalHeader() , SIGNAL( sectionClicked( int ) ) ,
           this , SLOT( slot_sort_data_table( int ) ) );

  connect( data_table_view_ , SIGNAL( doubleClicked( const QModelIndex & ) ) ,
           this , SLOT( slot_data_table_cell_double_clicked( const QModelIndex & ) ) );

  QWidget *hbox_wid = new QWidget;
  hbox_wid->setLayout( hbox );
  splitter->addWidget( hbox_wid );
  splitter->addWidget( data_table_view_ );
  setCentralWidget( splitter );

}

// ****************************************************************************
void SmiV::read_mol_file( const QString &filename ) {

  QFileInfo fi( filename );
  if( !fi.exists() ) {
    return;
  }

  last_mol_file_ = filename;
  last_dir_ = fi.absolutePath();

  if( filename.endsWith( ".smi") || filename.endsWith( ".smi.gz") ) {
    read_smiles_file( filename );
  } else {
    read_other_mol_file( filename );
  }

  show_all_molecules();
  update_status_count();

}

// *****************************************************************************
void SmiV::read_smiles_file( const QString &filename ) {

  QFileInfo fi( filename );
  if( !fi.exists() ) {
    return;
  }

  last_dir_ = fi.absolutePath();

  ifstream file( filename.toLocal8Bit().data() , ios_base::in | ios_base::binary );
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  if( filename.endsWith( ".smi.gz" ) ) {
    in.push( boost::iostreams::gzip_decompressor() );
  }
  in.push( file );

  istreambuf_iterator<char> init( &in ) , eos;

  vector<char> nextline;
  vector<pair<string,string> > in_smiles;

  while( init != eos ) {
    nextline.clear();
    while( init != eos && *init != '\n' ) {
      nextline.push_back( *init );
      ++init;
    }
    boost::trim( nextline );
    if( !nextline.empty() ) {
      vector<string> split_line;
      split( split_line , nextline , is_any_of( " ,\t" ) );
      if( !split_line.empty() ) {
        string mol_name;
        if( split_line.size() > 1 ) {
          mol_name = string( nextline.begin() + split_line[0].length() + 1 ,
                             nextline.end() );
          boost::trim( mol_name );
        } else {
          mol_name = string( "Mol" ) + lexical_cast<string>( smiv_recs_.size() + 1 );
        }
        SmiVRecord *smiv_rec = new SmiVRecord( split_line[0] , mol_name );
        smiv_recs_.push_back( pSmiVRec( smiv_rec ) );
      }
    }
    ++init; // get past '\n'
  }

}

// ****************************************************************************
void SmiV::read_other_mol_file( const QString &filename ) {

  QFileInfo fi( filename );
  if( !fi.exists() ) {
    return;
  }

  last_dir_ = fi.absolutePath();

  oemolistream ims( filename.toLocal8Bit().data() );
  OEMol mol;
  while( ims >> mol ) {
    DACLIB::apply_daylight_aromatic_model( mol );
    SmiVRecord *smiv_rec = new SmiVRecord( mol );
    smiv_recs_.push_back( pSmiVRec( smiv_rec ) );
  }

}

// ****************************************************************************
void SmiV::read_smarts_file( const QString &filename ) {

  QFileInfo fi( filename );
  if( !fi.exists() ) {
    return;
  }

  last_smarts_file_ = filename;
  last_dir_ = fi.absolutePath();

  try {
    DACLIB::read_smarts_file( filename.toLocal8Bit().data() ,
                              smarts_ , smarts_sub_defn_ );
  } catch( DACLIB::SMARTSSubDefnError &e ) {
    cout << e.what() << endl;
    QMessageBox::warning( this , "SMARTS file error" , e.what() );
    return;
  } catch( DACLIB::SMARTSFileError &e ) {
    cout << e.what() << endl;
    QMessageBox::warning( this , "SMARTS file error" , e.what() );
    return;
  }

  update_status_count();

}

// ****************************************************************************
void SmiV::read_mdl_query_file( const QString &filename ) {

  cout << "Reading MDL query file " << filename.toLocal8Bit().data() << endl;

  QFileInfo fi( filename );
  if( !fi.exists() ) {
    return;
  }

  last_dir_ = fi.absolutePath();

  ifstream ifs( filename.toLocal8Bit().data() );
  if( !ifs || !ifs.good() ) {
    QMessageBox::warning( this , "MDL query file error" ,
                          QString( "Couldn't open %1 for reading.").arg( filename ) );
    return;
  }

  int count = 0;
  string next_line;
  vector<string> next_query;
  while( 1 ) {
    getline( ifs , next_line );
    if( ifs.eof() || !ifs.good() ) {
      ++count;
      store_next_mdl_query( next_query , filename , count );
      break;
    }
    if( next_line == string( "$$$$") ) {
      ++count;
      store_next_mdl_query( next_query , filename , count );
      next_query.clear();
    } else {
      next_query.push_back( next_line );
    }

  }

  update_status_count();

}

// ****************************************************************************
void SmiV::read_data_file( const QString &filename ) {

  QFileInfo fi( filename );
  if( !fi.exists() ) {
    return;
  }

  last_dir_ = fi.absolutePath();

  data_table_->read_data_from_file( this , filename );

  statusBar()->showMessage( QString( "Read %1 lines of data." ).arg( data_table_->rowCount() ) , 2000 );

  data_table_view_->show();

}

// ****************************************************************************
void SmiV::write_smiles_file( const QString &filename ) {

  QFileInfo fi( filename );
  last_dir_ = fi.absolutePath();

  ofstream ofs( filename.toLocal8Bit().data() );
  if( !ofs.good() ) {
    QMessageBox::warning( this , "File Open Error" ,
			  QString( "Couldn't open file \n%1\nfor writing." ).arg( filename ) );
    return;
  }

  get_active_panel()->write_smiles_to_stream( ofs );

}

// ****************************************************************************
void SmiV::write_smarts_file( const QString &filename ) {

  QFileInfo fi( filename );
  last_dir_ = fi.absolutePath();

  ofstream ofs( filename.toLocal8Bit().data() );
  if( !ofs.good() ) {
    QMessageBox::warning( this , "File Open Error" ,
              QString( "Couldn't open file \n%1\nfor writing." ).arg( filename ) );
    return;
  }

  ofs << "#" << endl
      << "# SMARTS file written by SmiV" << endl
      << "#" << endl
      << "# Full Definitions" << endl
      << "#" << endl;

  for( int i = 0 , is = smarts_.size() ; i < is ; ++i ) {
    ofs << smarts_[i].first << "\t" << smarts_[i].second << "\t" << "1\t1" << endl;
  }

  ofs << "#" << endl
      << "# Sub-Definitions (Vector Bindings)" << endl
      << "#" << endl;

  for( int i = 0 , is = smarts_sub_defn_.size() ; i < is ; ++i ) {
    // don't write the sub def if it's also a SMARTS def, as we don't want things
    // in the file twice
    if( smarts_.end() == std::find_if( smarts_.begin() , smarts_.end() ,
                                       bind( std::equal_to<string>() ,
                                                    bind( &pair<string,string>::first , _1 ) ,
                                                    smarts_sub_defn_[i].first ) ) ) {
      ofs << smarts_sub_defn_[i].first << "\t" << smarts_sub_defn_[i].second << "\t" << "1\t0" << endl;
    }
  }

}

// ****************************************************************************
void SmiV::update_status_count() {

  QString msg = QString( "Now have %1 active molecules, %2 SMARTS definitions, %3 MDL queries." )
    .arg( smiv_recs_.size() ).arg( smarts_.size() ).arg( mdl_queries_.size() );
  statusBar()->showMessage( msg , 0 );

}

// ****************************************************************************
// reset to just left_panel_, showing all molecules
void SmiV::show_all_molecules() {

  right_panel_->hide();
  left_panel_->add_data( smiv_recs_ );
  left_panel_->set_title( QString( "All Molecules" ) );

}

// ****************************************************************************
// reset to just left_panel_, showing all molecules
void SmiV::show_mol_list( const QString &list_name ) {

  vector<pair<string,vector<pSmiVRec> > >::iterator p =
      find_if( rec_lists_.begin() , rec_lists_.end() ,
               bind( equal_to<string>() ,
                     bind( &pair<string,vector<pSmiVRec> >::first , _1 ) ,
                     list_name.toLocal8Bit().data() ) );
  right_panel_->hide();
  left_panel_->add_data( p->second );
  left_panel_->set_title( list_name );

}

// ****************************************************************************
// get name of molecule in row_num of data table and show whichever panel currently has it
void SmiV::show_molecule_from_data_table( int row_num ) {

  slot_molecule_search_name( data_table_->data( row_num , 0 ).toString() , 0 );

}

// ****************************************************************************
// do smarts matching on contents of left panel only
void SmiV::do_smarts_matching() {

  // first, select one or more SMARTS patterns to use
  vector<char> sel_smarts;
  get_query_to_use( sel_smarts , smarts_ , false );
  if( 0 == count( sel_smarts.begin() , sel_smarts.end() , 1 ) ) {
    show_all_molecules();
    return;
  }

  vector<pair<boost::shared_ptr<OESubSearch>,string> > sub_searches;
  QString smarts_list;
  build_sub_searches_from_smarts( sel_smarts , smarts_list , sub_searches );
  do_substructure_matching( sub_searches , smarts_list , true );

}

// ****************************************************************************
// do MDL query matching on contents of left panel only
void SmiV::do_mdl_query_matching() {

  cout << "SmiV::do_mdl_query_matching()" << endl;
  vector<char> sel_query;
  get_query_to_use( sel_query , mdl_queries_ , false );
  if( 0 == count( sel_query.begin() , sel_query.end() , 1 ) ) {
    show_all_molecules();
    return;
  }

  vector<pair<boost::shared_ptr<OESubSearch>,string> > sub_searches;
  QString query_list;
  build_sub_searches_from_mdl_queries( sel_query , query_list , sub_searches );
  do_substructure_matching( sub_searches , query_list , true );

}

// ****************************************************************************
void SmiV::do_substructure_matching( vector<pair<boost::shared_ptr<OESubSearch>,string> > &sub_searches ,
                                     const QString &list_name , bool show_non_matches ) {

  vector<pSmiVRec> ones_to_do = left_panel_->smiv_recs();

  QApplication::setOverrideCursor( Qt::WaitCursor );
  vector<pSmiVRec> left_list , right_list;
  for( int i = 0 , is = ones_to_do.size() ; i < is ; ++i ) {
    scoped_ptr<OEMolBase> mol( OENewMolBase( OEMolBaseType::OEDefault ) );
    OEParseSmiles( *mol ,  ones_to_do[i]->in_smi() );
    DACLIB::apply_daylight_aromatic_model( *mol );
    bool match_found( false );
    for( int j = 0 , js = sub_searches.size() ; j < js ; ++j ) {
      if( sub_searches[j].first->SingleMatch( *mol ) ) {
        match_found = true;
        break;
      }
    }
    if( match_found ) {
      left_list.push_back( ones_to_do[i] );
    } else {
      right_list.push_back( ones_to_do[i] );
    }
  }
  QApplication::restoreOverrideCursor();

  left_panel_->add_data( left_list );
  left_panel_->set_subsearches( sub_searches );
  QString title = "Matched : " + list_name;
  left_panel_->set_title( title );

  if( show_non_matches ) {
    right_panel_->show();
    right_panel_->add_data( right_list );
    title = "Didn't Match : " + list_name;
    right_panel_->set_title( title );
  } else {
    right_panel_->hide();
  }

}

// ****************************************************************************
void SmiV::get_query_to_use( vector<char> &sel_query ,
                             const vector<pair<string,string> > &query_set ,
                             bool single_sel ) {

  sel_query = vector<char>( query_set.size() , 0 );
  vector<QString> query_names;
  typedef pair<string,string> QP;
  BOOST_FOREACH( const QP query , query_set ) {
    query_names.push_back( QString( query.first.c_str() ) );
  }

  QT4SelectItems si( "Which query to Use?" , query_names , sel_query , single_sel , this );
  si.exec();

}

// ****************************************************************************
void SmiV::build_sub_searches_from_smarts( const vector<char> &sel_smarts ,
                                           QString &smarts_list ,
                                           vector<pair<boost::shared_ptr<OESubSearch>,string> > &sub_searches ) {

  smarts_list = "";
  for( int i = 0 , is = sel_smarts.size() ; i < is ; ++i ) {
    if( !sel_smarts[i] ) {
      continue;
    }
    if( smarts_list.isEmpty() ) {
      smarts_list = smarts_[i].first.c_str();
    } else {
      smarts_list += QString( "|%1" ).arg( smarts_[i].first.c_str() );
    }

    OESubSearch *subs;
    try {
      subs = DACLIB::create_eosubsearch( smarts_[i].second , smarts_[i].first ,
                                         false , smarts_sub_defn_ );
    } catch( DACLIB::SMARTSSubDefnError &e ) {
      QMessageBox::warning( this , "SMARTS Error" , e.what() );
      continue;
    } catch( DACLIB::SMARTSDefnError &e ) {
      QMessageBox::warning( this , "SMARTS Error" , e.what() );
      continue;
    }
    sub_searches.push_back( make_pair( boost::shared_ptr<OESubSearch>( subs ) ,
                                       smarts_[i].first ) );
  }

}

// ****************************************************************************
void SmiV::build_sub_searches_from_mdl_queries( const vector<char> &sel_mdl_queries ,
                                                QString &mdl_list ,
                                                vector<pair<boost::shared_ptr<OESubSearch>,string> > &sub_searches ) {

  mdl_list = "";
  for( int i = 0 , is = sel_mdl_queries.size() ; i < is ; ++i ) {
    if( !sel_mdl_queries[i] ) {
      continue;
    }
    if( mdl_list.isEmpty() ) {
      mdl_list = mdl_queries_[i].first.c_str();
    } else {
      mdl_list += QString( "|%1" ).arg( mdl_queries_[i].first.c_str() );
    }

    oemolistream ims;
    ims.openstring( mdl_queries_[i].second );
    unsigned int aromodel = OEIFlavor::Generic::OEAroModelDaylight;
    unsigned int qflavor  = ims.GetFlavor( OEFormat::MDL );
    ims.SetFlavor( OEFormat::MDL , qflavor|aromodel );

    OEQMol qmol;
    OEReadMDLQueryFile( ims , qmol , OEMDLQueryOpts::Default );

    OESubSearch *subs = new OESubSearch( qmol );
    sub_searches.push_back( make_pair( boost::shared_ptr<OESubSearch>( subs ) ,
                                       mdl_queries_[i].first ) );
  }

}

// ****************************************************************************
void SmiV::add_smarts_definition( const QString &smarts_name , const QString &smarts_def ) {

  string smt_name( smarts_name.toLocal8Bit().data() );
  bool over_write = false;
  bool found_in_full = check_existing_smarts( smt_name , smarts_ , over_write );
  if( found_in_full ) {
    cout << "check_existing_smarts 1 true" << endl;
    if( !over_write ) {
      return;
    }
  }
  bool found_in_sub = true; // if it's found_in_full, it must be found_in_sub
  if( !found_in_full ) {
    found_in_sub = check_existing_smarts( smt_name , smarts_sub_defn_ , over_write );
  }

  if( (found_in_full || found_in_sub ) && !over_write ) {
    cout << "not over-writing" << endl;
    return;
  }

  if( !found_in_full && !found_in_sub ) {
    // it's all new, so put it in both. All smarts_ entries should also be in smarts_sub_defn_.
    smarts_.push_back( make_pair( smarts_name.toLocal8Bit().data() , smarts_def.toLocal8Bit().data() ) );
    smarts_sub_defn_.push_back( smarts_.back() );
    return;
  }

  vector<pair<string,string> >::iterator ssdp =
      find_if( smarts_sub_defn_.begin() , smarts_sub_defn_.end() ,
               bind( std::equal_to<string>() ,
                     bind( &pair<string,string>::first , _1 ) , smt_name ) );

  if( found_in_sub && !found_in_full ) {
    // if only in smarts_sub_defn, just replace what's there
    ssdp->second = smarts_def.toLocal8Bit().data();
  } else if( found_in_full && found_in_sub ) {
    // if in both, replace both
    vector<pair<string,string> >::iterator fsp =
        find_if( smarts_.begin() , smarts_.end() ,
                 bind( std::equal_to<string>() ,
                       bind( &pair<string,string>::first , _1 ) , smt_name ) );
    fsp->second = smarts_def.toLocal8Bit().data();
    ssdp->second = smarts_def.toLocal8Bit().data();
  }

}

// ****************************************************************************
void SmiV::add_smarts_definition( QTSmartsEditDialog &sed ) {

  QString smt_string = sed.get_smarts();
  QString smt_name = sed.get_smarts_name();
#ifdef NOTYET
  cout << "SMARTS : " << smt_name.toStdString() << " = "
      << smt_string.toStdString() << endl;
#endif
  add_smarts_definition( smt_name , smt_string );

}

// ****************************************************************************
void SmiV::store_next_mdl_query( const vector<string>  &next_query ,
                                 const QString &filename , int count ) {

  // the query needs a little TLC, as if it has R Groups in it, it may have
  // extraneous guff that OEChem doesn't seem to deal with
  // first, strip out everything that starts with a $. Then, as soon as
  // we hit M  END, stop as we don't want anything after that.
  string query_string;
  BOOST_FOREACH( string next_line , next_query ) {
    if( '$' != next_line[0] ) {
      query_string += next_line + '\n';
    }
    if( string( "M  END") == next_line.substr( 0 , 6 ) ) {
      cout << "breaking on M  END" << endl;
      break;
    }

  }

  QString short_name( filename.section( "/" , -1 ) );
  string query_name = string( short_name.toLocal8Bit().data() ) + "_" +
                      lexical_cast<string>( count );
  mdl_queries_.push_back( make_pair( query_name , query_string ) );

}

// ****************************************************************************
// sees if the SMARTS name is already used. If it does, see what the user
// wants to do about it - returns true if user wants to leave it as it is,
// false if it's to be overwritten or it wasn't found
bool SmiV::check_existing_smarts( const string &smarts_name ,
                                  const vector<pair<string,string> > &smarts_defs ,
                                  bool &over_write ) {

  vector<pair<string,string> >::const_iterator p =
      find_if( smarts_defs.begin() , smarts_defs.end() ,
               bind( std::equal_to<string>() ,
                     bind( &pair<string,string>::first , _1 ) , smarts_name ) );
  if( p == smarts_defs.end() ) {
    over_write = false;
    return false;
  }

  QString msg = QString( "SMARTS named %1 already exists.\n Overwrite?").arg( smarts_name.c_str() );
  over_write = bool( QMessageBox::Yes == QMessageBox::question( this , "Existing SMARTS" , msg ,
                                                               QMessageBox::Yes | QMessageBox::No ,
                                                               QMessageBox::Yes ) );
  return true;

}

// ****************************************************************************
string SmiV::get_active_smiles() const {

  string act_smi;
  SmiVPanel *panel = get_active_panel();
  pSmiVRec smiv_rec = panel->current_smiv_rec();
  if( smiv_rec ) {
    act_smi = smiv_rec->in_smi();
  }

  return act_smi;

}

// ****************************************************************************
SmiVPanel *SmiV::get_active_panel() const {

  SmiVPanel *panel = left_panel_;
  if( !right_panel_->isHidden() && right_panel_->is_selected() ) {
    panel = right_panel_;
  }

  return panel;

}

// ****************************************************************************
SmiVPanel *SmiV::get_inactive_panel() const {

  SmiVPanel *panel = get_active_panel();

  if( panel == left_panel_ ) {
    return right_panel_;
  } else {
    return left_panel_;
  }

}

// ****************************************************************************
void SmiV::update_smiv_recs( const string &new_smiles , const string &new_name ) {

  pSmiVRec new_rec( new SmiVRecord( new_smiles , new_name ) );
  vector<pSmiVRec>::iterator p = find_if( smiv_recs_.begin() , smiv_recs_.end() ,
                                          bind( equal_to<string>() ,
                                                bind( &SmiVRecord::smi_name , _1 ) ,
                                                new_rec->smi_name() ) );
  if( p != smiv_recs_.end() ) {
    bool ok = true;
    QString msg = QString( "%1 already used. Please supply a new one or Ok to over-write.").arg( new_name.c_str() );
    QString nm = QInputDialog::getText( this , "New Molecule Name" , msg , QLineEdit::Normal , new_name.c_str() , &ok );
    if( !ok ) {
      QString nm1 = QInputDialog::getText( this , "New Molecule Name" , "New name?" , QLineEdit::Normal , nm , &ok );
      if( !ok ) {
        return;
      }
      if( nm1 != nm )  {
        // make sure this new name's ok
        update_smiv_recs( new_smiles , nm1.toLocal8Bit().data() );
        return;
      }
    } else {
      string nms( nm.toLocal8Bit().data() );
      if( nms != new_name )  {
        // make sure this new name's ok
        update_smiv_recs( new_smiles , nms );
        return;
      } else {
        new_rec->set_smi_name( nms );
      }
    }
    *p = new_rec;
  } else {
    smiv_recs_.push_back( new_rec );
  }

  if( right_panel_->isHidden() ) {
    left_panel_->add_data( smiv_recs_);
  }

}

// ****************************************************************************
void SmiV::add_mol_list( const string &list_name ) {

  SmiVPanel *sp = get_active_panel();
  add_mol_list( list_name , sp->smiv_recs() );

}

// ****************************************************************************
void SmiV::add_mol_list( const string &list_name ,
                         const vector<pSmiVRec> &new_recs ) {

  rec_lists_.push_back( make_pair( list_name , new_recs ) );
  QAction *list_action = new QAction( list_name.c_str() , this );
  mol_lists_menu_->insertAction( mol_list_separator_ , list_action );
  connect( list_action , SIGNAL( triggered() ) , this , SLOT( slot_show_mol_list() ) );

}

// ****************************************************************************
void SmiV::new_mol_list( QString list_name ) {

  bool ok;
  list_name = QInputDialog::getText( this , "List Name" , "What name for list?" ,
                                             QLineEdit::Normal , list_name , &ok );
  if( !ok || list_name.isEmpty() ) {
    return;
  }

  string sln( list_name.toLocal8Bit().data() );
  vector<pair<string,vector<pSmiVRec> > >::iterator p =
      find_if( rec_lists_.begin() , rec_lists_.end() ,
               bind( equal_to<string>() ,
                     bind( &pair<string,vector<pSmiVRec> >::first , _1 ) ,
                     sln ) );
  if( p == rec_lists_.end() ) {
    add_mol_list( sln );
  } else {
    if( QMessageBox::Ok == QMessageBox::question( this , "List Name Already Used" ,
                                                  "That name is already in use. Over-write?" ) ) {
      SmiVPanel *sp = get_active_panel();
      p->second = sp->smiv_recs();
    }
  }

}

// ****************************************************************************
// for the special case when the data file that has been read into the table contained the columns
// CoreSmiles and CoreSmarts - extract the SMILES/SMARTS strings into the relevant lists
void SmiV::build_core_smiles_list() {

  int smiles_col = get_data_table_column_number( QString( "CoreSmiles" ) );
  if( -1 == smiles_col ) {
    // need an error message
    return;
  }

  vector<pair<string,vector<pSmiVRec> > >::iterator p =
      find_if( rec_lists_.begin() , rec_lists_.end() ,
               bind( equal_to<string>() ,
                     bind( &pair<string,vector<pSmiVRec> >::first , _1 ) ,
                     string( "CoreSmiles" ) ) );
  if( p == rec_lists_.end() ) {
    add_mol_list( "CoreSmiles" , vector<pSmiVRec>() );
    p = rec_lists_.begin() + ( rec_lists_.size() - 1 );
  }

  for( int i = 0 , is = data_table_->rowCount() ; i < is ; ++i ) {
    SmiVRecord *smiv_rec = new SmiVRecord( data_table_->data( i , smiles_col ).toString().toLocal8Bit().data() ,
                                           data_table_->data( i , 0 ).toString().toLocal8Bit().data() );
    smiv_recs_.push_back( pSmiVRec( smiv_rec ) );
    p->second.push_back( smiv_recs_.back() );
  }

}

// ****************************************************************************
void SmiV::build_core_smarts_list() {

  int smarts_col = get_data_table_column_number( QString( "CoreSmarts" ) );
  if( -1 == smarts_col ) {
    // need an error message
    return;
  }

  for( int i = 0 , is = data_table_->rowCount() ; i < is ; ++i ) {
    add_smarts_definition( data_table_->data( i , 0 ).toString() ,
                           data_table_->data( i , smarts_col ).toString() );
  }

}

// ****************************************************************************
void SmiV::show_core_smiles_molecule( int mol_num ) {

  show_mol_list( QString( "CoreSmiles" ) );
  string mol_name( data_table_->data( mol_num , 0 ).toString().toLocal8Bit().data() );
  left_panel_->go_to_first_mol();
  left_panel_->show_molecule( mol_name , 0 );

}

// ****************************************************************************
void SmiV::search_with_core_smarts( int smarts_num ) {

  string smarts_name( data_table_->data( smarts_num , 0 ).toString().toLocal8Bit().data() );
  vector<char> smarts_to_use( smarts_.size() , 0 );
  vector<pair<string,string> >::iterator p =
      find_if( smarts_.begin() , smarts_.end() ,
               bind( std::equal_to<string>() ,
                     bind( &pair<string,string>::first , _1 ) , smarts_name ) );
  if( p == smarts_.end() ) {
    return; // something's badly awry, however
  }

  smarts_to_use[distance( smarts_.begin() , p )] = 1;
  vector<pair<boost::shared_ptr<OESubSearch>,string> > sub_searches;
  QString smarts_list;
  build_sub_searches_from_smarts( smarts_to_use , smarts_list , sub_searches );
  show_all_molecules();
  do_substructure_matching( sub_searches , smarts_list , false );
  left_panel_->go_to_first_mol();

  new_mol_list( QString( smarts_name.c_str() ) );

}

// ****************************************************************************
// for all the CoreSmarts column in the data_table_, run against the full molecule list and
// calculate the number of R Group subst points round each core that are used, and the number
// of unique R Groups, feeding the results back into the columns RgroupPositions and
// NumberOfUniqueRgroups
void SmiV::do_rgroup_analysis_of_core_smarts() {

  vector<char> smarts_to_use( smarts_.size() , 1 );
  vector<pair<boost::shared_ptr<OESubSearch>,string> > sub_searches;
  QString smarts_list;
  build_sub_searches_from_smarts( smarts_to_use , smarts_list , sub_searches );
  vector<set<int> > rgroup_pos( smarts_.size() , set<int>() );
  vector<set<string> > unique_rgroups( smarts_.size() , set<string>() );
  vector<int> core_counts( smarts_.size() , 0 );

  BOOST_FOREACH( pSmiVRec rec , smiv_recs_ ) {
#ifdef NOTYET
      cout << "doing molecule " << rec->smi_name() << " : " << rec->in_smi() << endl;
#endif
      OEGraphMol mol;
      OEParseSmiles( mol , rec->in_smi() );
      bool core_mol( rec->smi_name().substr( 0 , 4 ) == string( "core" ) );
      for( int i = 0 , is = sub_searches.size() ; i < is ; ++i ) {
#ifdef NOTYET
        cout << "doing smarts " << i << " : " << sub_searches[i].second << " : " << smarts_[i].first << " : " << smarts_[i].second << endl;
#endif
        if( sub_searches[i].second.substr( 0 , 4 ) == string( "core" ) ) {
          if( core_mol ) {
            rgroup_core_count( mol , sub_searches[i].first , core_counts[i] );
          } else {
            rgroup_counts_and_strip( mol , sub_searches[i].first , rgroup_pos[i] , unique_rgroups[i] );
          }
        }
      }
  }

  update_rgroup_position_counts( rgroup_pos , unique_rgroups , core_counts );

}

// ****************************************************************************

void SmiV::rgroup_counts_and_strip( OEGraphMol &mol , boost::shared_ptr<OESubSearch> &sub ,
                                    set<int> &rgroup_pos , set<string> &unique_rgroups ) {

  // just want the first match, I think.  Symmetry issues aren't of interest.
  OEIter<OEMatchBase> match = sub->Match( mol , true );
  if( !match ) {
    return;
  }
  vector<unsigned int> in_core( DACLIB::max_atom_index( mol ) , 0 );
  for( OEIter<OEAtomBase> atom = match->GetTargetAtoms() ; atom ; ++atom ) {
    in_core[DACLIB::atom_index( *atom )] = 1;
  }
  int i = 0;
  for( OEIter<OEAtomBase> atom = match->GetTargetAtoms() ; atom ; ++atom , ++i ) {
    for( OEIter<OEAtomBase> conn = atom->GetAtoms() ; conn ; ++conn ) {
      if( !in_core[DACLIB::atom_index( *conn )] ) {
        rgroup_pos.insert( i );
        unique_rgroups.insert( generate_rgroup_smiles( mol , conn , in_core ) );
      }
    }
  }

}

// ****************************************************************************
string SmiV::generate_rgroup_smiles( OEGraphMol &mol , OEAtomBase &first_atom ,
                                     const vector<unsigned int> &in_core ) {

  string rgroup_smi;

  OEGraphMol mol_copy( mol );
  vector<unsigned int> ats_to_keep( DACLIB::max_atom_index( mol_copy ) , 0 );
  list<OEAtomBase *> to_do( 1 , &first_atom );
  vector<unsigned int> done_atoms( in_core );
  while( !to_do.empty() ) {
    OEAtomBase *next_at = to_do.front();
    done_atoms[DACLIB::atom_index( *next_at )] = 1;
    to_do.pop_front();
    ats_to_keep[DACLIB::atom_index( *next_at )] = 1;
    for( OEIter<OEAtomBase> conn_ats = next_at->GetAtoms() ; conn_ats ; ++conn_ats ) {
      if( !done_atoms[DACLIB::atom_index( *conn_ats )] ) {
        to_do.push_back( conn_ats );
      }
    }
  }

  for( OEIter<OEAtomBase> atom = mol_copy.GetAtoms() ; atom ; ++atom ) {
    // if the atom going is attached to a keep atom, then put something on to show this.
    // if the keep atom is the first atom, then put attach an Xe, otherwise attach
    // a Y.
    if( !ats_to_keep[DACLIB::atom_index( *atom )] ) {
      for( OEIter<OEAtomBase> conn_at = atom->GetAtoms() ; conn_at ; ++conn_at ) {
        if( ats_to_keep[DACLIB::atom_index( *conn_at )] ) {
          unsigned int tag_atom = OEElemNo::H;
          if( DACLIB::atom_index( *conn_at ) == DACLIB::atom_index( first_atom ) ) {
            tag_atom = OEElemNo::Xe;
          } else {
            tag_atom = OEElemNo::Y;
          }
          OEAtomBase *new_at = mol_copy.NewAtom( tag_atom );
          mol_copy.NewBond( new_at , conn_at , mol_copy.GetBond( atom , conn_at )->GetOrder() );
        }
      }
      mol_copy.DeleteAtom( atom );
    }
  }
  OECreateCanSmiString( rgroup_smi , mol_copy );

  return rgroup_smi;

}

// ****************************************************************************
void SmiV::update_rgroup_position_counts( const vector<set<int> > &rgroup_pos ,
                                          const vector<set<string> > &unique_rgroups ,
                                          const vector<int> &core_counts ) {

  int pos_count_col = get_data_table_column_number( "RgroupPositions" );
  if( -1 == pos_count_col ) {
    QMessageBox::warning( this , "No column" , "No column named RgroupPositions for R Group counts, so skipping." );
    return;
  }
  int uniq_count_col = get_data_table_column_number( "NumberOfUniqueRgroups" );
  if( -1 == uniq_count_col ) {
    QMessageBox::warning( this , "No column" , "No column named NumberOfUniqueRgroups for R Group counts, so skipping." );
    return;
  }
  int core_count_col = get_data_table_column_number( "CoreContainedIn" );
  if( -1 == core_count_col ) {
    QMessageBox::warning( this , "No column" , "No column named CoreContainedIn for R Group counts, so skipping." );
    return;
  }

  for( int i = 0 , is = rgroup_pos.size() ; i < is ; ++i ) {
    if( rgroup_pos[i].size() ) {
      data_table_->change_data( i , pos_count_col , QVariant( static_cast<int>( rgroup_pos[i].size() ) ) );
      data_table_->change_data( i , uniq_count_col , QVariant( static_cast<int>( unique_rgroups[i].size() ) ) );
      data_table_->change_data( i , core_count_col , QVariant( core_counts[i] ) );
    }
  }

}

// ****************************************************************************
// when molecule and subsearch are both cores, count whether mol contains sub
void SmiV::rgroup_core_count( OEChem::OEGraphMol &mol , boost::shared_ptr<OESubSearch> &sub ,
                              int &core_count ) {

  if( sub->SingleMatch( mol ) ) {
    ++core_count;
  }

}

// ****************************************************************************
// return the number of the named column, or -1 if it wasn't found
int SmiV::get_data_table_column_number( const QString &col_name ) const {

  for( int i = 0 , is = data_table_->columnCount() ; i < is ; ++i ) {
    if( col_name == data_table_->headerData( i , Qt::Horizontal ).toString() ) {
      return i;
    }
  }

  return -1;

}
