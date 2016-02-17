//
// file SmiV.H
// David Cosgrove
// AstraZeneca
// 9th February 2009.

#ifndef DAC_SMIV
#define DAC_SMIV

#include <set>
#include <string>
#include <vector>

#include <QMainWindow>
#include <QModelIndex>
#include <QString>

#include <boost/shared_ptr.hpp>

// *************************************************************************

class SmivDataTable;
class SmiVFindMoleculeDialog;
class SmiVPanel;
class SmiVRecord;
class QTSmartsEditDialog; // one of mine, not Qt's

class QAction;
class QLabel;
class QLineEdit;
class QMenu;
class QSlider;
class QString;
class QTableView;

namespace OEChem {
  class OEAtomBase;
  class OEGraphMol;
  class OESubSearch;
}

namespace DACLIB {
  class QTSmilesEditDialog;
}

typedef boost::shared_ptr<SmiVRecord> pSmiVRec;

// *************************************************************************

class SmiV : public QMainWindow {

  Q_OBJECT

  signals :

  void atom_selected( unsigned int oe_ind );

public slots :

  void slot_new_smiles( QString new_smiles , QString new_name );

private slots :

  void slot_quit();
  void slot_file_read_mols();
  void slot_file_reread_mols();
  void slot_file_read_smarts();
  void slot_file_reread_smarts();
  void slot_file_read_mdl_query();
  void slot_file_read_data();
  void slot_write_smarts();
  void slot_clear_smarts();
  void slot_write_smiles();
  void slot_clear_molecules();
  void slot_find_mol();
  void slot_smarts_match();
  void slot_smarts_edit();
  void slot_smarts_input_int_pick();
  void slot_mdl_query_match();
  void slot_show_about_box();
  void slot_molecule_search_name( QString search_name , int search_mode );
  void slot_panel_selection_changed( QWidget *panel );
  void slot_input_smiles();
  void slot_edit_smiles();
  void slot_full_list();
  void slot_save_mol_list();
  void slot_new_mol_list();
  void slot_show_mol_list();
  void slot_sort_data_table( int col_num );
  void slot_data_table_cell_double_clicked( const QModelIndex &ind );
  void slot_data_table_show_row( QString row_name );

public :

  SmiV();
  ~SmiV();

  void parse_args( int argc , char **argv );

private :

  QAction *file_read_mol_ , *file_reread_mol_ , *file_read_smarts_ , *file_reread_smarts_;
  QAction *file_read_mdl_query_ , *file_read_data_;
  QAction *file_write_smiles_ , *file_quit_;
  QAction *find_mol_ , *input_smiles_ , *edit_smiles_ , *clear_mols_;
  QAction *full_list_ , *save_list_ , *new_list_ , *mol_list_separator_;
  QAction *smarts_match_ , *smarts_input_edit_ , *smarts_input_int_pick_ , *smarts_write_ , *clear_smarts_;
  QAction *mdl_query_match_;
  QAction *help_show_about_;
  QMenu *mol_lists_menu_;
  QLineEdit *in_smiles_;
  QTableView *data_table_view_;
  SmivDataTable *data_table_;

  SmiVPanel *left_panel_ , *right_panel_;
  SmiVFindMoleculeDialog *find_mol_dialog_;
  DACLIB::QTSmilesEditDialog *smiles_edit_dialog_;

  QString last_dir_ , last_mol_file_ , last_smarts_file_;

  // SMILES records
  std::vector<pSmiVRec> smiv_recs_;
  std::vector<std::pair<std::string,std::vector<pSmiVRec> > > rec_lists_;

  // SMARTS records
  std::vector<std::pair<std::string,std::string> > smarts_;
  std::vector<std::pair<std::string,std::string> > smarts_sub_defn_;

  // MDL Query file records - 1 string per query mol, file can contain multiple queries
  std::vector<std::pair<std::string,std::string> > mdl_queries_;

  std::string usage_text_;

  void build_actions();
  void build_file_actions();
  void build_smarts_actions();
  void build_mdl_query_actions();
  void build_molecule_actions();
  void build_help_actions();

  void build_menubar();
  void build_widget();

  void read_mol_file( const QString &filename );
  void read_smiles_file( const QString &filename );
  void read_other_mol_file( const QString &filename );
  void read_smarts_file( const QString &filename );
  void read_mdl_query_file( const QString &filename );
  void read_data_file( const QString &filename );
  void write_smiles_file( const QString &filename );
  void write_smarts_file( const QString &filename );

  void update_status_count();

  // reset to just left_panel_, showing all molecules
  void show_all_molecules();
  void show_mol_list( const QString &list_name );
  // get name of molecule in row_num of data table and show whichever panel currently has it
  void show_molecule_from_data_table( int row_num );
  // do smarts matching on contents of left panel only
  void do_smarts_matching();
  // do MDL query matching on contents of left panel only
  void do_mdl_query_matching();
  void do_substructure_matching( std::vector<std::pair<boost::shared_ptr<OEChem::OESubSearch>,std::string> > &sub_searches ,
                                 const QString &list_name , bool show_non_matches );

  void get_query_to_use( std::vector<char> &sel_smarts ,
                         const std::vector<std::pair<std::string,std::string> > &query_set ,
                         bool single_sel );
  void build_sub_searches_from_smarts( const std::vector<char> &sel_smarts ,
                                       QString &smarts_list ,
                                       std::vector<std::pair<boost::shared_ptr<OEChem::OESubSearch>,std::string> > &sub_searches );
  void build_sub_searches_from_mdl_queries( const std::vector<char> &sel_mdl_queries ,
                                            QString &mdl_list ,
                                            std::vector<std::pair<boost::shared_ptr<OEChem::OESubSearch>,std::string> > &sub_searches );
  void add_smarts_definition( const QString &smarts_name , const QString &smarts_def );
  void add_smarts_definition( QTSmartsEditDialog &sed );

  void store_next_mdl_query( const std::vector<std::string> &next_query ,
                             const QString &filename , int count );

  // sees if the SMARTS name is already used. If it does, see what the user
  // wants to do about it - returns over_write = true if user wants to replace it
  bool check_existing_smarts( const std::string &smarts_name ,
                              const std::vector<std::pair<std::string,std::string> > &smarts_defs ,
                              bool &over_write );

  std::string get_active_smiles() const;

  SmiVPanel *get_active_panel() const;
  SmiVPanel *get_inactive_panel() const;

  void update_smiv_recs( const std::string &new_smiles , const std::string &new_name );
  void add_mol_list( const std::string &list_name );
  void add_mol_list( const std::string &list_name ,
                     const std::vector<pSmiVRec> &new_recs );
  void new_mol_list( QString list_name );

  // for the special case when the data file that has been read into the table contained the columns
  // CoreSmiles and CoreSmarts - extract the SMILES/SMARTS strings into the relevant lists
  void build_core_smiles_list();
  void build_core_smarts_list();

  // taking stuff out of the data_table_ - the mol_num or smarts_num refers to the row in the table.
  void show_core_smiles_molecule( int mol_num );
  void search_with_core_smarts( int smarts_num );
  // for all the CoreSmarts column in the data_table_, run against the full molecule list and
  // calculate the number of R Group subst points round each core that are used, and the number
  // of unique R Groups, feeding the results back into the columns RgroupPositions and
  // NumberOfUniqueRgroups
  void do_rgroup_analysis_of_core_smarts();
  void rgroup_counts_and_strip( OEChem::OEGraphMol &mol , boost::shared_ptr<OEChem::OESubSearch> &sub ,
                                std::set<int> &rgroup_pos , std::set<std::string> &unique_rgroups );
  std::string generate_rgroup_smiles( OEChem::OEGraphMol &mol , OEChem::OEAtomBase &first_atom ,
                                      const std::vector<unsigned int> &in_core );
  void update_rgroup_position_counts( const std::vector<std::set<int> > &rgroup_pos ,
                                      const std::vector<std::set<std::string> > &unique_rgroups ,
                                      const std::vector<int> &core_counts );
  // when molecule and subsearch are both cores, count whether mol contains sub
  void rgroup_core_count( OEChem::OEGraphMol &mol , boost::shared_ptr<OEChem::OESubSearch> &sub ,
                          int &core_count );

  // return the number of the named column, or -1 if it wasn't found
  int get_data_table_column_number( const QString &col_name ) const;

};

#endif
