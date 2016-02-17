//
// file SmivDataTable.H
// David Cosgrove
// AstraZeneca
// 11th Feb 2011
//
// This is the declaration of the class SmivDataTable, derived from QAbstractTableModel.
// It holds the data read from a file which will be displayed in SmiV.

#ifndef SMIVDATATABLE_H
#define SMIVDATATABLE_H

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <QAbstractTableModel>

// *****************************************************************************************

class QString;

typedef boost::shared_ptr<QVariant> pQVar;

// *****************************************************************************************

class SmivDataTable : public QAbstractTableModel {

public :

  SmivDataTable( QObject *parent = 0 );
  int rowCount( const QModelIndex &parent = QModelIndex() ) const;
  int columnCount( const QModelIndex &parent = QModelIndex() ) const;
  QVariant data( const QModelIndex &index , int role = Qt::DisplayRole ) const;
  QVariant data( int row_num , int col_num ) const;
  QVariant headerData( int section , Qt::Orientation orientation ,
                       int role = Qt::DisplayRole ) const;
  Qt::ItemFlags flags( const QModelIndex &index ) const {
    return Qt::ItemIsSelectable | Qt::ItemIsEnabled;
  }

  void sort( int column , Qt::SortOrder order = Qt::AscendingOrder );

  void read_data_from_file( QWidget *parent_widget , const QString &filename );
  Qt::SortOrder last_sort_order() const { return last_sort_order_; }

  void change_data( int row_num , int col_num , const QVariant &new_val );

private :

  std::vector<std::vector<pQVar> > data_;
  std::vector<QString> col_names_;
  std::vector<int> sort_order_;
  Qt::SortOrder last_sort_order_;

  int get_sorted_row_number( int raw_row_num ) const;
  void sort_by_number( int col_num , Qt::SortOrder order );
  void sort_by_string( int col_num , Qt::SortOrder order );

};

#endif // SMIVDATATABLE_H
