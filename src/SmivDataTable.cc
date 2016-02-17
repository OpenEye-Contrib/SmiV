//
// file SmivDataTable.cc
// David Cosgrove
// AstraZeneca
// 11th Feb 2011
//

#include "SmivDataTable.H"

#include <QFileInfo>
#include <QMessageBox>
#include <QString>

#include <fstream>
#include <iostream>
#include <limits>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

using namespace boost;
using namespace std;

// *****************************************************************************************
SmivDataTable::SmivDataTable( QObject *parent ) :
    QAbstractTableModel( parent ) ,
    last_sort_order_( Qt::AscendingOrder ) {

}

// *****************************************************************************************
int SmivDataTable::rowCount( const QModelIndex &parent ) const {

  return data_.size();

}

// *****************************************************************************************
int SmivDataTable::columnCount( const QModelIndex &parent ) const {

  return col_names_.size();

}

// *****************************************************************************************
QVariant SmivDataTable::data( const QModelIndex &index , int role ) const {

  if( !index.isValid() || role != Qt::DisplayRole ) {
     return QVariant();
   }

  // convert the index number into the right index if the table is sorted,
  // otherwise returns index.row().
  int row_num = get_sorted_row_number( index.row() );

  if( index.column() < columnCount() && row_num < rowCount() ) {
    return *data_[row_num][index.column()];
  } else {
    return QVariant();
  }

  return QVariant();

}

// *****************************************************************************************
QVariant SmivDataTable::data( int row_num , int col_num ) const {

  // convert the index number into the right index if the table is sorted,
  // otherwise returns index.row().
  int sort_row_num = get_sorted_row_number( row_num );

  if( col_num < columnCount() && sort_row_num < rowCount() ) {
    return *data_[sort_row_num][col_num];
  } else {
    return QVariant();
  }

  return QVariant();

}

// *****************************************************************************************
QVariant SmivDataTable::headerData( int section , Qt::Orientation orientation ,
                                    int role ) const {

  if( role != Qt::DisplayRole ) {
    return QVariant();
  }

  if( orientation == Qt::Horizontal ) {
    if( section < columnCount() ) {
      return col_names_[section];
    } else {
      return QVariant();
    }
  } else {
    return QString( "Row %1" ).arg( section + 1 );
  }

}

// *****************************************************************************************
void SmivDataTable::sort( int column , Qt::SortOrder order ) {

  if( column >= columnCount() || column < 0 ) {
    return;
  }

  // see if it's numbers or strings
  bool is_num( true );
  for( int i = 0 , is = rowCount() ; i < is ; ++i ) {
    // Originally I did data_[i][column]->canConvert( Qt::Double ) but this always gave
    // true for a string, with a double value of 0.0. Made the whole canConvert business
    // seem a bit pointless, somehow!
    data_[i][column]->toString().toDouble( & is_num );
    if( !is_num ) {
      break;
    }
  }

  if( is_num ) {
    sort_by_number( column , order );
  } else {
    sort_by_string( column , order );
  }

  last_sort_order_ = order;

  emit layoutChanged(); // tell the view

}

// *****************************************************************************************
void SmivDataTable::read_data_from_file( QWidget *parent_widget , const QString &filename ) {

  ifstream ifs( filename.toLocal8Bit().data() );
  if( !ifs || !ifs.good() ) {
    QMessageBox::warning( parent_widget , "Data file error" ,
                          QString( "Couldn't open %1 for reading.").arg( filename ) );
    return;
  }

  string first_line;
  getline( ifs , first_line );
  vector<string> split_line;
  split( split_line , first_line , is_any_of( " ,\t") );
  beginInsertColumns( QModelIndex() , columnCount() , columnCount() + split_line.size() );
  BOOST_FOREACH( string &header , split_line ) {
    col_names_.push_back( header.c_str() );
  }
  endInsertColumns();

  string next_line;
  int line_count = 0;
  while( 1 ) {
    getline( ifs , next_line );
    if( ifs.eof() || !ifs.good() ) {
      break;
    }
    ++line_count;
    beginInsertRows( QModelIndex() , rowCount() , rowCount() );
    data_.push_back( vector<pQVar>() );
    split_line.clear();
    split( split_line , next_line , is_any_of( "," ) );
    // initialise CoreContainedIn column
    if( int( split_line.size() ) != columnCount() ) {
      QMessageBox msg( parent_widget );
      msg.setText( QString( "Line %1 has different number of columns from rest ").arg( line_count ) );
      msg.setInformativeText( QString( "First line had %1 columns, this has %2." ).arg( columnCount() ).arg( split_line.size() ) );
      msg.setStandardButtons( QMessageBox::Ignore | QMessageBox::Abort );
      msg.setDefaultButton( QMessageBox::Abort );
      int ret_val = msg.exec();
      if( QMessageBox::Abort == ret_val ) {
        break;
      } else if( QMessageBox::Ignore == ret_val ) {
        continue;
      }
    }
    BOOST_FOREACH( string &next_data , split_line ) {
      data_.back().push_back( pQVar( new QVariant( next_data.c_str() ) ) );
    }
    endInsertRows();

  }

  // make a new unsorted order
  sort_order_.clear();
  sort_order_.reserve( data_.size() );
  for( int i = 0 , is = data_.size() ; i < is ; ++i ) {
    sort_order_.push_back( i );
  }

}

// *****************************************************************************************
void SmivDataTable::change_data( int row_num , int col_num , const QVariant &new_val ) {

  if( row_num < rowCount() && col_num < columnCount() ) {
    data_[row_num][col_num] = pQVar( new QVariant( new_val ) );
  }
}

// *****************************************************************************************
int SmivDataTable::get_sorted_row_number( int raw_row_num ) const {

  if( raw_row_num < static_cast<int>( sort_order_.size() ) ) {
    return sort_order_[raw_row_num];
  } else {
    return numeric_limits<int>::max();
  }

}

// *****************************************************************************************
// assumes col_num and numerical type of column has already been verified
void SmivDataTable::sort_by_number( int col_num , Qt::SortOrder order ) {

  vector<pair<int,double> > vals;
  for( int i = 0 , is = rowCount() ; i < is ; ++i ) {
    vals.push_back( make_pair( i , data_[i][col_num]->toDouble() ) );
  }

  if( order == Qt::AscendingOrder ) {
    std::sort( vals.begin() , vals.end() , bind( less<double>() ,
                                                 bind( &pair<int,double>::second , _1 ) ,
                                                 bind( &pair<int,double>::second , _2 ) ) );
  } else {
    std::sort( vals.begin() , vals.end() , bind( greater<double>() ,
                                                 bind( &pair<int,double>::second , _1 ) ,
                                                 bind( &pair<int,double>::second , _2 ) ) );
  }

  sort_order_.clear();
  transform( vals.begin() , vals.end() , back_inserter( sort_order_ ) ,
             bind( &pair<int,double>::first , _1 ) );

}

// *****************************************************************************************
void SmivDataTable::sort_by_string( int col_num , Qt::SortOrder order ) {

  vector<pair<int,string> > vals;
  for( int i = 0 , is = rowCount() ; i < is ; ++i ) {
    vals.push_back( make_pair( i , data_[i][col_num]->toString().toStdString() ) );
  }

  if( order == Qt::AscendingOrder ) {
    std::sort( vals.begin() , vals.end() , bind( less<string>() ,
                                                 bind( &pair<int,string>::second , _1 ) ,
                                                 bind( &pair<int,string>::second , _2 ) ) );
  } else {
    std::sort( vals.begin() , vals.end() , bind( greater<string>() ,
                                                 bind( &pair<int,string>::second , _1 ) ,
                                                 bind( &pair<int,string>::second , _2 ) ) );
  }

  sort_order_.clear();
  transform( vals.begin() , vals.end() , back_inserter( sort_order_ ) ,
             bind( &pair<int,string>::first , _1 ) );

}
