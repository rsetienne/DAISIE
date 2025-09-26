#pragma once
#include <vector>

enum species_type {I, A, C}; // I, A, C
enum extinction_type {clado_extinct, immig_parent, none};

struct island_spec_row {
  
  double id = -1.0; // 1
  
  double parent = -1.0; // 2
  
  double colonisation_time = -1.0;  // 3
  
  species_type type_species;   // 4
  
  std::string anc_type;  // 5
  
  double extinction_time = -1.0;  // 6
  
  extinction_type ext_type;   // 7
  
  island_spec_row() {
    parent = -1;
    id = -1;
    colonisation_time = -1.0;
    extinction_time = -1.0;
    ext_type = none;
  }
  
  island_spec_row(double colonist, 
                  double timeval, 
                  species_type st) {
    parent = colonist;
    id = colonist;
    colonisation_time = timeval;
    extinction_time = -1;
    type_species = st;
    ext_type = none;
  }
  
  island_spec_row(int id_,
                  double parent_,
                  double timeval, species_type st) {
    parent = parent_;
    id = id_;
    colonisation_time = timeval;
    extinction_time = -1;
    type_species = st;
    ext_type = none;
  }
  
  island_spec_row(int id_,
                  double parent_,
                  double timeval,
                  species_type st,
                  std::string anc,
                  double ext_t) {
    parent = parent_;
    id = id_;
    colonisation_time = timeval;
    type_species = st;
    anc_type = anc;
    extinction_time = ext_t;
    ext_type = none;
  }   
};

struct island_spec {
  std::vector< island_spec_row > data_;
  
  bool empty() const {
    return data_.empty();
  }
  
  size_t size() const {
    return data_.size();
  }
  
  void clear() {
    data_.clear();
  }
  
  island_spec_row& operator[](size_t index) {
    return data_[index];
  }
  
  island_spec_row operator[](size_t index) const {
    return data_[index];
  }
  
  size_t find_species(const size_t colonist) {
    size_t index = 0;
    for (; index < data_.size(); ++index) {
      if (data_[index].id == colonist) {
        return index;
      }
    }
    return data_.size();
  }
  
  void remove(size_t index) {
    // data_.erase(data_.begin() + index);
    data_[index] = data_.back();
    data_.pop_back();
  }
  
  void anagenesis(size_t index, size_t new_id) {
    data_[index].id = new_id;
    data_[index].type_species = species_type::A;
    data_[index].ext_type = immig_parent;
  }
  
  void clado_genesis_c(size_t index,
                       size_t max_spec_id,
                       double t) {
    data_[index].type_species = species_type::C;
    data_[index].id = max_spec_id + 1;
    auto old_anc_type = data_[index].anc_type;
    old_anc_type.push_back('B');
    data_[index].anc_type.push_back('A');
    data_[index].ext_type = none;
    
    data_.emplace_back(max_spec_id + 2,
                       data_[index].parent,
                       data_[index].colonisation_time,
                       species_type::C,
                       old_anc_type,
                       t);
  }
  
  void clado_genesis_not_c(const size_t& index,
                           const size_t& max_spec_id,
                           const double& t) {
    data_[index].type_species = species_type::C;
    data_[index].id = static_cast<double>(max_spec_id + 1);
    data_[index].anc_type.clear();
    data_[index].anc_type.push_back('A');
    data_[index].extinction_time = data_[index].colonisation_time;
    data_[index].ext_type = none;
    
    data_.emplace_back(max_spec_id + 2,
                       data_[index].parent,
                       data_[index].colonisation_time,
                       species_type::C,
                       std::string({'B'}),
                       t);
  }
  
  
  void cladogenesis(size_t index, size_t *max_spec_id, double t) {
    if (data_[index].type_species == species_type::C) {
      clado_genesis_c(index, *max_spec_id, t);
    } else {
      clado_genesis_not_c(index, *max_spec_id, t);
    }
    *max_spec_id += 2;
  }
};
