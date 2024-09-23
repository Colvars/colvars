// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_REPLICAS_H
#define COLVARPROXY_REPLICAS_H


/// \brief Methods for multiple-replica communication
class colvarproxy_replicas {

public:

  /// Constructor
  colvarproxy_replicas();

  /// Destructor
  virtual ~colvarproxy_replicas();

  /// Indicate if multi-replica support is available and active
  virtual int check_replicas_enabled();

  /// Index of this replica
  virtual int replica_index();

  /// Total number of replicas
  virtual int num_replicas();

  /// Synchronize replica with others
  virtual void replica_comm_barrier();

  /// Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);

  /// Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);
};

#endif
